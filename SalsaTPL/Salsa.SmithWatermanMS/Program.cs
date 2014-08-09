using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using Bio;
using Bio.Algorithms.Alignment;
using Bio.IO.FastA;
using Bio.IO.GenBank;
using Bio.SimilarityMatrices;
using Salsa.Core;
using Salsa.Core.Bio.Algorithms;
using Salsa.Core.Blas;
using Salsa.Core.Configuration;


#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using TDistance = System.Int16;
#else
using TDistance = System.Double;
#endif

namespace Salsa.SmithWatermanMS
{
    class Program
    {
        private static string _projectName;

        private static int _sequenceCount = 0;
        private static int _nodeCount = 1;
        private static int _processPerNodeCount = 1;
        private const int ThreadPerProcessCount = 1;

        private static long _myOverallTime = 0;
        private static long _myOverallComputeTime = 0;
        private static long _myOverallMPIScatterTime = 0;
        private static long _myOverallDiskIOTime = 0;

        private static long _totalPairCount = 0;

        private static string _jobId = string.Empty;
        private static string _taskId = string.Empty;

        private static string _fastaFile = string.Empty;
        private static string _distanceFile = string.Empty;
        private static string _indexFile = string.Empty;
        private static string _timingFile = string.Empty;
        private static string _summaryFile = string.Empty;
        private static string _configFile = string.Empty;
        private static bool _writeFullResults = true;

        private static bool _writeAlignments = false;
        private static string _writeAlignmentFile = string.Empty;
        private static StreamWriter _writeAlignmentsWriter = null;

        private static int _blockWriteFrequency = 0;
        private static string _blockDir;
        private static int _logWriteFrequency = 0;
        private static string _logDir;
        private static string _restartFile;
        private static bool _isRestart = false;
        private static int _lastSavedBlock = -1;

        private static int _gapOpen = 14;
        private static int _gapExtension = 4;
        private static string _scoringMatrixName = "EDNAFULL";
        private static MoleculeType _moleculeType = MoleculeType.NA;
        private static DistanceFunctionType _distanceFunction = DistanceFunctionType.PercentIdentity;
        private static string _emailResultsTo = "";

        private static IList<ISequence> _sequences;
        private static ConfigurationMgr _mgr;
        private static long _averageTime;
        private static long _averageComputeTime;
        private static long _averageScatterTime;
        private static long[] _computeTimes;

        private static bool _isOriginalRank0;
        private static int _rank;
        private static int _worldSize;
        private static string _doneStatusDir;


        static void Main(string[] args)
        {
            // Load the command line args into our helper class which allows us to name arguments
            Arguments pargs = new Arguments(args);
            pargs.Usage = Messages.InfoUsage;

            if (pargs.CheckRequired(new string[] {Messages.ParameterConfigFile, Messages.ParameterNodeCount}) == false)
            {
                Console.WriteLine(pargs.Usage);
                return;
            }


            /* Create MPI environment and dispose it after retrieving rank and world size */
            MPI.Environment env = new MPI.Environment(ref args);
            _rank = MPI.Communicator.world.Rank;
            _isOriginalRank0 = _rank == 0;
            _worldSize = MPI.Communicator.world.Size;

            try
            {
                if (env != null)
                {
                    env.Dispose();
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(string.Format(Messages.ErrorUnableToDisposeEnv, _rank, e.Message));
                throw;
            }

            /* Main logic starts here */
            try
            {
                Stopwatch overallTimer = Stopwatch.StartNew();

                /* Load configuration */
                _configFile = pargs.GetValue<string>(Messages.ParameterConfigFile);
                _mgr = ConfigurationMgr.LoadConfiguration(_configFile, true);
                ReadConfiguration(_mgr);

                _jobId = Environment.GetEnvironmentVariable(Messages.ParameterJobId);
                _taskId = Environment.GetEnvironmentVariable(Messages.ParameterTaskId);
                _nodeCount = pargs.GetValue<int>(Messages.ParameterNodeCount);
                _processPerNodeCount = _worldSize/_nodeCount;

                /* Fake rank and world size if restart */
                if (_isRestart)
                {
                    int tempOriginalRank = _rank;
                    ReadRestartFile();
                    Console.WriteLine(Messages.Restart, tempOriginalRank, _rank, _worldSize, _lastSavedBlock);
                }

                /* Read the fasta input file using FASTA parser */
                FastAParser parser = new FastAParser(_fastaFile);
                _sequences = parser.Parse().ToList();
                _sequenceCount = _sequences.Count;
                Console.WriteLine(Messages.InfoReadFasta, _rank, _sequenceCount, _fastaFile);

                /* Prepare write alignments */
                if (_writeAlignments)
                {
                    PrepareAlignmentWriter();
                }


                RowBlocksMetaData rowBlocksMetaData = new RowBlocksMetaData(_projectName, Environment.MachineName,
                                                                            _fastaFile, _blockWriteFrequency, _blockDir);
                MPIComputeLowerTriangularBlocked computeAlignments = new MPIComputeLowerTriangularBlocked(
                    _sequenceCount, ComputeBlock, rowBlocksMetaData, _rank, _worldSize);

                /* Block Compute Section */
                Stopwatch timer = Stopwatch.StartNew();
                Console.WriteLine(Messages.InfoBeginBlocksCompute, _rank, DateTime.Now);
                computeAlignments.Compute(_lastSavedBlock);
                Console.WriteLine(Messages.InfoEndBlocksCompute, _rank, DateTime.Now);
                timer.Stop();
                // Includes disk I/O time for partial files, status files, and alignment files
                _myOverallComputeTime = timer.ElapsedMilliseconds;

                /* Finalize alignment write */
                if (_writeAlignments)
                {
                    FinalizeAlignmentWrite();
                }

                overallTimer.Stop();
                _myOverallTime = overallTimer.ElapsedMilliseconds;

                if (_isOriginalRank0)
                {
                    WriteIndexFile(_indexFile, _sequences);
                    //WriteTimingFile(_timingFile);
                    //WriteSummaryFile(_summaryFile);
                    WriteConfiguration(_mgr);
                    _mgr.SaveAs(_configFile);
                }

                // Those who come here must write the message I'm done to see, which ones are completed.
                string imdonefile = Path.Combine(Path.GetDirectoryName(_indexFile), _doneStatusDir, (_rank + ".txt"));
                try
                {
                    using (StreamWriter writer = new StreamWriter(imdonefile))
                    {
                        string txt = string.Format(Messages.DoneMsg, _rank);
                        writer.WriteLine(txt);
                        Console.WriteLine(txt);
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(Messages.ErrorUnableToWriteIamDone, _rank, e.Message);
                }

            }
            catch (Exception e)
            {
                Console.WriteLine(Messages.ErrorGeneral, _rank, e.Message);
            }
        }

        private static void FinalizeAlignmentWrite()
        {
            if(_writeAlignmentsWriter != null)
            {
                _writeAlignmentsWriter.Flush();
                _writeAlignmentsWriter.Dispose();
            }
        }

        private static void PrepareAlignmentWriter()
        {
            if (!string.IsNullOrWhiteSpace(_writeAlignmentFile))
            {
                string fileName = string.Format("{0}_{1}{2}", Path.GetFileNameWithoutExtension(_writeAlignmentFile),
                                                _rank,
                                                Path.GetExtension(_writeAlignmentFile));
                string filePath = Path.Combine(Path.GetDirectoryName(_writeAlignmentFile), fileName);
                try
                {
                    _writeAlignmentsWriter = File.CreateText(filePath);
                }
                catch(Exception e)
                {
                    Console.WriteLine(Messages.ErrorCreatingFile, filePath, e.Message);
                    throw;
                }
            }
            else
            {
                Console.WriteLine(Messages.WarningIgnoreAlignmentWrite, _writeAlignmentFile);
                _writeAlignments = false;
            }
            
        }

        private static void ReadRestartFile()
        {
            /* Assumes the number of MPI processes is equal to the number of 
             * to be completed ranks */
            using (StreamReader reader = new StreamReader(_restartFile))
            {
                _worldSize = int.Parse(reader.ReadLine().Trim());
                // Skip lines till you find your rank
                for (int i = 0; i < _rank; i++)
                {
                    reader.ReadLine();
                }

                string [] splits = reader.ReadLine().Trim().Split(new[]{'\t',' '});
                _rank = int.Parse(splits[0]);
                _lastSavedBlock = int.Parse(splits[1]);
            }
        }


        static int ComputeBlock(Block block, PartialMatrix<TDistance> matrix, bool isDiagonal, int rank)
        {
            block.IsDiagonal = isDiagonal;

            // actually computed number of pairs
            int count = 0;

            string logFile = Path.Combine(_logDir, "log_rank_" + rank + ".txt");
            string nl = Environment.NewLine;
            
            SmithWatermanAligner aligner = new SmithWatermanAligner();
            SimilarityMatrix scoringMatrix = _mgr.SmithWatermanMS.LoadSimilarityMatrix(_scoringMatrixName, _moleculeType);

            int totalPairsToCompute;
            int computedPairCount = 0; // in the case of diagonal block actually computed number of pairs is half of this
            int logCounter = 0;

            Range range = block.RowRange;
            range.StartSeqName = _sequences[range.StartIndex].ID;
            range.EndSeqName = _sequences[range.EndIndex].ID;

            range = block.ColumnRange;
            range.StartSeqName = _sequences[range.StartIndex].ID;
            range.EndSeqName = _sequences[range.EndIndex].ID;

            ISequence si, sj;
            // Moving "if" inside loop is expensive, hence the justification for duplicate code.
            // Yet, the logWriting has to be done inside.
            /*if (isDiagonal)
            {
                totalPairsToCompute = (block.RowRange.Length*(block.ColumnRange.Length - 1))/2;

                for (int i = block.RowRange.StartIndex; i <= block.RowRange.EndIndex; i++)
                {
                    si = _sequences[i]; // ith row seqeunces in blcok
                    for (int j = block.ColumnRange.StartIndex; j < i; j++)
                    {
                        sj = _sequences[j]; // jth column seqeunce in block

                        // Do log writing if _logWriteFrequence > 0 && logCounter == _logWriteFrequency
                        if (_logWriteFrequency > 0 && logCounter == _logWriteFrequency)
                        {
                            WriteLog(rank, block, i, j, logFile, nl, computedPairCount, totalPairsToCompute, si.ID, sj.ID, false);
                            logCounter = 0;
                        }

                        IList<IPairwiseSequenceAlignment> psas = aligner.Align(scoringMatrix, _gapOpen, _gapExtension, si, sj);
                        IPairwiseSequenceAlignment psa = psas[0]; // Take the first alignment
                        IList<PairwiseAlignedSequence> pass = psa.PairwiseAlignedSequences;
                        PairwiseAlignedSequence pas = pass[0]; // Take the first PairwisedAlignedSequence

                        if ((_writeAlignments) && (_writeAlignmentsWriter != null))
                        {
//                            WriteAlignment(_writeAlignmentsWriter, pas, si, sj);
                            WriteAlignment(_writeAlignmentsWriter, pas, i, j, si, sj, scoringMatrix);
                        }

                        matrix[i, j] = matrix[j, i] = ComputeDistance(pas, scoringMatrix, si, sj);

                        computedPairCount += 2;
                        logCounter++;
                        count++;
                    }
                }
                // Write the done status
                 if (_logWriteFrequency > 0)
                 {
                     int i = block.RowRange.EndIndex;
                     int j = block.RowRange.EndIndex - 1;
                     WriteLog(rank, block, i, j, logFile, nl, computedPairCount, 
                         totalPairsToCompute, _sequences[i].ID, _sequences[j].ID, true);
                     
                 }
            }
            else
            {*/
                totalPairsToCompute = block.RowRange.Length * block.ColumnRange.Length;
                
                for (int i = block.RowRange.StartIndex; i <= block.RowRange.EndIndex; i++)
                {
                    si = _sequences[i]; // ith row seqeunces in blcok

                    for (int j = block.ColumnRange.StartIndex; j <= block.ColumnRange.EndIndex; j++)
                    {
                        sj = _sequences[j]; // jth column seqeunce in block

                        // Do log writing if _logWriteFrequence > 0 && logCounter == _logWriteFrequency
                        if (_logWriteFrequency > 0 && logCounter == _logWriteFrequency)
                        {
                            WriteLog(rank, block, i, j, logFile, nl, computedPairCount, totalPairsToCompute, si.ID, sj.ID, false);
                            logCounter = 0;
                        }

                        IList<IPairwiseSequenceAlignment> psas = aligner.Align(scoringMatrix, _gapOpen, _gapExtension, si, sj);
                        IPairwiseSequenceAlignment psa = psas[0]; // Take the first alignment
                        IList<PairwiseAlignedSequence> pass = psa.PairwiseAlignedSequences;
                        PairwiseAlignedSequence pas = pass[0]; // Take the first PairwisedAlignedSequence

                        if ((_writeAlignments) && (_writeAlignmentsWriter != null))
                        {
//                            WriteAlignment(_writeAlignmentsWriter, pas, si, sj);
                            WriteAlignment(_writeAlignmentsWriter, pas, i, j, si, sj, scoringMatrix);
                        }

                        matrix[i, j] = ComputeDistance(pas, scoringMatrix, si, sj);

                        computedPairCount++;
                        logCounter++;
                        count++;
                    }
                }
                // Write the done status
                 if (_logWriteFrequency > 0)
                 {
                     int i = block.RowRange.EndIndex;
                     int j = block.ColumnRange.EndIndex;
                     WriteLog(rank, block, i, j, logFile, nl, computedPairCount, 
                         totalPairsToCompute, _sequences[i].ID, _sequences[j].ID, true);
                     
                 }
            /*}*/
            return count;
        }

        private static void WriteLog(int rank, Block block, int i, int j, string logFile, string nl, 
            int computedPairCount, int totalPairsToCompute, string siID, string sjID, bool isComplete)
        {
            try
            {
                using (StreamWriter writer = new StreamWriter(File.Create(logFile)))
                {
                    string log = "MPIRank: " + rank + nl + nl +
                                 "RowBlock#:     " + block.RowBlockNumber + nl +
                                 "ColBlock#:     " + block.ColumnBlockNumber + nl + nl +
                                 "RowRange:      " + block.RowRange + nl +
                                 "ColRange:      " + block.ColumnRange + nl + nl +
                                 "RowSequence:   " + i + "  " + siID + nl +
                                 "ColSequence:   " + j + "  " + sjID + nl + nl +
                                 "ComputedPairs: " + computedPairCount + nl +
                                 "TotalPairs:    " + totalPairsToCompute + nl +
                                 "Percentage:    " + ((computedPairCount*100)/totalPairsToCompute) + nl +
                                 "Completed:     " + isComplete;
                    writer.WriteLine(log);
                }
            }
            catch(Exception e)
            {
                Console.WriteLine("Rank " + rank + " ignoring log write because: " + e.Message);
                return;
            }
        }

#if USE_UINT16 || USE_INT16
        static TDistance ComputeDistance(PairwiseAlignedSequence pas, SimilarityMatrix scoringMatrix, ISequence si, ISequence sj)
        {
            TDistance result = 0;

            if (_moleculeType != MoleculeType.Protein && _moleculeType != MoleculeType.Invalid)
            {
                switch (_distanceFunction)
                {
                    case DistanceFunctionType.PercentIdentity:
                        result = (TDistance)((1.0f - ComputePercentIdentity(pas)) * TDistance.MaxValue);
                        break;
                    case DistanceFunctionType.Kimura2:
                        result = (TDistance)(ComputeKimuraDistance(pas) * TDistance.MaxValue);
                        break;
                    case DistanceFunctionType.JukesCantor:
                        result = (TDistance)(ComputeJukesCantorDistance(pas) * TDistance.MaxValue);
                        break;
                    case DistanceFunctionType.MinMaxNormScore:
                        result = (TDistance)((1.0 - ComputeMinMaxNormScore(pas, si, sj, scoringMatrix)) * TDistance.MaxValue);
                        break;
                }
            }
            else
            {
                result = (TDistance)((1.0f - ComputePercentSimilarity(pas, scoringMatrix)) * TDistance.MaxValue);
            }

            return result;
        }

#else
        static TDistance ComputeDistance(PairwiseAlignedSequence pas, SimilarityMatrix scoringMatrix)
        {
            TDistance result = 0;

            if (_moleculeType != MoleculeType.Protein && _moleculeType != MoleculeType.Invalid)
            {
                switch (_distanceFunction)
                {
                    case DistanceFunctionType.PercentIdentity:
                        result = (TDistance)(1.0f - ComputePercentIdentity(pas));
                        break;
                    case DistanceFunctionType.Kimura2:
                        result = (TDistance)(ComputeKimuraDistance(pas));
                        break;
                    case DistanceFunctionType.JukesCantor:
                        result = (TDistance)(ComputeJukesCantorDistance(pas));
                        break;
                }
            }
            else
            {
                result = (TDistance)(1.0f - ComputePercentSimilarity(pas, scoringMatrix));
            }

            return result;
        }
#endif

        #region Distance Calculations

        static float ComputePercentIdentity(PairwiseAlignedSequence pas)
        {
            ISequence alignedSeqA = pas.FirstSequence;
            ISequence alignedSeqB = pas.SecondSequence;
            float identity = 0.0f;
            for (int i = 0; i < alignedSeqA.Count; i++)
            {
                char ca = Char.ToUpper((char)alignedSeqA[i]);
                char cb = Char.ToUpper((char)alignedSeqB[i]);
                if (ca == cb)
                {
                    identity++;
                }

            }
            return identity / alignedSeqA.Count;
        }

        private static double ComputePercentIdentityNonGap(PairwiseAlignedSequence pas)
        {
            byte gap = (byte) '-';
            ISequence alignedSeqA = pas.FirstSequence;
            ISequence alignedSeqB = pas.SecondSequence;
            int count = 0;
            float identity = 0.0f;
            for (int i = 0; i < alignedSeqA.Count; i++)
            {
                char ca = Char.ToUpper((char)alignedSeqA[i]);
                char cb = Char.ToUpper((char)alignedSeqB[i]);
                if (ca == cb) // ca and cb both will not be gaps
                {
                    ++identity;
                }
                if (!((((byte)ca) == gap) || (((byte)cb) == gap)))
                {
                    ++count;
                }

            }
            return identity / count;
        }

        static float ComputePercentSimilarity(PairwiseAlignedSequence pas, SimilarityMatrix scoringMatrix)
        {
            ISequence alignedSeqA = pas.FirstSequence;
            ISequence alignedSeqB = pas.SecondSequence;
            float similarity = 0.0f;
            for (int i = 0; i < alignedSeqA.Count; i++)
            {
                byte itemA = (byte)char.ToUpper((char)alignedSeqA[i]);
                byte itemB = (byte)char.ToUpper((char)alignedSeqB[i]);
                // Similarity pairs of residues are those with positive score in the scroing matrix
                // See Chapter 2 in Biological Sequence Analysis; Durbin, Eddy, Krogh and Mitchison; Cambridge Press; 1998. 

                IAlphabet alphabet = alignedSeqA.Alphabet;
                if (!alphabet.CheckIsGap(itemA) && !alphabet.CheckIsGap(itemB) && scoringMatrix[itemA, itemB] > 0)
                {
                    similarity++;
                }

            }
            return similarity / alignedSeqA.Count;
        }

        /* Obviously NOT GOOD */
        private static double ComputeMinMaxNormScore(PairwiseAlignedSequence pas, ISequence si, ISequence sj, SimilarityMatrix mat)
        {
            long s = pas.Score;
            int scoreSi = si.Sum(b => mat[b, b]);
            int scoreSj = sj.Sum(b => mat[b, b]);
            int maxS = Math.Max(scoreSi, scoreSj);
            // For SWG min score is zero
            const int minS = 0;
            return (s - minS) * 1.0 / (maxS - minS);
        }
        
        static double ComputeKimuraDistance(PairwiseAlignedSequence pas)
        {
            double length = 0;
            double gapCount = 0;
            double transitionCount = 0;    // P = A -> G | G -> A | C -> T | T -> C
            double transversionCount = 0;  // Q = A -> C | A -> T | C -> A | C -> G | T -> A  | T -> G | G -> T | G -> C

            ISequence alignedSeqA = pas.FirstSequence;
            ISequence alignedSeqB = pas.SecondSequence;

            byte itemA;
            byte itemB;

            IAlphabet alphabet = alignedSeqA.Alphabet;

            for (int i = 0; i < alignedSeqA.Count; i++)
            {
                length++;
                itemA = alignedSeqA[i]; //nucleotide 1
                itemB = alignedSeqB[i]; //nucleotide 2

                if (itemA != itemB)
                {
                    // Don't consider gaps at all in this computation;
                    if (alphabet.CheckIsGap(itemA) || alphabet.CheckIsGap(itemB))
                    {
                        gapCount++;
                    }
                    else if ((((char)itemA) == 'A' && ((char)itemB) == 'G') ||
                        (((char)itemA) == 'G' && ((char)itemB) == 'A') ||
                        (((char)itemA) == 'C' && ((char)itemB) == 'T') ||
                        (((char)itemA) == 'T' && ((char)itemB) == 'C'))
                    {
                        transitionCount++;
                    }
                    else
                    {
                        transversionCount++;
                    }
                }
            }

            double P = transitionCount / (length - gapCount);
            double Q = transversionCount / (length - gapCount);

            double artificialDistance = 10;
            if (1.0 - (2.0 * P + Q) <= double.Epsilon || 1.0 - (2.0 * Q) <= double.Epsilon)
            {
                PrintArtificialDistanceAlignments(alignedSeqA, alignedSeqB, artificialDistance, "Kimura2");
                return artificialDistance;
            }

            double d = (-0.5 * Math.Log(1.0 - 2.0 * P - Q) - 0.25 * Math.Log(1.0 - 2.0 * Q));
            if (d > artificialDistance)
            {
                PrintArtificialDistanceAlignments(alignedSeqA, alignedSeqB, artificialDistance, "Kimura2");
                return artificialDistance;
            }
            return d;
        }

        static void PrintArtificialDistanceAlignments(ISequence alignedSeqA, ISequence alignedSeqB, double aDistance, string distanceType)
        {
            Console.WriteLine("*******************************************");
            Console.WriteLine(alignedSeqA.ToString());
            Console.WriteLine(alignedSeqB.ToString());
            Console.WriteLine("Artificial " + distanceType + "Distance: " + aDistance);
            Console.WriteLine("*******************************************");
        }

        static double ComputeJukesCantorDistance(PairwiseAlignedSequence pas)
        {
            double length = 0;
            double gapCount = 0;
            double differenceCount = 0;

            ISequence alignedSeqA = pas.FirstSequence;
            ISequence alignedSeqB = pas.SecondSequence;

            byte itemA;
            byte itemB;
            IAlphabet alphabet = alignedSeqA.Alphabet;
            for (int i = 0; i < alignedSeqA.Count; i++)
            {
                length++;
                itemA = alignedSeqA[i]; //nucleotide 1
                itemB = alignedSeqB[i]; //nucleotide 2

                if (itemA != itemB)
                {
                    // Don't consider gaps at all in this computation;
                    if (alphabet.CheckIsGap(itemA) || alphabet.CheckIsGap(itemB))
                    {
                        gapCount++;
                    }
                    else
                    {
                        differenceCount++;
                    }
                }
            }

            double d = 1.0 - ((4.0 / 3.0) * (differenceCount / (length - gapCount)));

            // Todo: saliya - add an artificial value here
            double artificialDistance = 10;
            if (d <= double.Epsilon)
            {
                PrintArtificialDistanceAlignments(alignedSeqA, alignedSeqB, artificialDistance, "JukesCantor");
                return artificialDistance;
            }

            return (-0.75 * Math.Log(d));
        }

        #endregion

        #region Write Outputs (Alignment, Index, Timing, Summary, email)

        static void WriteAlignment(StreamWriter writer, PairwiseAlignedSequence pas, ISequence si, ISequence sj)
        {
            ISequence shortSeq = si.Count <= sj.Count ? si : sj;
            ISequence otherSeq = shortSeq.Equals(si) ? sj : si;

            writer.WriteLine(shortSeq.ID + "\t" + otherSeq.ID + "\t" + shortSeq.Count + "\t" + pas.FirstSequence.Count + "\t" + pas.SecondSequence.Count);
        }

        static void WriteAlignment(StreamWriter writer, PairwiseAlignedSequence pas, int globalRow, int globalCol, ISequence origFirstSeq, ISequence origSecondSeq, SimilarityMatrix mat)
        {
            double pid = ComputePercentIdentity(pas);
            double pidng = ComputePercentIdentityNonGap(pas);
            ISequence alignedFirstSeq = pas.FirstSequence;
            ISequence alignedSecondSeq = pas.SecondSequence;

            writer.WriteLine("#Row:{0}\tCol:{1}", globalRow, globalCol);
            //FirstLength\tSecondLength\tAlignmentLength
            writer.WriteLine("#FL:{0}\tSL:{1}\tAL:{2}", origFirstSeq.Count, origSecondSeq.Count, alignedFirstSeq.Count);
            //PID\tScore\tPIDNG\tFirstSelfLocalScore\tSecondSelfLocalScore\tFirstSelfGlobalScore\tSecondSelfGlobalScore
            writer.WriteLine("#PID:{0}\tScore:{1}\tPIDNG:{2}\tFSLS:{3}\tSSLS:{4}\tFSGS:{5}\tSSGS:{6}",
                             pid.ToString("F5"), pas.Score.ToString("F5"), pidng.ToString("F5"),
                             SelfAlignedScore(alignedFirstSeq, mat), SelfAlignedScore(alignedSecondSeq, mat),
                             SelfAlignedScore(origFirstSeq, mat), SelfAlignedScore(origSecondSeq, mat));
            writer.WriteLine(">" + alignedFirstSeq.ID);
            writer.WriteLine(SeqToString(alignedFirstSeq));
            writer.WriteLine(">" + alignedSecondSeq.ID);
            writer.WriteLine(SeqToString(alignedSecondSeq));
            writer.WriteLine();
        }

        static int SelfAlignedScore(IEnumerable<byte> seq, SimilarityMatrix mat)
        {
            const char gap = '-';
            return seq.Select(b => Char.ToUpper((char)b)).Where(c => c != gap).Sum(c => mat[(byte)c, (byte)c]);
        }

        static string SeqToString(IEnumerable<byte> seq)
        {
            ASCIIEncoding encoding = new ASCIIEncoding();
            return encoding.GetString(seq.ToArray()).ToUpper();
        }

        static void WriteIndexFile(string fileName, IEnumerable<Bio.ISequence> sequences)
        {
            using (StreamWriter writer = new StreamWriter(fileName))
            {
                int index = 0;

                foreach (Bio.ISequence sequence in sequences)
                {
                    writer.WriteLine("{0}\t{1}\t{2}", index, sequence.ID, sequence.Count);
                    index++;
                }

                writer.Flush();
                writer.Close();
            }
        }

        static void WriteSummaryFile(string fileName)
        {
            StringBuilder sb = new StringBuilder();
            string pattern = string.Format("{0}x{1}x{2}", ThreadPerProcessCount, _processPerNodeCount, _nodeCount);
            int parallelism = _nodeCount * _processPerNodeCount * ThreadPerProcessCount;
            sb.AppendFormat("Run Date: {0}", DateTime.Now.ToLongDateString());
            sb.AppendLine();
            sb.AppendFormat("Machine: {0}", System.Environment.MachineName);
            sb.AppendLine();
            sb.AppendFormat("Job Id: {0}", _jobId);
            sb.AppendLine();
            sb.AppendFormat("Task Id: {0}", _taskId);
            sb.AppendLine();
            sb.AppendFormat("Average Time: {0} ms, {1}", _averageTime, TimeSpan.FromMilliseconds(_averageTime).ToString());
            sb.AppendLine();
            sb.AppendFormat("Average Compute Time: {0} ms, {1}", _averageComputeTime, TimeSpan.FromMilliseconds(_averageComputeTime).ToString());
            sb.AppendLine();
            
            
//            sb.AppendFormat("Average MPI Scatter Time: {0} ms, {1}", _averageScatterTime, TimeSpan.FromMilliseconds(_averageScatterTime).ToString());
//            sb.AppendLine();
//            sb.AppendFormat("Rank0 Disk IO Time: {0} ms, {1}", _myOverallDiskIOTime, TimeSpan.FromMilliseconds(_myOverallDiskIOTime).ToString());
//            sb.AppendLine();
            
            
            sb.AppendFormat("Sequence Count: {0}", _sequenceCount);
            sb.AppendLine();
            sb.AppendFormat("Pairwise Alignments Count: {0}", _totalPairCount);
            sb.AppendLine();
            sb.AppendFormat("Pattern: {0}", pattern);
            sb.AppendLine();
            sb.AppendFormat("Threads Per Process Count: {0}", ThreadPerProcessCount);
            sb.AppendLine();
            sb.AppendFormat("Processes Per Node Count: {0}", _processPerNodeCount);
            sb.AppendLine();
            sb.AppendFormat("Node count: {0}", _nodeCount);
            sb.AppendLine();
            sb.AppendFormat("Parallelism: {0}", parallelism);
            sb.AppendLine();
            sb.AppendFormat("Gap Open Penalty: {0}", _gapOpen);
            sb.AppendLine();
            sb.AppendFormat("Gap Extension Penalty: {0}", _gapExtension);
            sb.AppendLine();
            sb.AppendFormat("Index File: {0}", _indexFile);
            sb.AppendLine();
            sb.AppendFormat("MoleculeType: {0}", _moleculeType);
            sb.AppendLine();
            sb.AppendFormat("DistanceFunction: {0}", _distanceFunction);
            sb.AppendLine();
            sb.AppendFormat("ScoringMatrixName: {0}", _scoringMatrixName);
            sb.AppendLine();
            sb.AppendFormat("DistanceMatrixFormat: {0}", typeof(TDistance));

            Console.WriteLine(sb.ToString());
            File.WriteAllText(fileName, sb.ToString());

            if (string.IsNullOrWhiteSpace(_emailResultsTo) == false)
            {
                EmailResults(_emailResultsTo, _jobId, _taskId, pattern, sb.ToString());
            }

        }

        static void WriteTimingFile(string fileName)
        {
            string pattern = string.Format("{0}x{1}x{2}", ThreadPerProcessCount, _processPerNodeCount, _nodeCount);
            int parallelism = _nodeCount * _processPerNodeCount * ThreadPerProcessCount;

            StringBuilder sb = new StringBuilder();
            sb.Append("Job Id: ");
            sb.AppendLine(_jobId);
            sb.Append("TaskId: ");
            sb.AppendLine(_taskId);
            sb.Append("Pattern (ThreadsPerProc x ProcPerNode x NodeCount): ");
            sb.AppendLine(pattern);
            sb.Append("Parallelism: ");
            sb.AppendLine(parallelism.ToString());
            sb.Append("Pairs Aligned: ");
            sb.AppendLine(_totalPairCount.ToString());
            sb.Append("Average Time (ms): ");
            sb.AppendLine(_averageTime.ToString());
            sb.Append("Average Compute Time (ms): ");
            sb.AppendLine(_averageComputeTime.ToString());


//            sb.Append("Average Scatter Time (ms)");
//            sb.AppendLine(_averageScatterTime.ToString());
//            sb.Append("Rank0 Disk I/O Time (ms) ");
//            sb.AppendLine(_myOverallDiskIOTime.ToString());
            
            
            sb.Append("Threads: ");
            sb.AppendLine(ThreadPerProcessCount.ToString());
            sb.Append("Processes per Node: ");
            sb.AppendLine(_processPerNodeCount.ToString());
            sb.Append("Node Count: ");
            sb.AppendLine(_nodeCount.ToString());
            sb.Append("Machine Name: ");
            sb.AppendLine(Environment.MachineName);
            sb.Append("Run Date: ");
            sb.AppendLine(DateTime.Now.ToString());

            long sumOfComputeTimes = _computeTimes.Sum();
            for (int i = 0; i < _computeTimes.Length; i++)
            {
                sb.Append("Rank ");
                sb.Append(i);
                sb.Append(" fraction of compute time: ");
                sb.Append((100*_computeTimes[i]/sumOfComputeTimes));
                sb.AppendLine("%");
            }

            using (StreamWriter writer = new StreamWriter(fileName))
            {
                writer.WriteLine(sb.ToString());
            }
        }

        static void EmailResults(string to, string jobId, string taskId, string pattern, string results)
        {
            try
            {
                string subjectLine = string.Format("SmithWatermanMS Timing: Machine {0} Job {1} Task {2} Pattern {3}", System.Environment.MachineName, jobId.ToString(), taskId.ToString(), pattern);
                System.Net.Mail.SmtpClient client = new System.Net.Mail.SmtpClient("mail-relay.iu.edu");
                client.EnableSsl = true;
                client.DeliveryMethod = System.Net.Mail.SmtpDeliveryMethod.Network;
                client.Send("SmithWatermanMS@indiana.edu", to, subjectLine, results);
            }
            catch
            {
            }
        }

        #endregion

        #region Configuration Read/Write

        static void ReadConfiguration(ConfigurationMgr mgr)
        {
            var swgmsSection = mgr.SmithWatermanMS;
            _fastaFile = swgmsSection.FastaFile;
            _distanceFile = swgmsSection.DistanceMatrixFile;
            _indexFile = swgmsSection.IndexFile;
            _timingFile = swgmsSection.TimingFile;
            _summaryFile = swgmsSection.SummaryFile;
            _writeFullResults = swgmsSection.WriteFullMatrix;
            _writeAlignments = swgmsSection.WriteAlignments;
            _writeAlignmentFile = swgmsSection.WriteAlignmentsFile;

            _blockWriteFrequency = swgmsSection.BlockWriteFrequency;
            _blockDir = swgmsSection.BlockDir;
            _logWriteFrequency = swgmsSection.LogWriteFrequency;
            _logDir = swgmsSection.LogDir;
            _restartFile = swgmsSection.RestartFile;
            _isRestart = !string.IsNullOrEmpty(_restartFile);
            _doneStatusDir = swgmsSection.DoneStatusDir;

            _nodeCount = swgmsSection.NodeCount;
            _processPerNodeCount = swgmsSection.ProcessPerNodeCount;

            _sequenceCount = 0;
            _gapOpen = swgmsSection.GapOpenPenalty;
            _gapExtension = swgmsSection.GapExtensionPenalty;

            _moleculeType = swgmsSection.MoleculeType;
            _distanceFunction = swgmsSection.DistanceFunctionType;
            _scoringMatrixName = swgmsSection.ScoringMatrixName;

            _projectName = swgmsSection.ProjectName;

            _emailResultsTo = string.Join(",", mgr.GlobalSection.EmailAddresses);
        }

        static void WriteConfiguration(ConfigurationMgr mgr)
        {
            var swgmsSection = mgr.SmithWatermanMS;
            swgmsSection.ProjectName = _projectName;
            swgmsSection.FastaFile = _fastaFile;
            swgmsSection.DistanceMatrixFile = _distanceFile;
            swgmsSection.IndexFile = _indexFile;
            swgmsSection.TimingFile = _timingFile;
            swgmsSection.SummaryFile = _summaryFile;
            swgmsSection.WriteFullMatrix = _writeFullResults;
            swgmsSection.WriteAlignments = _writeAlignments;
            swgmsSection.WriteAlignmentsFile = _writeAlignmentFile;

            swgmsSection.BlockWriteFrequency = _blockWriteFrequency;
            swgmsSection.BlockDir = _blockDir;
            swgmsSection.LogWriteFrequency = _logWriteFrequency;
            swgmsSection.LogDir = _logDir;
            swgmsSection.RestartFile = _restartFile;
            swgmsSection.DoneStatusDir = _doneStatusDir;

            swgmsSection.NodeCount = _nodeCount;
            swgmsSection.ProcessPerNodeCount = _processPerNodeCount;

            swgmsSection.SequenceCount = _sequenceCount;
            swgmsSection.GapOpenPenalty = _gapOpen;
            swgmsSection.GapExtensionPenalty = _gapExtension;

            swgmsSection.MoleculeType = _moleculeType;
            swgmsSection.DistanceFunctionType = _distanceFunction;
            swgmsSection.ScoringMatrixName = _scoringMatrixName;

            // Changes the MDS config
            mgr.ManxcatSection.DataPoints = _sequenceCount;
            mgr.ManxcatSection.DistanceMatrixFile = _distanceFile;
            mgr.ManxcatSection.IndexFile = _indexFile;

            // Changes the Pairwise config
            mgr.PairwiseSection.DataPoints = _sequenceCount;
            mgr.PairwiseSection.DistanceMatrixFile = _distanceFile;
            mgr.PairwiseSection.IndexFile = _indexFile;
        }

        #endregion

    }

    public class RowBlocksMetaData
    {
        private byte [] _projectName;
        private byte [] _sysName;
        private byte[] _fastaName;
        private int _blockWriteFrequency;
        private string _blockDir;
        
        public RowBlocksMetaData(string projectName, string sysName, string fastaName, int blockWriteFrequency, string blockDir)
        {
            ASCIIEncoding enc = new ASCIIEncoding();
            _projectName = enc.GetBytes(projectName);
            _sysName = enc.GetBytes(sysName);
            _fastaName = enc.GetBytes(fastaName);
            _blockWriteFrequency = blockWriteFrequency;
            _blockDir = blockDir;
        }

        public byte[] ProjectName
        {
            get { return _projectName; }
            set { _projectName = value; }
        }

        public byte[] SysName
        {
            get { return _sysName; }
            set { _sysName = value; }
        }

        public byte[] FastaName
        {
            get { return _fastaName; }
            set { _fastaName = value; }
        }

        public UInt16 ProjectNameLen
        {
            get { return (UInt16) _projectName.Length; }
        }

        public UInt16 SysNameLen
        {
            get { return (UInt16) _sysName.Length; }
        }

        public UInt16 FastaNameLen
        {
            get { return (UInt16) _fastaName.Length; }
        }

        public int BlockWriteFrequency
        {
            get { return _blockWriteFrequency; }
            set { _blockWriteFrequency = value; }
        }

        public string BlockDir
        {
            get { return _blockDir; }
            set { _blockDir = value; }
        }
    }

    public class MPIComputeLowerTriangularBlocked
    {
        public delegate int ComputeBlockHandler(Block block, PartialMatrix<TDistance> matrix, bool isDiagonal, int rank);

        private readonly int _size = 0;
        private readonly int _rank = 0;
        private long _blockComputeTime = 0;
        private long _mpiScatterTime = 0;
        private long _cellsComputed = 0;
        private PartialMatrix<TDistance> _partialMatrix;
        private readonly Range[] _processRowRanges;
        private readonly Block[][] _processBlocks;
        private readonly Range _processRowRange;
        private readonly ComputeBlockHandler _blockHander;
        private readonly RowBlocksMetaData _rowBlocksMetaData;
        private readonly int _worldSize;
        public MPIComputeLowerTriangularBlocked(int size, ComputeBlockHandler blockHandler, 
            RowBlocksMetaData rowBlocksMetaData, int rank, int worldsize)
        {
            _size = size;
            _worldSize = worldsize;
            _blockHander = blockHandler;
            _rank = rank;
            _processBlocks = BlockPartitioner.Partition(size, size, worldsize, worldsize);
            _processRowRanges = RangePartitioner.Partition(size, worldsize);
            _processRowRange = _processRowRanges[rank];
            _partialMatrix = new PartialMatrix<TDistance>(_processRowRange.StartIndex, _processRowRange.EndIndex, 0, size - 1);

            _rowBlocksMetaData = rowBlocksMetaData;
        }

        public PartialMatrix<TDistance> PartialMatrix
        {
            get
            {
                return _partialMatrix;
            }
        }

        public long CellsComputed
        {
            get
            {
                return _cellsComputed;
            }
        }

        public long BlockComputeTime
        {
            get
            {
                return _blockComputeTime;
            }
        }

        public long MPIScatterTime
        {
            get
            {
                return _mpiScatterTime;
            }
        }
        

        // Performs block output
        public void WriteBlocks(Block[] blocks, int lastSavedBlock, int lastUnsavedBlock)
        {
            ASCIIEncoding enc = new ASCIIEncoding();
            int firstUnsavedBlock = lastSavedBlock + 1;

            string blockFile;
            Block block;
            Range rowRange, colRange;
            for (int j = lastUnsavedBlock; j >= firstUnsavedBlock; j--)
            {
                // See if this block is actually computed, otherwise ignore it
                if ((_rank == j) || (_rank > j && !IsOdd(_rank + j)) || (_rank < j && IsOdd(_rank + j)))
                {
                    blockFile = Path.Combine(_rowBlocksMetaData.BlockDir, _rank + "-" + j + ".bin");
                    block = blocks[j];
                    rowRange = block.RowRange;
                    colRange = block.ColumnRange;

                    using (BinaryWriter writer = new BinaryWriter(File.Create(blockFile)))
                    {
                        WriteRowBlockMetaData(block, writer, rowRange, colRange, enc);
                        for (int m = rowRange.StartIndex; m <=rowRange.EndIndex; m++)
                        {
                            for (int n = colRange.StartIndex; n <= colRange.EndIndex; n++)
                            {
                                writer.Write(_partialMatrix[m,n]);
                            }
                        }
                    }
                }
                
            }

        }

        private void WriteRowBlockMetaData(Block block, BinaryWriter writer, Range rowRange, Range colRange, ASCIIEncoding enc)
        {
            // Project Name
            writer.Write(_rowBlocksMetaData.ProjectNameLen);
            writer.Write(_rowBlocksMetaData.ProjectName);
            
            // Rank
            writer.Write(((UInt16) _rank));

            // System Name
            writer.Write(_rowBlocksMetaData.SysNameLen);
            writer.Write(_rowBlocksMetaData.SysName);

            // FASTA Name
            writer.Write(_rowBlocksMetaData.FastaNameLen);
            writer.Write(_rowBlocksMetaData.FastaName);
            // RB#
            writer.Write(((UInt16) block.RowBlockNumber));
            // CB#
            writer.Write(((UInt16) block.ColumnBlockNumber));
            // Row Range Start Index
            writer.Write(((UInt32) block.RowRange.StartIndex));
            // Row Range End Index
            writer.Write(((UInt32) block.RowRange.EndIndex));
            // Col Range Start Index
            writer.Write(((UInt32) block.ColumnRange.StartIndex));
            // Col Range End Index
            writer.Write(((UInt32) block.ColumnRange.EndIndex));
            // Row Range Start Seq Name Length
            writer.Write(((UInt16) rowRange.StartSeqName.Length));
            // Row Range Start Seq Name
            writer.Write(enc.GetBytes(rowRange.StartSeqName));
            // Row Range End Seq Name Len
            writer.Write(((UInt16) rowRange.EndSeqName.Length));
            // Row Range End Seq Name
            writer.Write(enc.GetBytes(rowRange.EndSeqName));
            // Col Range Start Seq Name Len
            writer.Write(((UInt16) colRange.StartSeqName.Length));
            // Col Range Start Seq Name
            writer.Write(enc.GetBytes(colRange.StartSeqName));
            // Col Range End Seq Name Len
            writer.Write(((UInt16) colRange.EndSeqName.Length));
            // Col Range End Seq Name
            writer.Write(enc.GetBytes(colRange.EndSeqName));
        }

        // Perform the block computations

        public void Compute(int lastSavedBlock)
        {
            int totalComputableBlocks = IsOdd(_worldSize)
                                            ? ((_worldSize/2) + 1)
                                            : (IsOdd(_rank) ? (_worldSize/2) : ((_worldSize/2) + 1));
            int initiallyCompletedBlockCount = lastSavedBlock == -1
                                                   ? 0
                                                   : ((!IsOdd(_rank) && IsOdd(lastSavedBlock))
                                                          ? ((lastSavedBlock/2) + 2)
                                                          : ((lastSavedBlock/2) + 1));
            Console.WriteLine(Messages.StatusInitiallyComputed, _rank, initiallyCompletedBlockCount, totalComputableBlocks, (totalComputableBlocks - initiallyCompletedBlockCount));
            totalComputableBlocks -= initiallyCompletedBlockCount;
            int completedBlockPercentage;
            int completedBlockCount = 0;
            int unsavedBlockCount = 0;

            Stopwatch computeTimer = new Stopwatch();
            for (int j = (lastSavedBlock + 1); j < _processBlocks[_rank].Length; j++)
            {
                if ((_rank == j) || (_rank > j && !IsOdd(_rank + j)) || (_rank < j && IsOdd(_rank + j)))
                {
                    /* Working on a computable block */
                    Block block = _processBlocks[_rank][j];
                    // Set row and column block numbers for the block
                    block.SetIndex(_rank, j);
                    // Completed block percentage
                    completedBlockPercentage = 100 * completedBlockCount / totalComputableBlocks;

                    if (completedBlockCount % 5 == 0)
                    {
                        Console.WriteLine(Messages.StatusComputingBlocks, _rank, completedBlockPercentage,
                                          completedBlockCount + 1, totalComputableBlocks, block.RowBlockNumber,
                                          block.ColumnBlockNumber);
                    }

                    computeTimer.Start();
                    // _rank == j is the diagonal test
                    _cellsComputed += _blockHander(block, _partialMatrix, (_rank == j), _rank);
                    computeTimer.Stop();
                    _blockComputeTime += computeTimer.ElapsedMilliseconds;
                    computeTimer.Reset();

                    completedBlockCount++;
                    unsavedBlockCount++;

                    if (_rowBlocksMetaData.BlockWriteFrequency > 0 && unsavedBlockCount == _rowBlocksMetaData.BlockWriteFrequency)
                    {
                        WriteBlocks(_processBlocks[_rank], lastSavedBlock, j);
                        unsavedBlockCount = 0;
                        lastSavedBlock = j;
                    }
                }
            }

            Console.WriteLine(Messages.StatusComputingBlocksCompleted, _rank, completedBlockCount, totalComputableBlocks,
                              TimeSpan.FromMilliseconds(_blockComputeTime), DateTime.Now);

            // Flush any remaining unsaved blocks
            if (_rowBlocksMetaData.BlockWriteFrequency > 0 && unsavedBlockCount > 0)
            {
                WriteBlocks(_processBlocks[_rank], lastSavedBlock, _processBlocks[_rank].Length - 1);
            }

        }


        // Perform the MPI scatter
        public void Scatter()
        {
            Stopwatch scatterTimer = Stopwatch.StartNew();

            TDistance[][][] sentBlockValues = new TDistance[MPI.Intracommunicator.world.Size][][];

            _partialMatrix = _partialMatrix.Transpose();
            for (int k = 0; k < _processBlocks[_rank].Length; k++)
            {
                Block block = _processBlocks[_rank][k].Transpose();
                sentBlockValues[k] = _partialMatrix.GetBlockValues(block);
            }
            _partialMatrix = _partialMatrix.Transpose();

            for (int receivedRank = 0; receivedRank < MPI.Intracommunicator.world.Size; receivedRank++)
            {
                TDistance[][] transposedBlockValues = MPI.Intracommunicator.world.Scatter<TDistance[][]>(
                    sentBlockValues, receivedRank);
                int rowIndex = _rank;
                int columnIndex = receivedRank;

                if ((rowIndex < columnIndex) && !IsOdd(rowIndex + columnIndex))
                {
                    _partialMatrix.SetBlockValues(_processBlocks[rowIndex][columnIndex], transposedBlockValues);
                }

                if ((rowIndex > columnIndex) && IsOdd(rowIndex + columnIndex))
                {
                    _partialMatrix.SetBlockValues(_processBlocks[rowIndex][columnIndex], transposedBlockValues);
                }
            }

            _mpiScatterTime = scatterTimer.ElapsedMilliseconds;
        }

        public Matrix<TDistance> GetFullMatrix(int destinationRank)
        {
            Matrix<TDistance> result = null;

            if (_rank == destinationRank)
            {
                result = new Matrix<TDistance>(_size, _size);
            }

            // for each row
            for (int i = 0; i < _size; i++)
            {
                // for each process
                for (int sourceRank = 0; sourceRank < MPI.Intracommunicator.world.Size; sourceRank++)
                {
                    if (_processRowRanges[sourceRank].Contains(i))
                    {
                        if ((_rank == sourceRank) && (_rank == destinationRank))
                        {
                            // I'm the source and destination rank
                            TDistance[] values = _partialMatrix.GetRowValues(i);
                            result.SetRowValues(i, values);
                        }
                        else if ((_rank == sourceRank) && (_rank != destinationRank))
                        {
                            // I'm the source rank and I'm not the destination rank
                            TDistance[] values = _partialMatrix.GetRowValues(i);
                            MPI.Intracommunicator.world.Send<TDistance[]>(values, destinationRank, 100);
                        }
                        else if ((_rank == destinationRank) && (_rank != sourceRank))
                        {
                            // I'm the destination rank and I'm not the source rank
                            TDistance[] values = MPI.Intracommunicator.world.Receive<TDistance[]>(sourceRank, 100);
                            result.SetRowValues(i, values);
                        }
                    }
                }
            }

            return result;
        }

        public void WriteFullMatrixOnRank0(string fileName)
        {
            FileStream fileStream = null;
            BinaryWriter writer = null;

            int a = _size / MPI.Intercommunicator.world.Size;
            int b = _size % MPI.Intercommunicator.world.Size;

            /*
             * A note on row ranges and assigned process numbers.
             * First b number of process will have (a + 1) number of rows each.
             * The rest will have only 'a' number of rows. So if a row number, j,
             * falls inside the first set, i.e. j < (b * (a + 1)), then the rank 
             * of the process that handles this row is equal to the integer division
             * of j / (a + 1). Else, i.e. j >= (b * (a + 1)) then that row is 
             * in the second set of processes. Thus, the rank of the process handling
             * this row is equal to the integer calculation of b + [(j - (b * (a + 1)) / a]
             */

            int numOfRowsPerReceive = a;

            Range nextRowRange = null;

            if (_rank == 0)
            {
                fileStream = File.Create(fileName, 4194304);
                writer = new BinaryWriter(fileStream);

                // I am rank0 and I am the one who will fill the fullMatrix. So let's fill what I have already.
                for (int i = _partialMatrix.GlobalRowStartIndex; i <= _partialMatrix.GlobalRowEndIndex; i++)
                {
                    TDistance[] values = _partialMatrix.GetRowValues(i);
                    foreach (TDistance value in values)
                    {
                        writer.Write(value);
                    }
                }
            }

            // For all the remaining rows that rank0 does not have receive in blocks of rows
            for (int i = _processRowRanges[0].EndIndex + 1; i < _size; )
            {
                if (_rank == 0)
                {
                    // I am rank0 and let's declare the next row range that I want to receive.
                    int end = i + numOfRowsPerReceive - 1;
                    end = end >= _size ? _size - 1 : end;
                    nextRowRange = new Range(i, end);
                }

                // Announce everyone about the next row ranges that rank0 has declared.
                MPI.Intercommunicator.world.Broadcast<Range>(ref nextRowRange, 0);

                if (_rank == 0)
                {
                    /* I am rank0 and now let's try to receive the declared next row range from others */

                    // A variable to hold the rank of the process, which has the row that I am (rank0) going to receive
                    int processRank;

                    TDistance[] values;
                    for (int j = nextRowRange.StartIndex; j <= nextRowRange.EndIndex; j++)
                    {
                        // Let's find the process that has the row j.
                        processRank = j < (b * (a + 1)) ? j / (a + 1) : b + ((j - (b * (a + 1))) / a);

                        // For each row that I (rank0) require I will receive from the process, which has that row.
                        values = MPI.Intercommunicator.world.Receive<TDistance[]>(processRank, 100);

                        // Set the received values in the fullMatrix
                        foreach (TDistance value in values)
                        {
                            writer.Write(value);
                        }
                    }
                }
                else
                {
                    /* I am just an ordinary process and I am ready to give rank0 whatever the row it requests if I have that row */

                    // find the intersection of the row ranges of what I (the ordinary process) have and what rank0 wants and then send those rows to rank0
                    if (_processRowRange.IntersectsWith(nextRowRange))
                    {
                        Range intersection = _processRowRange.GetIntersectionWith(nextRowRange);
                        for (int k = intersection.StartIndex; k <= intersection.EndIndex; k++)
                        {
                            MPI.Intercommunicator.world.Send<TDistance[]>(_partialMatrix.GetRowValues(k), 0, 100);
                        }
                    }
                }

                i += numOfRowsPerReceive;
            }

            // I am rank0 and I came here means that I wrote full matrix to disk. So I will clear the writer and stream.
            if (_rank == 0)
            {
                writer.Flush();
                fileStream.Close();
                writer.Close();
            }

        }

        public void WriteFullMatrix(int destinationRank, string fileName)
        {
            if (destinationRank == _rank)
            {
                using (FileStream fileStream = File.Create(fileName, 4194304))
                {
                    using (BinaryWriter writer = new BinaryWriter(fileStream))
                    {
                        WriteFullMatrixInternal(destinationRank, writer);
                    }
                }
            }
            else
            {
                WriteFullMatrixInternal(destinationRank, null);
            }
        }

        public void WriteFullMatrixInternal(int destinationRank, BinaryWriter writer)
        {
            // for each row
            for (int i = 0; i < _size; i++)
            {
                // for each process
                for (int sourceRank = 0; sourceRank < MPI.Intracommunicator.world.Size; sourceRank++)
                {
                    if (_processRowRanges[sourceRank].Contains(i))
                    {
                        if ((_rank == sourceRank) && (_rank == destinationRank))
                        {
                            // I'm the source and destination rank
                            TDistance[] values = _partialMatrix.GetRowValues(i);

                            for (int index = 0; index < values.Length; index++)
                            {
                                writer.Write(values[index]);
                            }
                        }
                        else if ((_rank == sourceRank) && (_rank != destinationRank))
                        {
                            // I'm the source rank and I'm not the destination rank
                            TDistance[] values = _partialMatrix.GetRowValues(i);
                            MPI.Intracommunicator.world.Send<TDistance[]>(values, destinationRank, 100);
                        }
                        else if ((_rank == destinationRank) && (_rank != sourceRank))
                        {
                            // I'm the destination rank and I'm not the source rank
                            TDistance[] values = MPI.Intracommunicator.world.Receive<TDistance[]>(sourceRank, 100);

                            for (int index = 0; index < values.Length; index++)
                            {
                                writer.Write(values[index]);
                            }
                        }
                    }
                }
            }
        }

        private static bool IsOdd(int value)
        {
            return (value & 1) == 1;
        }
    }
    public static class MPIExtensions
    {
        public static int GetNumberOfNodes(this MPI.Intracommunicator communicator)
        {
            string[] processors = communicator.Allgather<string>(MPI.Environment.ProcessorName).Distinct().ToArray();
            return processors.Length;
        }
        public static int GetMinRankOnNode(this MPI.Intracommunicator communicator, string processorName)
        {
            string[] processers = communicator.Allgather<string>(MPI.Environment.ProcessorName);

            for (int i = 0; i < processers.Length; i++)
            {
                if (processers[i] == processorName)
                {
                    return i;
                }
            }

            throw new ApplicationException("An error occured finding the min rank on processor");
        }
        public static int GetMaxRankOnNode(this MPI.Intracommunicator communicator, string processorName)
        {
            string[] processers = communicator.Allgather<string>(MPI.Environment.ProcessorName);

            for (int i = processers.Length - 1; i >= 0; i--)
            {
                if (processers[i] == processorName)
                {
                    return i;
                }
            }

            throw new ApplicationException("An error occured finding the max rank on processor");
        }
        public static int[] GetMinRankOnAllNodes(this MPI.Intracommunicator communicator)
        {
            string[] processors = communicator.Allgather<string>(MPI.Environment.ProcessorName).Distinct().ToArray();
            int[] minRanks = new int[processors.Length];

            for (int i = 0; i < processors.Length; i++)
            {
                minRanks[i] = GetMinRankOnNode(communicator, processors[i]);
            }

            return minRanks;
        }
        public static int[] GetMaxRankOnAllNodes(this MPI.Intracommunicator communicator)
        {
            string[] processors = communicator.Allgather<string>(MPI.Environment.ProcessorName).Distinct().ToArray();
            int[] maxRanks = new int[processors.Length];

            for (int i = 0; i < processors.Length; i++)
            {
                maxRanks[i] = GetMaxRankOnNode(communicator, processors[i]);
            }

            return maxRanks;
        }
    }
}
