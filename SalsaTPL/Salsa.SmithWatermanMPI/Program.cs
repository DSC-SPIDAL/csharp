using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Text;
using System.Linq;
using Salsa.Core;
using Salsa.Core.Blas;
using Salsa.Core.Configuration;
using Salsa.Core.Bio.Algorithms;
using Salsa.Core.Bio.Algorithms.SmithWaterman;
using Salsa.Core.Bio.IO;
using System.Threading.Tasks;

#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using TDistance = System.Int16;
#else
using TDistance = System.Double;
#endif

namespace Salsa.SmithWatermanTPL
{
    public class Program
    {
        private static int _sequenceCount = 0;
        private static int _nodeCount = 1;
        private static int _processPerNodeCount = 1;
        private static int _threadPerProcessCount = 1;

        private static long _overallTime = 0;
        private static long _overallComputeTime = 0;
        private static long _overallMPIScatterTime = 0;
        private static long _overallDiskIOTime = 0;

        private static long _totalPairCount = 0;
        private static long _totalComputeTime = 0;
        private static long _totalMPIScatterTime = 0;

        private static string _jobId = string.Empty;
        private static string _taskId = string.Empty;

        private static string _fastaFile = string.Empty;
        private static string _distanceFile = string.Empty;
        private static string _indexFile = string.Empty;
        private static string _timingFile = string.Empty;
        private static string _summaryFile = string.Empty;
        private static string _configFile = string.Empty;
        private static bool _writeFullResults = true;
        private static bool _writePartialResults = false;

        private static bool _writeAlignments = false;
        private static string _writeAlignmentFile = string.Empty;
        private static StreamWriter _writeAlignmentsWriter = null;

        private static float _gapOpen = 14;
        private static float _gapExtension = 4;
        private static string _scoringMatrixName = "EDNAFULL";
        private static AlignmentType _alignmentType = AlignmentType.Nucleic;
        private static DistanceFunctionType _distanceFunction = DistanceFunctionType.PercentIdentity;
        private static string _emailResultsTo = "";

        private static IList<Sequence> _sequences;
        private static ConfigurationMgr _mgr;

        static void Main(string[] args)
        {
            // Load the command line args into our helper class which allows us to name arguments
            Arguments pargs = new Arguments(args);
            pargs.Usage = "Usage: Salsa.SmithWatermanMPI.exe /configFile=<string> /nodeCount=<int>";

            if (pargs.CheckRequired(new string[] { "configFile", "nodeCount" }) == false)
            {
                Console.WriteLine(pargs.Usage);
                return;
            }

            using (MPI.Environment env = new MPI.Environment(ref args))
            {
                Stopwatch overallTimer = Stopwatch.StartNew();
                _configFile = pargs.GetValue<string>("configFile");

                _mgr = ConfigurationMgr.LoadConfiguration(_configFile, true);
                ReadConfiguration(_mgr);

                _jobId = Environment.GetEnvironmentVariable("%CCP_JOBID%");
                _taskId = Environment.GetEnvironmentVariable("%CCP_TASKID%");
                _nodeCount = pargs.GetValue<int>("nodeCount");
                _processPerNodeCount = MPI.Intracommunicator.world.Size / _nodeCount;

                // Read the fasta input file using our parser.      
                FastaSequenceParser parser = new FastaSequenceParser();
                _sequences = parser.Parse(_fastaFile);
                _sequenceCount = _sequences.Count;

                Console.WriteLine("Read {0} sequence from Fasta File {1}", _sequenceCount, _fastaFile);

                if ((_writeAlignments == true) && (!string.IsNullOrWhiteSpace(_writeAlignmentFile)))
                {
                    string fileName = string.Format("{0}_{1}{2}", Path.GetFileNameWithoutExtension(_writeAlignmentFile), MPI.Intracommunicator.world.Rank, Path.GetExtension(_writeAlignmentFile));
                    string filePath = Path.Combine(Path.GetDirectoryName(_writeAlignmentFile), fileName);
                    _writeAlignmentsWriter = File.CreateText(filePath);
                }

                MPIComputeLowerTriangularBlocked computeAlignments = new MPIComputeLowerTriangularBlocked(_sequenceCount, ComputeBlock);

                // Block Compute Section
                MPI.Intracommunicator.world.Barrier();
                Stopwatch computeTimer = Stopwatch.StartNew();

                Console.WriteLine("Rank {0} Beginning Blocks Compute: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());
                computeAlignments.Compute();
                Console.WriteLine("Rank {0} Finished Blocks Compute: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());

                MPI.Intracommunicator.world.Barrier();
                _overallComputeTime = computeTimer.ElapsedMilliseconds;
                computeTimer.Stop();

                if (_writeAlignmentsWriter != null)
                {
                    _writeAlignmentsWriter.Flush();
                    _writeAlignmentsWriter.Dispose();
                    return;
                }

                // MPI Scatter Section
                MPI.Intracommunicator.world.Barrier();
                Stopwatch mpiScatterTimer = Stopwatch.StartNew();

                Console.WriteLine("Rank {0} Beginning Scatter Blocks: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());
                computeAlignments.Scatter();
                Console.WriteLine("Rank {0} Finished Scatter Blocks: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());

                MPI.Intracommunicator.world.Barrier();

                _overallMPIScatterTime = mpiScatterTimer.ElapsedMilliseconds;
                mpiScatterTimer.Stop();

                Stopwatch ioTimer = Stopwatch.StartNew();
                Console.WriteLine("Rank {0} Beginning File IO: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());
                {

                    if (_writePartialResults)
                    {
                        string partialDistanceFile = string.Format("{0}_{1}{2}", Path.GetFileNameWithoutExtension(_distanceFile), MPI.Intercommunicator.world.Rank, Path.GetExtension(_distanceFile));
                        partialDistanceFile = Path.Combine(Path.GetDirectoryName(_distanceFile), partialDistanceFile);
                        PartialMatrixBinaryWriter partialWriter = new PartialMatrixBinaryWriter();
                        partialWriter.Write(computeAlignments.PartialMatrix, partialDistanceFile);
                    }


                    // todo: saliya - Workaround code for issue with Tempest. Remove this once Tempest is healthy.
                    if (_writeFullResults)
                    {
                        computeAlignments.WriteFullMatrixOnRank0(_distanceFile);
                    }

                    // todo: saliya - This code doesn't work on Tempest with all cores for large data sizes. 
                    //                Something is wrong with MPI in Tempest probably. Till then we will use the 
                    //                code above.
                    /*
                    if (_writeFullResults)
                    {
                        foreach (int minRankOnNode in MPI.Intracommunicator.world.GetMinRankOnAllNodes())
                        {
                            computeAlignments.WriteFullMatrix(minRankOnNode, _distanceFile);
                            //    Matrix<TDistance> fullMatrix = computeAlignments.GetFullMatrix(minRankOnNode);

                            //    if (MPI.Intracommunicator.world.Rank == minRankOnNode)
                            //    {
                            //        MatrixBinaryWriter fullWriter = new MatrixBinaryWriter();
                            //        fullWriter.Write(fullMatrix, _distanceFile);
                            //        Console.WriteLine("\tWrote Distance Matrix: {0}\t{1}", Environment.MachineName, _distanceFile);
                            //    }
                        }
                    }*/
                }
                Console.WriteLine("Rank {0} Finished File IO: {1}", MPI.Intracommunicator.world.Rank, DateTime.Now.ToString());
                MPI.Intracommunicator.world.Barrier();
                _overallDiskIOTime = ioTimer.ElapsedMilliseconds;
                ioTimer.Stop();

                _overallTime = overallTimer.ElapsedMilliseconds;
                overallTimer.Stop();

                _totalPairCount = MPI.Intercommunicator.world.Reduce<long>(computeAlignments.CellsComputed, MPI.Operation<long>.Add, 0);
                _totalComputeTime = MPI.Intracommunicator.world.Reduce<long>(computeAlignments.BlockComputeTime, MPI.Operation<long>.Add, 0);
                _totalMPIScatterTime = MPI.Intracommunicator.world.Reduce<long>(computeAlignments.MPIScatterTime, MPI.Operation<long>.Add, 0);

                if (MPI.Intracommunicator.world.Rank == 0)
                {
                    WriteIndexFile(_indexFile, _sequences);
                    WriteTimingFile(_timingFile);
                    WriteSummaryFile(_summaryFile);
                    WriteConfiguration(_mgr);
                    _mgr.SaveAs(_configFile);
                }
            }
        }

        static int ComputeBlock(Block block, PartialMatrix<TDistance> matrix, bool isDiagonal)
        {
            int count = 0;
            SmithWatermanGotoh aligner = new SmithWatermanGotoh(_alignmentType);
            ScoringMatrix scoringMatrix = ScoringMatrix.Load(_scoringMatrixName);

            if (isDiagonal)
            {
                for (int i = block.RowRange.StartIndex; i <= block.RowRange.EndIndex; i++)
                {
                    Sequence si = _sequences[i];

                    for (int j = block.ColumnRange.StartIndex; j < i; j++)
                    {
                        Sequence sj = _sequences[j];

                        PairwiseAlignment alignment = aligner.Align(si, sj, scoringMatrix, _gapOpen, _gapExtension);

                        if ((_writeAlignments) && (_writeAlignmentsWriter != null))
                        {
                            WriteAlignment(_writeAlignmentsWriter, alignment);
                        }

                        matrix[i, j] = matrix[j, i] = ComputeDistance(alignment);

                        count++;
                    }
                }
            }
            else
            {
                for (int i = block.RowRange.StartIndex; i <= block.RowRange.EndIndex; i++)
                {
                    Sequence si = _sequences[i];

                    for (int j = block.ColumnRange.StartIndex; j <= block.ColumnRange.EndIndex; j++)
                    {
                        Sequence sj = _sequences[j];

                        PairwiseAlignment alignment = aligner.Align(si, sj, scoringMatrix, _gapOpen, _gapExtension);

                        if ((_writeAlignments) && (_writeAlignmentsWriter != null))
                        {
                            WriteAlignment(_writeAlignmentsWriter, alignment);
                        }

                        matrix[i, j] = ComputeDistance(alignment);

                        count++;
                    }
                }
            }

            return count;
        }

#if USE_UINT16 || USE_INT16
        static TDistance ComputeDistance(PairwiseAlignment alignment)
        {
            TDistance result = 0;

            if (_alignmentType == AlignmentType.Nucleic)
            {
                switch (_distanceFunction)
                {
                    case DistanceFunctionType.PercentIdentity:
                        result = (TDistance)((1.0f - alignment.PercentIdentity) * TDistance.MaxValue);
                        break;
                    case DistanceFunctionType.Kimura2:
                        result = (TDistance)(alignment.ComputeKimuraDistance() * TDistance.MaxValue);
                        break;
                    case DistanceFunctionType.JukesCantor:
                        result = (TDistance)(alignment.ComputeJukesCantorDistance() * TDistance.MaxValue);
                        break;
                }
            }
            else
            {
                result = (TDistance)((1.0f - alignment.PercentSimilarity) * TDistance.MaxValue);
            }

            return result;
        }

#else
        static TDistance ComputeDistance(PairwiseAlignment alignment)
        {
            TDistance result = 0;

            if (_alignmentType == AlignmentType.Nucleic)
            {
                switch (_distanceFunction)
                {
                    case DistanceFunctionType.PercentIdentity:
                        result = (TDistance)(1.0f - alignment.PercentIdentity);
                        break;
                    case DistanceFunctionType.Kimura2:
                        result = (TDistance)(alignment.ComputeKimuraDistance());
                        break;
                    case DistanceFunctionType.JukesCantor:
                        result = (TDistance)(alignment.ComputeJukesCantorDistance());
                        break;
                }
            }
            else
            {
                result = (TDistance)(1.0f - alignment.PercentSimilarity);
            }

            return result;
        }
#endif

        static void WriteIndexFile(string fileName, IEnumerable<Sequence> sequences)
        {
            using (StreamWriter writer = new StreamWriter(fileName))
            {
                int index = 0;

                foreach (Sequence sequence in sequences)
                {
                    writer.WriteLine("{0}\t{1}\t{2}", index, sequence.Label, sequence.Residues.Length);
                    index++;
                }

                writer.Flush();
                writer.Close();
            }
        }

        static void WriteSummaryFile(string fileName)
        {
            StringBuilder sb = new StringBuilder();
            string pattern = string.Format("{0}x{1}x{2}", _threadPerProcessCount, _processPerNodeCount, _nodeCount);
            int parallelism = _nodeCount * _processPerNodeCount * _threadPerProcessCount;
            sb.AppendFormat("Run Date: {0}", DateTime.Now.ToLongDateString());
            sb.AppendLine();
            sb.AppendFormat("Machine: {0}", System.Environment.MachineName);
            sb.AppendLine();
            sb.AppendFormat("Job Id: {0}", _jobId);
            sb.AppendLine();
            sb.AppendFormat("Task Id: {0}", _taskId);
            sb.AppendLine();
            sb.AppendFormat("Overall Time: {0} ms, {1}", _overallTime, TimeSpan.FromMilliseconds(_overallTime).ToString());
            sb.AppendLine();
            sb.AppendFormat("Overall Block Compute Time: {0} ms, {1}", _overallComputeTime, TimeSpan.FromMilliseconds(_overallComputeTime).ToString());
            sb.AppendLine();
            sb.AppendFormat("Overall MPI Scatter Time: {0} ms, {1}", _overallMPIScatterTime, TimeSpan.FromMilliseconds(_overallMPIScatterTime).ToString());
            sb.AppendLine();
            sb.AppendFormat("Overall Disk IO Time: {0} ms, {1}", _overallDiskIOTime, TimeSpan.FromMilliseconds(_overallDiskIOTime).ToString());
            sb.AppendLine();
            sb.AppendFormat("Total Block Compute Duration (ms): {0}", _totalComputeTime);
            sb.AppendLine();
            sb.AppendFormat("Total MPI Scatter Duration (ms): {0}", _totalMPIScatterTime);
            sb.AppendLine();
            sb.AppendFormat("Sequence Count: {0}", _sequenceCount);
            sb.AppendLine();
            sb.AppendFormat("Pairwise Alignments Count: {0}", _totalPairCount);
            sb.AppendLine();
            sb.AppendFormat("Pattern: {0}", pattern);
            sb.AppendLine();
            sb.AppendFormat("Threads Per Process Count: {0}", _threadPerProcessCount);
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
            sb.AppendFormat("AlignmentType: {0}", _alignmentType);
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
            string pattern = string.Format("{0}x{1}x{2}", _threadPerProcessCount, _processPerNodeCount, _nodeCount);
            int parallelism = _nodeCount * _processPerNodeCount * _threadPerProcessCount;

            if (File.Exists(fileName))
            {
                using (StreamWriter writer = new StreamWriter(fileName, true))
                {
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                        _jobId,
                        _taskId,
                        pattern,
                        parallelism,
                        _totalPairCount,
                        _overallTime,
                        _overallComputeTime,
                        _overallMPIScatterTime,
                        _overallDiskIOTime,
                        _totalComputeTime,
                        _totalMPIScatterTime,
                        _threadPerProcessCount.ToString(),
                        _processPerNodeCount.ToString(),
                        _nodeCount.ToString(),
                        System.Environment.MachineName,
                        DateTime.Now.ToString()
                        );
                }
            }
            else
            {
                using (StreamWriter writer = new StreamWriter(fileName, true))
                {
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                    "Job Id",                                   //0
                    "TaskId",                                   //1
                    "Pattern",                                  //2
                    "Parallelism",                              //3
                    "Pairs Aligned",                            //4
                    "Overall Duration (ms)",                    //5
                    "Overall Block Compute Duration (ms)",      //6
                    "Overall MPI Scatter Duration (ms)",        //7
                    "Overall DiskIO Duration (ms)",             //8
                    "Total Block Compute Duration (ms)",        //9
                    "Total MPI Scatter Duration (ms)",          //10
                    "Threads",                                  //11
                    "Processes",                                //12
                    "Nodes",                                    //13
                    "Machine Name",                             //14
                    "Run Date"                                  //15
                    );

                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                        _jobId,
                        _taskId,
                        pattern,
                        parallelism,
                        _totalPairCount,
                        _overallTime,
                        _overallComputeTime,
                        _overallMPIScatterTime,
                        _overallDiskIOTime,
                        _totalComputeTime,
                        _totalMPIScatterTime,
                        _threadPerProcessCount.ToString(),
                        _processPerNodeCount.ToString(),
                        _nodeCount.ToString(),
                        System.Environment.MachineName,
                        DateTime.Now.ToString()
                        );
                }
            }
        }

        static void EmailResults(string to, string jobId, string taskId, string pattern, string results)
        {
            try
            {
                string subjectLine = string.Format("SmithWatermanTPL Timing: Machine {0} Job {1} Task {2} Pattern {3}", System.Environment.MachineName, jobId.ToString(), taskId.ToString(), pattern);
                System.Net.Mail.SmtpClient client = new System.Net.Mail.SmtpClient("mail-relay.iu.edu");
                client.EnableSsl = true;
                client.DeliveryMethod = System.Net.Mail.SmtpDeliveryMethod.Network;
                client.Send("SmithWatermanTPL@indiana.edu", to, subjectLine, results);
            }
            catch
            {
            }
        }

        static void ReadConfiguration(ConfigurationMgr mgr)
        {
            _fastaFile = mgr.SmithWatermanSection.FastaFile;
            _distanceFile = mgr.SmithWatermanSection.DistanceMatrixFile;
            _indexFile = mgr.SmithWatermanSection.IndexFile;
            _timingFile = mgr.SmithWatermanSection.TimingFile;
            _summaryFile = mgr.SmithWatermanSection.SummaryFile;
            _writeFullResults = mgr.SmithWatermanSection.WriteFullMatrix;
            _writePartialResults = mgr.SmithWatermanSection.WritePartialMatrix;
            _writeAlignments = mgr.SmithWatermanSection.WriteAlignments;
            _writeAlignmentFile = mgr.SmithWatermanSection.WriteAlignmentsFile;

            _nodeCount = mgr.SmithWatermanSection.NodeCount;
            _processPerNodeCount = mgr.SmithWatermanSection.ProcessPerNodeCount;

            _sequenceCount = 0;
            _gapOpen = mgr.SmithWatermanSection.GapOpenPenalty;
            _gapExtension = mgr.SmithWatermanSection.GapExtensionPenalty;

            _alignmentType = mgr.SmithWatermanSection.AlignmentType;
            _distanceFunction = mgr.SmithWatermanSection.DistanceFunctionType;
            _scoringMatrixName = mgr.SmithWatermanSection.ScoringMatrixName;

            _emailResultsTo = string.Join(",", mgr.GlobalSection.EmailAddresses);
        }

        static void WriteConfiguration(ConfigurationMgr mgr)
        {
            mgr.SmithWatermanSection.FastaFile = _fastaFile;
            mgr.SmithWatermanSection.DistanceMatrixFile = _distanceFile;
            mgr.SmithWatermanSection.IndexFile = _indexFile;
            mgr.SmithWatermanSection.TimingFile = _timingFile;
            mgr.SmithWatermanSection.SummaryFile = _summaryFile;
            mgr.SmithWatermanSection.WriteFullMatrix = _writeFullResults;
            mgr.SmithWatermanSection.WritePartialMatrix = _writePartialResults;
            mgr.SmithWatermanSection.WriteAlignments = _writeAlignments;
            mgr.SmithWatermanSection.WriteAlignmentsFile = _writeAlignmentFile;

            mgr.SmithWatermanSection.NodeCount = _nodeCount;
            mgr.SmithWatermanSection.ProcessPerNodeCount = _processPerNodeCount;

            mgr.SmithWatermanSection.SequenceCount = _sequenceCount;
            mgr.SmithWatermanSection.GapOpenPenalty = _gapOpen;
            mgr.SmithWatermanSection.GapExtensionPenalty = _gapExtension;

            mgr.SmithWatermanSection.AlignmentType = _alignmentType;
            mgr.SmithWatermanSection.DistanceFunctionType = _distanceFunction;
            mgr.SmithWatermanSection.ScoringMatrixName = _scoringMatrixName;

            // Changes the MDS config
            mgr.ManxcatSection.DataPoints = _sequenceCount;
            mgr.ManxcatSection.DistanceMatrixFile = _distanceFile;
            mgr.ManxcatSection.IndexFile = _indexFile;

            // Changes the Pairwise config
            mgr.PairwiseSection.DataPoints = _sequenceCount;
            mgr.PairwiseSection.DistanceMatrixFile = _distanceFile;
            mgr.PairwiseSection.IndexFile = _indexFile;
        }

        static void WriteAlignment(StreamWriter writer, PairwiseAlignment alignment)
        {
            double kimeraDistance = alignment.ComputeKimuraDistance();
            double jukesCantorDistance = alignment.ComputeJukesCantorDistance();
            double percentIdentityDistance = alignment.PercentIdentity;

            writer.WriteLine(">{0} Kimera: {1}, JukesCantor: {2}, PercentIdentity: {3}", alignment.Name1, kimeraDistance.ToString("F5"), jukesCantorDistance.ToString("F5"), alignment.PercentIdentity.ToString("F5"));
            writer.WriteLine(alignment.Sequence1);
            writer.WriteLine(">{0} Kimera: {1}, JukesCantor: {2}, PercentIdentity: {3}", alignment.Name2, kimeraDistance.ToString("F5"), jukesCantorDistance.ToString("F5"), alignment.PercentIdentity.ToString("F5"));
            writer.WriteLine(alignment.Sequence2);
            writer.WriteLine();
        }
    }

    public class MPIComputeLowerTriangularBlocked
    {
        public delegate int ComputeBlockHandler(Block block, PartialMatrix<TDistance> matrix, bool isDiagonal);

        private int _size = 0;
        private int _rank = 0;
        private long _blockComputeTime = 0;
        private long _mpiScatterTime = 0;
        private long _cellsComputed = 0;
        private PartialMatrix<TDistance> _partialMatrix;
        private Range[] _processRowRanges;
        private Block[][] _processBlocks;
        private Range _processRowRange;
        private ComputeBlockHandler _blockHander;

        public MPIComputeLowerTriangularBlocked(int size, ComputeBlockHandler blockHandler)
        {
            _size = size;
            _blockHander = blockHandler;
            _rank = MPI.Intracommunicator.world.Rank;
            _processBlocks = BlockPartitioner.Partition(size, size, MPI.Intracommunicator.world.Size, MPI.Intracommunicator.world.Size);
            _processRowRanges = RangePartitioner.Partition(size, MPI.Intracommunicator.world.Size);
            _processRowRange = _processRowRanges[MPI.Intracommunicator.world.Rank];
            _partialMatrix = new PartialMatrix<TDistance>(_processRowRange.StartIndex, _processRowRange.EndIndex, 0, size - 1);
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

        // Perform the block computations
        public void Compute()
        {
            Stopwatch computeTimer = Stopwatch.StartNew();
            for (int j = 0; j < _processBlocks[_rank].Length; j++)
            {
                if ((j % 100) == 0)
                {
                    Console.WriteLine("\tRank {0} is Computing Block {1} of {2} - {3}", MPI.Intracommunicator.world.Rank, j, _processBlocks[_rank].Length, _processBlocks[_rank][j].ToString());
                }

                if (_rank == j)
                {
                    _cellsComputed += _blockHander(_processBlocks[_rank][j], _partialMatrix, true);
                }
                else if (_rank > j && !IsOdd(_rank + j))
                {
                    _cellsComputed += _blockHander(_processBlocks[_rank][j], _partialMatrix, false);
                }
                else if (_rank < j && IsOdd(_rank + j))
                {
                    _cellsComputed += _blockHander(_processBlocks[_rank][j], _partialMatrix, false);
                }
            }
            _blockComputeTime = computeTimer.ElapsedMilliseconds;
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
                TDistance[][] transposedBlockValues = MPI.Intracommunicator.world.Scatter<TDistance[][]>(sentBlockValues, receivedRank);
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
    }

    [Serializable]
    public struct RowData
    {
        int globalRrowIndex;
        TDistance[] rowValues;

        public RowData(int index, TDistance[] values)
        {
            globalRrowIndex = index;
            rowValues = values;
        }

        public int GlobalRowIndex
        {
            get
            {
                return globalRrowIndex;
            }
        }

        public TDistance[] RowValues
        {
            get
            {
                return rowValues;
            }
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

