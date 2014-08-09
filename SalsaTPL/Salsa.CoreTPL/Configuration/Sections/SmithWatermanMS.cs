using System;
using System.ComponentModel;
using System.IO;
using System.Reflection;
using System.Text;
using Bio.IO.GenBank;
using Bio.SimilarityMatrices;
using Salsa.Core.Bio.Algorithms;

namespace Salsa.Core.Configuration.Sections
{
    public class SmithWatermanMS : Section
    {
        #region Members

        private string _distanceFile = @"distance.bin";
        private DistanceFunctionType _distanceFunction = DistanceFunctionType.PercentIdentity;
        private string _fastaFile = "$(InputRootPath)";
        private int _gapExtension = 4;
        private int _gapOpen = 14;
        private string _indexFile = @"$(OutputRootPath)\index.txt";
        private MoleculeType _moleculeType = MoleculeType.DNA;
        private int _nodeCount;
        private int _processPerNodeCount;
        private string _scoringMatrixName = "EDNAFULL";
        private int _sequenceCount;
        private string _summaryFile = @"$(OutputRootPath)\swms_summary.txt";
        private string _timingFile = @"$(OutputRootPath)\swms_timing.txt";
        private string _writeAlignmentFile = @"$(OutputRootPath)\swms_alignments.txt";
        private bool _writeFullMatrix = true;

        /* Note. Saliya - Adding Nblock save and logs as suggested by Prof. Fox on 9/12/2011. See mail "64-bit MBF" Search for "introduce Nbackup set to" */
        // The number of column blocks computed for each partial write. 
        // Default = 0, indicates no partial writing, is not adviced for large runs

        // Directory to write blocks from each process. Adviced is to use a directory local to each process

        // The number of pairs to be calculated before writing status.
        // Default = 0, indicates no writing of status rankwsie status, is not adviced for large runs

        // Directory to write logs for each process. Adviced is to use a directory local to each process

        /* -- End Nblock save and log  members --*/

        #endregion

        #region Properties

        public string ProjectName { get; set; }

        [Category("I/O")]
        public string DoneStatusDir { get; set; }

        [Category("I/O")]
        public string RestartFile { get; set; }

        [Category("I/O")]
        public int BlockWriteFrequency { get; set; }

        [Category("I/O")]
        public int LogWriteFrequency { get; set; }

        [Category("I/O")]
        public string BlockDir { get; set; }

        [Category("I/O")]
        public string LogDir { get; set; }

        [Category("I/O")]
        public string FastaFile
        {
            get { return _fastaFile; }
            set
            {
                _fastaFile = value;
                OnPropertyChanged("FastaFile");
            }
        }

        [Category("I/O")]
        [Description("Full Path to the location where the distance file is to be written to. Macros will be expanded.")]
        public string DistanceMatrixFile
        {
            get { return _distanceFile; }
            set
            {
                _distanceFile = value;
                OnPropertyChanged("DistanceFile");
            }
        }

        [Category("I/O")]
        [Description("Full Path to the location where the index file is to be written to. Macros will be expanded.")]
        public string IndexFile
        {
            get { return _indexFile; }
            set
            {
                _indexFile = value;
                OnPropertyChanged("IndexFile");
            }
        }

        [Category("I/O")]
        public string TimingFile
        {
            get { return _timingFile; }
            set
            {
                _timingFile = value;
                OnPropertyChanged("TimingFile");
            }
        }

        [Category("I/O")]
        public string SummaryFile
        {
            get { return _summaryFile; }
            set
            {
                _summaryFile = value;
                OnPropertyChanged("SummaryFile");
            }
        }

        [Category("I/O")]
        [Description("The file where the alignments are written to.")]
        public string WriteAlignmentsFile
        {
            get { return _writeAlignmentFile; }
            set
            {
                _writeAlignmentFile = value;
                OnPropertyChanged("WriteAlignmentsFile");
            }
        }

        [TypeConverter(typeof (ScoringMatrixTypeConverter))]
        [Description("The name of scoring matrix to use")]
        public string ScoringMatrixName
        {
            get { return _scoringMatrixName; }
            set
            {
                _scoringMatrixName = value;
                OnPropertyChanged("ScoringMatrixName");
            }
        }

        [Description("The type of molecule")]
        public MoleculeType MoleculeType
        {
            get { return _moleculeType; }
            set { _moleculeType = value; }
        }

        [Description("The distance function to use when coverting alignments to distance measures.")]
        public DistanceFunctionType DistanceFunctionType
        {
            get { return _distanceFunction; }
            set
            {
                _distanceFunction = value;
                OnPropertyChanged("DistanceFunctionType");
            }
        }

        [Browsable(false)]
        public int NodeCount
        {
            get { return _nodeCount; }
            set
            {
                _nodeCount = value;
                OnPropertyChanged("NodeCount");
            }
        }

        [Browsable(false)]
        public int ProcessPerNodeCount
        {
            get { return _processPerNodeCount; }
            set
            {
                _processPerNodeCount = value;
                OnPropertyChanged("ProcessPerNodeCount");
            }
        }

        [Browsable(false)]
        public int SequenceCount
        {
            get { return _sequenceCount; }
            set
            {
                _sequenceCount = value;
                OnPropertyChanged("SequenceCount");
            }
        }

        [Description(
            "If true the full matrix is written at the specified location by the minimum rank mpi process on a given node."
            )]
        public bool WriteFullMatrix
        {
            get { return _writeFullMatrix; }
            set
            {
                _writeFullMatrix = value;
                OnPropertyChanged("WriteFullMatrix");
            }
        }

        [Description(
            "If true the alignment for each pair is written to the file specified in the WriteAlignmentsPath property.")
        ]
        public bool WriteAlignments { get; set; }

        [Description("The cost associated with opening a gap in the alignment.")]
        public int GapOpenPenalty
        {
            get { return _gapOpen; }
            set
            {
                _gapOpen = value;
                OnPropertyChanged("GapOpen");
            }
        }

        [Description("The cost associated with extending a gap in the alignment.")]
        public int GapExtensionPenalty
        {
            get { return _gapExtension; }
            set
            {
                _gapExtension = value;
                OnPropertyChanged("GapExtensionPenalty");
            }
        }

        #endregion

        internal override void ExpandMacro(MacroReplacement macroReplacement)
        {
            base.ExpandMacro(macroReplacement);
            FastaFile = macroReplacement.Expand(FastaFile);
            DistanceMatrixFile = macroReplacement.Expand(DistanceMatrixFile);
            IndexFile = macroReplacement.Expand(IndexFile);
            TimingFile = macroReplacement.Expand(TimingFile);
            SummaryFile = macroReplacement.Expand(SummaryFile);
            //WriteAlignmentsFile = macroReplacement.Expand(WriteAlignmentsFile);
        }

        /// <summary>
        /// Loads a scoring matrix from the predefined set of matrices inside Salsa.Core.Bio.Algorithms.Matrices
        /// </summary>
        /// <param name="matrixName">The name of the matrix</param>
        /// <param name="moleculeType">
        /// Type of molecule for which this matrix is designed. 
        /// Must be DNA, RNA (may have variants like tRNA, mRNA, etc.) or Protein</param>
        /// <returns>The SimilarityMatrix</returns>
        public SimilarityMatrix LoadSimilarityMatrix(string matrixName, MoleculeType moleculeType)
        {
            /*
             * MBF 2.0 requires the format of the matrix to be as (without angle brackets), 
             * <Name>
             * <MoleculeType>
             * <Alphabet>
             * <ScoreRow0>
             * <ScoreRow1>
             * ...
             * ...
             * <ScoreRowN>
             */

            if (moleculeType == MoleculeType.DNA || moleculeType == MoleculeType.mRNA ||
                moleculeType == MoleculeType.RNA ||
                moleculeType == MoleculeType.rRNA || moleculeType == MoleculeType.snoRNA ||
                moleculeType == MoleculeType.snRNA ||
                moleculeType == MoleculeType.tRNA || moleculeType == MoleculeType.uRNA ||
                moleculeType == MoleculeType.Protein)
            {
                using (
                    Stream stream =
                        Assembly.GetExecutingAssembly().GetManifestResourceStream(
                            "Salsa.Core.Bio.Algorithms.Matrices." + matrixName))
                {
                    using (var reader = new StreamReader(stream))
                    {
                        char commentStarter = '#';
                        string line;
                        // Skip comments
                        while ((line = reader.ReadLine()) != null && line.Trim()[0] == commentStarter)
                            ;

                        var sb = new StringBuilder();
                        sb.AppendLine(matrixName); // Matrix name

                        string mt = moleculeType.ToString();
                        sb.AppendLine((moleculeType == MoleculeType.Protein) ? mt : mt.Substring(mt.Length - 3));
                            // Molecule Type

                        sb.AppendLine(line); // Alphabet line

                        while ((line = reader.ReadLine()) != null)
                        {
                            sb.AppendLine(line.Substring(1).Trim());
                                // ScoreRow i (ignores the first symbol in current file format)
                        }
                        return new SimilarityMatrix(new StringReader(sb.ToString()));
                    }
                }
            }
            else
            {
                throw new Exception("Unsupported molecule type: " + moleculeType);
            }
        }
    }
}