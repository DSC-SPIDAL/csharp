using System.ComponentModel;
using Salsa.Core.Bio.Algorithms;

namespace Salsa.Core.Configuration.Sections
{
    public class ScoringMatrixTypeConverter : StringConverter
    {
        public override bool GetStandardValuesSupported(ITypeDescriptorContext context)
        {
            return true;
        }

        public override bool GetStandardValuesExclusive(ITypeDescriptorContext context)
        {
            return true;
        }

        public override StandardValuesCollection GetStandardValues(ITypeDescriptorContext context)
        {
            return new StandardValuesCollection(ScoringMatrix.MatrixNames);
        }
    }

    public class SmithWatermanSection : Section
    {
        #region Constructors

        #endregion

        #region Members

        private AlignmentType _alignmentType = AlignmentType.Nucleic;
        private string _distanceFile = @"distance.bin";
        private DistanceFunctionType _distanceFunction = DistanceFunctionType.PercentIdentity;
        private string _fastaFile = "$(InputRootPath)";
        private float _gapExtension = 4;
        private float _gapOpen = 14;
        private string _indexFile = @"$(OutputRootPath)\index.txt";
        private int _nodeCount;
        private int _processPerNodeCount;
        private string _scoringMatrixName = "EDNAFULL";
        private int _sequenceCount;
        private string _summaryFile = @"$(OutputRootPath)\swg_summary.txt";
        private string _timingFile = @"$(OutputRootPath)\swg_timing.txt";
        private string _writeAlignmentFile = @"$(OutputRootPath)\swg_alignments.txt";
        private bool _writeFullMatrix = true;
        private bool _writePartialMatrix = true;

        #endregion

        #region Properties

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

        [Description("The type of sequences in the fasta file")]
        public AlignmentType AlignmentType
        {
            get { return _alignmentType; }
            set
            {
                _alignmentType = value;

                if (_alignmentType == AlignmentType.Nucleic)
                {
                    ScoringMatrixName = "EDNAFULL";
                }
                else
                {
                    _distanceFunction = DistanceFunctionType.PercentIdentity;
                }
            }
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
            "If true the partial matrix for each process is written to file.  Good idea to set this to true for long running jobs."
            )]
        public bool WritePartialMatrix
        {
            get { return _writePartialMatrix; }
            set
            {
                _writePartialMatrix = value;
                OnPropertyChanged("WritePartialMatrix");
            }
        }

        [Description(
            "If true the alignment for each pair is written to the file specified in the WriteAlignmentsPath property.")
        ]
        public bool WriteAlignments { get; set; }

        [Description("The cost associated with opening a gap in the alignment.")]
        public float GapOpenPenalty
        {
            get { return _gapOpen; }
            set
            {
                _gapOpen = value;
                OnPropertyChanged("GapOpen");
            }
        }

        [Description("The cost associated with extending a gap in the alignment.")]
        public float GapExtensionPenalty
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
            WriteAlignmentsFile = macroReplacement.Expand(WriteAlignmentsFile);
        }
    }
}