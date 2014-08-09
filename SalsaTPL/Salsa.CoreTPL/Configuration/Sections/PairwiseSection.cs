using System.ComponentModel;

namespace Salsa.Core.Configuration.Sections
{
    public class PairwiseSection : Section
    {
        #region Constructors

        #endregion

        #region Properties

        [Category("I/O")]
        [Description("Full Path to the location where the index file is to be read from. Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\index.txt")]
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
        [Description("Full Path to the location where the data file is to be read from. Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\distance.txt")]
        public string DistanceMatrixFile
        {
            get { return _dataFileName; }
            set { _dataFileName = value; }
        }

        [Category("I/O")]
        [Description("Full Path to the location where the timing file is to be written to.  Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\pwc_timing.txt")]
        public string TimingFile
        {
            get { return _timingOutputFileName; }
            set { _timingOutputFileName = value; }
        }

        [Category("I/O")]
        [Description("Full Path to the location where the summary file is to be written to.  Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\pwc_summary.txt")]
        public string SummaryFile
        {
            get { return _summaryFileName; }
            set { _summaryFileName = value; }
        }

        [Category("I/O")]
        [Description("Full Path to the location where the plot file is to be written to.  Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\pwc_plot.pviz")]
        public string CenterPlotFile { get; set; }

        [Category("I/O")]
        [Description("Full Path to the location where the cluster file is to be written to. Macros will be expanded.")]
        [DefaultValue(@"$(OutputRootPath)\cluster.txt")]
        public string ClusterFile
        {
            get { return _clusterFile; }
            set { _clusterFile = value; }
        }

        [Category("I/O")]
        public bool ReadPartialMatrix { get; set; }

        [Browsable(false)]
        public string Comments { get; set; }

        [Category("I/O")]
        [Description("point-num<TAB>name<TAB>length")]
        public string LabelFile { get; set; }

        [Category("I/O")]
        [Description("point-num<TAB>x<TAB>y<TAB>z<TAB>cluster-num")]
        public string AddMdsFile { get; set; }

        [Category("I/O")]
        [Description("Previous cluster assignmnet of point-num<TAB>clus-num")]
        public string ClusterNumberFile { get; set; }

        [Category("Parameters")]
        [Description("The number of data points or the dimension of the square distance matrix.")]
        public int DataPoints { get; set; }


        [Category("Parameters")]
        [Description("Description Needed")]
        public int ProcessingOption { get; set; }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int TransformDimension
        {
            get { return _transformDimension; }
            set { _transformDimension = value; }
        }

        [Category("Parameters")]
        [Description("The maximum number of center to consider when clustering.")]
        public int MaxNcent
        {
            get { return _maxNcent; }
            set { _maxNcent = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int Splitorexpandit
        {
            get { return _splitorexpandit; }
            set { _splitorexpandit = value; }
        }

        [Browsable(false)]
        public string Pattern
        {
            get { return _pattern; }
            set { _pattern = value; }
        }

        [Browsable(false)]
        public int ThreadCount
        {
            get { return _threadCount; }
            set { _threadCount = value; }
        }

        [Browsable(false)]
        public int NodeCount
        {
            get { return _nodeCount; }
            set { _nodeCount = value; }
        }

        [Browsable(false)]
        public int MPIperNodeCount
        {
            get { return _mPIperNodeCount; }
            set { _mPIperNodeCount = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int MPIIOStrategy { get; set; }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double ToosmalltoSplit
        {
            get { return _toosmalltoSplit; }
            set { _toosmalltoSplit = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double MinEigtest
        {
            get { return _minEigtest; }
            set { _minEigtest = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public bool ConvergeIntermediateClusters { get; set; }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int Waititerations
        {
            get { return _waititerations; }
            set { _waititerations = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double Epsi_max_change
        {
            get { return _epsi_max_change; }
            set { _epsi_max_change = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double InitialCoolingFactor
        {
            get { return _initialCoolingFactor; }
            set { _initialCoolingFactor = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double FineCoolingFactor
        {
            get { return _fineCoolingFactor; }
            set { _fineCoolingFactor = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double Eigenvaluechange
        {
            get { return _eigenvaluechange; }
            set { _eigenvaluechange = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double Eigenvectorchange
        {
            get { return _eigenvectorchange; }
            set { _eigenvectorchange = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int Iterationatend
        {
            get { return _iterationatend; }
            set { _iterationatend = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int ConvergenceLoopLimit
        {
            get { return _convergenceLoopLimit; }
            set { _convergenceLoopLimit = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public double FreezingLimit
        {
            get { return _freezingLimit; }
            set { _freezingLimit = value; }
        }

        [Category("Parameters")]
        [Description("Description Needed")]
        public int PowerIterationLimit
        {
            get { return _powerIterationLimit; }
            set { _powerIterationLimit = value; }
        }

        [Category("Debug")]
        [Description("0 = None, 1 = Full, 2 = Summary")]
        [DefaultValue(1)]
        public int DebugPrintOption
        {
            get { return _debugPrintOption; }
            set { _debugPrintOption = value; }
        }

        [Category("Debug")]
        [DefaultValue(true)]
        public bool ConsoleDebugOutput
        {
            get { return _consoleDebugOutput; }
            set { _consoleDebugOutput = value; }
        }

        [Category("Parameters")]
        [Description("If true use the Ken Rose Continuous Clustering")]
        public bool ContinuousClustering { get; set; }

        [Category("Parameters")]
        [Description("Specify MDS versions of center finding; = 0 ignore")]
        public int AddMds
        {
            get { return _addMDS; }
            set { _addMDS = value; }
        }

        [Category("Parameters")]
        [Description("Comma separated list of bucket fractions")]
        public string BucketFractions
        {
            get { return _bucketFractions; }
            set { _bucketFractions = value; }
        }

        [Category("Parameters")]
        [Description("Number of centers to be found with each method")]
        public int NumberOfCenters
        {
            get { return _numberOfCenters; }
            set { _numberOfCenters = value; }
        }

        [Category("Parameters")]
        [Description("Number of centers to include in each center type of the output plot")]
        public int CenterPointsPerCenterTypeInOuput
        {
            get { return _centerPointsPerCenterTypeInOuput; }
            set { _centerPointsPerCenterTypeInOuput = value; }
        }

        #endregion

        #region Members

        private int _addMDS = 1;
        private string _bucketFractions = "0.15,0.4,0.75";
        private int _centerPointsPerCenterTypeInOuput = 3;
        private string _clusterFile = @"$(OutputRootPath)\cluster.txt";
        private bool _consoleDebugOutput = true;
        private int _convergenceLoopLimit = 2000;
        private string _dataFileName = @"distance.bin";
        private int _debugPrintOption = 1;
        private double _eigenvaluechange = 0.001;
        private double _eigenvectorchange = 0.001;
        private double _epsi_max_change = 0.001;
        private double _fineCoolingFactor = 0.99;
        private double _freezingLimit = 0.002;
        private string _indexFile = @"$(OutputRootPath)\index.txt";
        private double _initialCoolingFactor = 0.9;
        private int _iterationatend = 2000;
        private int _mPIperNodeCount = 24;
        private int _maxNcent = 10;
        private double _minEigtest = -0.01;
        private int _nodeCount = 32;
        private int _numberOfCenters = 8;
        private string _pattern = string.Empty;
        private int _powerIterationLimit = 200;
        private int _splitorexpandit = 1;
        private string _summaryFileName = @"$(OutputRootPath)\pwc_summary.txt";
        private int _threadCount = 1;
        private string _timingOutputFileName = @"$(OutputRootPath)\pwc_timing.txt";
        private double _toosmalltoSplit = 50;
        private int _transformDimension = 4;
        private int _waititerations = 10;

        #endregion

        internal override void ExpandMacro(MacroReplacement macroReplacement)
        {
            base.ExpandMacro(macroReplacement);
            DistanceMatrixFile = macroReplacement.Expand(DistanceMatrixFile);
            TimingFile = macroReplacement.Expand(TimingFile);
            SummaryFile = macroReplacement.Expand(SummaryFile);
            ClusterFile = macroReplacement.Expand(ClusterFile);
        }
    }
}