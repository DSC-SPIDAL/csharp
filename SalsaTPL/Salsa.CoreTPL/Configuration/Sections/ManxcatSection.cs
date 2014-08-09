using System.ComponentModel;

namespace Salsa.Core.Configuration.Sections
{
    public class ManxcatSection : Section
    {
        #region Constructors

        #endregion

        #region Properties

        #region I/O Propertis

        [Category("I/O")]
        public string BaseResultDirectoryName
        {
            get { return _baseResultDirectoryName; }
            set { _baseResultDirectoryName = value; }
        }

        [Category("I/O")]
        public string ControlDirectoryName
        {
            get { return _controlDirectoryName; }
            set { _controlDirectoryName = value; }
        }

        [Category("I/O")]
        public string ClusterDirectory
        {
            get { return _clusterDirectory; }
            set { _clusterDirectory = value; }
        }

        [Category("I/O")]
        public string DistanceMatrixFile
        {
            get { return _dataFileName; }
            set { _dataFileName = value; }
        }

        [Category("I/O")]
        public string DataLabelsFileName
        {
            get { return _dataLabelsFileName; }
            set { _dataLabelsFileName = value; }
        }

        [Category("I/O")]
        public string ReducedVectorOutputFileName
        {
            get { return _reducedVectorOutputFileName; }
            set { _reducedVectorOutputFileName = value; }
        }

        [Category("I/O")]
        public string ResultDirectoryExtension
        {
            get { return _resultDirectoryExtension; }
            set { _resultDirectoryExtension = value; }
        }

        [Category("I/O")]
        public string TimingOutputFileName
        {
            get { return _timingOutputFileName; }
            set { _timingOutputFileName = value; }
        }

        [Category("I/O")]
        public string SummaryOutputFileName
        {
            get { return _summaryOutputFileName; }
            set { _summaryOutputFileName = value; }
        }

        [Category("I/O")]
        public string InitializationFileName
        {
            get { return _initializationFileName; }
            set { _initializationFileName = value; }
        }

        [Category("I/O")]
        public string WeightingFileName
        {
            get { return _weightingFileName; }
            set { _weightingFileName = value; }
        }

        [Category("I/O")]
        public string ScalingFileName
        {
            get { return _scalingFileName; }
            set { _scalingFileName = value; }
        }

        [Category("I/O")]
        public string Selectedvariedpointfile
        {
            get { return _selectedvariedpointfile; }
            set { _selectedvariedpointfile = value; }
        }

        [Category("I/O")]
        public string Selectedfixedpointfile
        {
            get { return _selectedfixedpointfile; }
            set { _selectedfixedpointfile = value; }
        }

        [Category("I/O")]
        public string RotationLabelsFileName
        {
            get { return _rotationLabelsFileName; }
            set { _rotationLabelsFileName = value; }
        }

        [Category("I/O")]
        public bool ReadPartialMatrix
        {
            get { return _readPartialMatrix; }
            set { _readPartialMatrix = value; }
        }


        [Category("I/O")]
        [Description("Frequency, in iterations, with which to write out MDS coordinates.")]
        [DefaultValue(80)]
        public int CoordinateWriteFrequency { get; set; }


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

        #endregion

        public int DataPoints { get; set; }

        public bool CalcFixedCrossFixed
        {
            get { return _calcFixedCrossFixed; }
            set { _calcFixedCrossFixed = value; }
        }

        public int StoredDistanceOption
        {
            get { return _storedDistanceOption; }
            set { _storedDistanceOption = value; }
        }

        public int DiskDistanceOption
        {
            get { return _diskDistanceOption; }
            set { _diskDistanceOption = value; }
        }

        public int ProcessingOption { get; set; }

        public int DistanceProcessingOption { get; set; }

        public int InitializationOption { get; set; }

        public int WeightingOption { get; set; }

        public bool Write2Das3D
        {
            get { return _write2Das3D; }
            set { _write2Das3D = value; }
        }

        public int LocalVectorDimension
        {
            get { return _localVectorDimension; }
            set { _localVectorDimension = value; }
        }

        public string Selectedvariedpoints
        {
            get { return _selectedvariedpoints; }
            set { _selectedvariedpoints = value; }
        }

        public string VariedPointCriterion
        {
            get { return _variedPointCriterion; }
            set { _variedPointCriterion = value; }
        }

        public string Selectedfixedpoints
        {
            get { return _selectedfixedpoints; }
            set { _selectedfixedpoints = value; }
        }

        public string FixedPointCriterion
        {
            get { return _fixedPointCriterion; }
            set { _fixedPointCriterion = value; }
        }

        public string ConversionOption
        {
            get { return _conversionOption; }
            set { _conversionOption = value; }
        }

        public string ConversionInformation
        {
            get { return _conversionInformation; }
            set { _conversionInformation = value; }
        }

        public int RotationOption { get; set; }

        public int InitializationLoops
        {
            get { return _initializationLoops; }
            set { _initializationLoops = value; }
        }

        public int Chisqnorm
        {
            get { return _chisqnorm; }
            set { _chisqnorm = value; }
        }

        public int DistanceFormula
        {
            get { return _distanceFormula; }
            set { _distanceFormula = value; }
        }

        public int FullSecondDerivativeOption { get; set; }

        public double MinimumDistance
        {
            get { return _minimumDistance; }
            set { _minimumDistance = value; }
        }

        public int FunctionErrorCalcMultiplier
        {
            get { return _functionErrorCalcMultiplier; }
            set { _functionErrorCalcMultiplier = value; }
        }

        public int ChisqPrintConstant
        {
            get { return _chisqPrintConstant; }
            set { _chisqPrintConstant = value; }
        }

        public int Maxit
        {
            get { return _maxit; }
            set { _maxit = value; }
        }

        public int Nbadgo
        {
            get { return _nbadgo; }
            set { _nbadgo = value; }
        }

        public double ChisqChangePerPoint
        {
            get { return _chisqChangePerPoint; }
            set { _chisqChangePerPoint = value; }
        }

        public double FletcherRho
        {
            get { return _fletcherRho; }
            set { _fletcherRho = value; }
        }

        public double FletcherSigma
        {
            get { return _fletcherSigma; }
            set { _fletcherSigma = value; }
        }

        public double Omega
        {
            get { return _omega; }
            set { _omega = value; }
        }

        public int OmegaOption { get; set; }

        public double QHighInitialFactor
        {
            get { return _qHighInitialFactor; }
            set { _qHighInitialFactor = value; }
        }

        public double QgoodReductionFactor
        {
            get { return _qgoodReductionFactor; }
            set { _qgoodReductionFactor = value; }
        }

        public int QLimitscalculationInterval
        {
            get { return _qLimitscalculationInterval; }
            set { _qLimitscalculationInterval = value; }
        }

        public double Extraprecision
        {
            get { return _extraprecision; }
            set { _extraprecision = value; }
        }

        public int AddonforQcomputation
        {
            get { return _addonforQcomputation; }
            set { _addonforQcomputation = value; }
        }

        public int InitialSteepestDescents { get; set; }

        public int TimeCutmillisec
        {
            get { return _timeCutmillisec; }
            set { _timeCutmillisec = value; }
        }

        public double CGResidualLimit
        {
            get { return _cGResidualLimit; }
            set { _cGResidualLimit = value; }
        }

        public int PowerIterationLimit
        {
            get { return _powerIterationLimit; }
            set { _powerIterationLimit = value; }
        }

        public double Eigenvaluechange
        {
            get { return _eigenvaluechange; }
            set { _eigenvaluechange = value; }
        }

        public double Eigenvectorchange
        {
            get { return _eigenvectorchange; }
            set { _eigenvectorchange = value; }
        }

        public bool Derivtest { get; set; }

        public int RunNumber
        {
            get { return _runNumber; }
            set { _runNumber = value; }
        }

        public string RunSetLabel
        {
            get { return _runSetLabel; }
            set { _runSetLabel = value; }
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

        public int MPIIOStrategy { get; set; }

        public int HistogramBinCount
        {
            get { return _histogramBinCount; }
            set { _histogramBinCount = value; }
        }

        public string Extradata1
        {
            get { return _extradata1; }
            set { _extradata1 = value; }
        }

        public string Extradata2
        {
            get { return _extradata2; }
            set { _extradata2 = value; }
        }

        public string Extradata3
        {
            get { return _extradata3; }
            set { _extradata3 = value; }
        }

        public string Extradata4
        {
            get { return _extradata4; }
            set { _extradata4 = value; }
        }

        public int ExtraOption1 { get; set; }

        [Category("Debug")]
        public int DebugPrintOption
        {
            get { return _debugPrintOption; }
            set { _debugPrintOption = value; }
        }

        [Category("Debug")]
        public bool ConsoleDebugOutput
        {
            get { return _consoleDebugOutput; }
            set { _consoleDebugOutput = value; }
        }

        public string Comment
        {
            get { return _comment; }
            set { _comment = value; }
        }

        public double UndefinedDistanceValue
        {
            get { return _undefinedDistanceValue; }
            set { _undefinedDistanceValue = value; }
        }

        public string DistanceWeightsCuts
        {
            get { return _distanceWeightsCuts; }
            set { _distanceWeightsCuts = value; }
        }

        public float DistanceCut
        {
            get { return _distanceCut; }
            set { _distanceCut = value; }
        }

        public int LinkCut
        {
            get { return _linkCut; }
            set { _linkCut = value; }
        }

        public int TransformMethod { get; set; }


        public float TransformParameter
        {
            get { return _transformParameter; }
            set { _transformParameter = value; }
        }

        [Category("Description")]
        public string ManxcatRunDescription
        {
            get { return _manxcatRunDescription; }
            set { _manxcatRunDescription = value; }
        }

        [Category("Description")]
        public string ManxcatRunName
        {
            get { return _manxcatRunName; }
            set { _manxcatRunName = value; }
        }

        [Category("Density")]
        public string SelectedClusters { get; set; }

        [Category("Density")]
        public string ClusterFile { get; set; }

        [Category("Density")]
        public double Pcutf
        {
            get { return _pcutf; }
            set { _pcutf = value; }
        }

        [Category("Density")]
        public double Alpha
        {
            get { return _alpha; }
            set { _alpha = value; }
        }

        [Category("Density")]
        public int Yres
        {
            get { return _yres; }
            set { _yres = value; }
        }

        [Category("Density")]
        public int Xres
        {
            get { return _xres; }
            set { _xres = value; }
        }

        [Category("Density")]
        public double YmaxBound
        {
            get { return _ymaxBound; }
            set { _ymaxBound = value; }
        }

        [Category("Density")]
        public double XmaxBound
        {
            get { return _xmaxBound; }
            set { _xmaxBound = value; }
        }

        [Category("Density")]
        public bool Normalize
        {
            get { return _normalize; }
            set { _normalize = value; }
        }

        [Category("Publish Settings")]
        public string ServerUrlPrefix
        {
            get { return _serverUrlPrefix; }
            set { _serverUrlPrefix = value; }
        }

        #endregion

        #region Members

        private int _addonforQcomputation = 2;
        private double _alpha = 2;
        private string _baseResultDirectoryName = "$(OutputRootPath)";
        private double _cGResidualLimit = 1E-05;
        private bool _calcFixedCrossFixed = true;
        private double _chisqChangePerPoint = 0.001;
        private int _chisqPrintConstant = 1;
        private int _chisqnorm = 2;
        private string _clusterDirectory = string.Empty;
        private string _comment = string.Empty;
        private bool _consoleDebugOutput = true;
        private string _controlDirectoryName = "$(ConfigRootPath)";
        private string _conversionInformation = string.Empty;
        private string _conversionOption = string.Empty;
        private string _dataFileName = @"$(OutputRootPath)\distance.txt";
        private string _dataLabelsFileName = string.Empty;
        private int _debugPrintOption = 2;
        private int _diskDistanceOption = 2;
        private float _distanceCut = -1.0f; // -1 to indicate no cut on distance
        private int _distanceFormula = 1;

        private string _distanceWeightsCuts = string.Empty;
                       // List of Distance Cuts for Weights. These define upper limits of distance bins with "infinity" as upper and 0 of course as lower

        private double _eigenvaluechange = 0.001;
        private double _eigenvectorchange = 0.001;
        private string _extradata1 = string.Empty;
        private string _extradata2 = string.Empty;
        private string _extradata3 = string.Empty;
        private string _extradata4 = string.Empty;
        private double _extraprecision = 0.05;
        private string _fixedPointCriterion = "none";
        private double _fletcherRho = 0.25;
        private double _fletcherSigma = 0.75;

        private int _functionErrorCalcMultiplier = 10;
        private int _histogramBinCount = 100;
        private string _indexFile = @"$(OutputRootPath)\index.txt";
        private string _initializationFileName = string.Empty;
        private int _initializationLoops = 4;
        private int _linkCut = 5; // Delete from fit all points with <= LinkCut connections
        private int _localVectorDimension = 3;
        private int _mPIperNodeCount = 1;
        private string _manxcatRunDescription = "Unspecified Description";
        private string _manxcatRunName = "Unspecified Run";
        private int _maxit = 80;
        private double _minimumDistance = -0.001;
        private int _nbadgo = 6;
        private int _nodeCount = 30;
        private bool _normalize = true;
        private double _omega = 1.25;
        private string _pattern = string.Empty;
        private double _pcutf = 0.85;
        private int _powerIterationLimit = 200;
        private double _qHighInitialFactor = 0.01;
        private int _qLimitscalculationInterval = 1;
        private double _qgoodReductionFactor = 0.5;
        private bool _readPartialMatrix = true;
        private string _reducedVectorOutputFileName = @"$(OutputRootPath)\points.txt";
        private string _resultDirectoryExtension = string.Empty;
        private string _rotationLabelsFileName = string.Empty;
        private int _runNumber = 27;
        private string _runSetLabel = string.Empty;
        private string _scalingFileName = string.Empty;
        private string _selectedfixedpointfile = string.Empty;
        private string _selectedfixedpoints = string.Empty;
        private string _selectedvariedpointfile = string.Empty;
        private string _selectedvariedpoints = string.Empty;
        private string _serverUrlPrefix = "http://salsahpc.indiana.edu/manxcat";
        private int _storedDistanceOption = 2;
        private string _summaryOutputFileName = @"$(OutputRootPath)\mds_summary.txt";
        private int _threadCount = 24;
        private int _timeCutmillisec = -1;
        private string _timingOutputFileName = @"$(OutputRootPath)\mds_timings.txt";
        private float _transformParameter = 0.125f; // Ignored unless _transformParameter = 11

        private double _undefinedDistanceValue = -1.0;
                       // If positive replace undefined distances by this value; use -1.0 if want to unset these distances

        private string _variedPointCriterion = "all";

        private string _weightingFileName = string.Empty;
                       // the format of the file is row<tab>rowcount<tab>col<tab>colcount<tab>scalingfactor

        private bool _write2Das3D = true;

        /* Density and Web page creation members */
        private double _xmaxBound = 1;
        private int _xres = 50;
        private double _ymaxBound = 1;
        private int _yres = 50;

        #endregion

        internal override void ExpandMacro(MacroReplacement macroReplacement)
        {
            base.ExpandMacro(macroReplacement);
            BaseResultDirectoryName = macroReplacement.Expand(BaseResultDirectoryName);
            ControlDirectoryName = macroReplacement.Expand(ControlDirectoryName);
            ClusterDirectory = macroReplacement.Expand(ClusterDirectory);
            DistanceMatrixFile = macroReplacement.Expand(DistanceMatrixFile);
            DataLabelsFileName = macroReplacement.Expand(DataLabelsFileName);
            ReducedVectorOutputFileName = macroReplacement.Expand(ReducedVectorOutputFileName);
            TimingOutputFileName = macroReplacement.Expand(TimingOutputFileName);
            InitializationFileName = macroReplacement.Expand(InitializationFileName);
            WeightingFileName = macroReplacement.Expand(WeightingFileName);
            Selectedvariedpointfile = macroReplacement.Expand(Selectedvariedpointfile);
            Selectedfixedpointfile = macroReplacement.Expand(Selectedfixedpointfile);
            RotationLabelsFileName = macroReplacement.Expand(RotationLabelsFileName);
            SummaryOutputFileName = macroReplacement.Expand(SummaryOutputFileName);
        }
    }
}