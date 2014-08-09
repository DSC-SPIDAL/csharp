using System.ComponentModel;

namespace Salsa.Core.Configuration.Sections
{
    public class DAVectorSpongeSection : Section
    {
        #region Constructor

        public DAVectorSpongeSection()
        {
            UseSponge = false;
            SpongeFactor1 = 3.0;
            SpongeFactor2 = 3.0;
            SpongePOption = 1;
            SpongePWeight = 0.1;
            CreateSpongeScaledSquaredWidth = -1.0;
            ContinuousClustering = true;
            ParameterVectorDimension = 2;

            SpongeTemperature1 = -1.0;
            SpongeTemperature2 = -1.0;
            RestartTemperature = -1.0;

            NumberDataPoints = -1;
            SelectedInputLabel = 6;
            OutputFileType = 0;
            Replicate = 1;

            SigmaMethod = 0;
            FinalTargetTemperature = 3.0;
            FinalTargetSigma0 = 0.0;
            InitialSigma0 = 0.0;

            ClusterCountOutput = 0;
            NumberNearbyClusters = 5;
            NearbySpongePointLimit = -1.0;

            ProcessingOption = 0;

            CacheLineSize = 0;
            ClusterPrintNumber = 5;
            PrintInterval = 3;
            RemovalDiagnosticPrint = false;
            MagicTemperatures = new[] {4.0, 3.0, 2.0, 1.0, 0.5};
            MagicIndex = 0;

            MaxNcentPerNode = 0;
            MaxNcentTotal = 0;
            TargetNcentPerPoint = 20;
            TargetMinimumNcentPerPoint = 1;
            MaxNcentPerPoint = 25;

            MaxIntegerComponents = 2;
            MaxDoubleComponents = 3;
            MaxMPITransportBuffer = 500;
            MaxNumberAccumulationsPerNode = 30000;
            MaxTransportedClusterStorage = 500;

            ExpArgumentCut1 = 20.0;
            ExpArgumentCut2 = 40.0;
            ExpArgumentCut3 = 50.0;
            Tminimum = -1000.0;

            InitalNcent = 1;
            MinimumCountForClusterCk = 1.0;
            MinimumCountForClusterCkWithSponge = 1.5;
            MinimuCountForClusterPoints = 2;
            CountForClusterCkToBeZero = 0.001;
            AddSpongeScaledWidthSquared = -1.0;

            InitialCoolingFactor = 0.9;
            FineCoolingFactor = 0.99;
            WaitIterations = 10;

            IterationAtEnd = 2000;
            ConvergenceLoopLimit = 20;

            FreezingLimit = 0.002;
            MalphaMaxChange = 0.005;
            MaxNumberSplitClusters = 3;
            ConvergeIntermediateClusters = false;
            TooSmallToSplit = 4.0;
            ScaledWidthSquaredToSplit = 6.0;

            ClusterLimitForDistribution = -1;
            TemperatureLimitForDistribution = -1.0;

            DebugPrintOption = 1;
            ConsoleDebugOutput = true;
        }

        #endregion

        #region Properties

        #region Parameters

        [Category("Parameters")]
        public bool UseSponge { get; set; }

        [Category("Parameters")]
        public double SpongeFactor1 { get; set; }

        [Category("Parameters")]
        public double SpongeFactor2 { get; set; }

        [Category("Parameters")]
        public int SpongePOption { get; set; }

        [Category("Parameters")]
        public double SpongePWeight { get; set; }

        [Category("Parameters")]
        public double CreateSpongeScaledSquaredWidth { get; set; }

        [Category("Parameters")]
        public bool ContinuousClustering { get; set; }

        [Category("Parameters")]
        public int ParameterVectorDimension { get; set; }

        [Category("Parameters")]
        public double SpongeTemperature1 { get; set; }

        [Category("Parameters")]
        public double SpongeTemperature2 { get; set; }

        [Category("Parameters")]
        public double RestartTemperature { get; set; }

        [Category("Parameters")]
        public int NumberDataPoints { get; set; }

        [Category("Parameters")]
        public int SelectedInputLabel { get; set; }

        [Category("Parameters")]
        public int OutputFileType { get; set; }

        [Category("Parameters")]
        public int Replicate { get; set; }

        [Category("Parameters")]
        public int SigmaMethod { get; set; }

        [Category("Parameters")]
        public double FinalTargetTemperature { get; set; }

        [Category("Parameters")]
        public double FinalTargetSigma0 { get; set; }

        [Category("Parameters")]
        public double InitialSigma0 { get; set; }

        [Category("Parameters")]
        public int ClusterCountOutput { get; set; }

        [Category("Parameters")]
        public int NumberNearbyClusters { get; set; }

        [Category("Parameters")]
        public double NearbySpongePointLimit { get; set; }


        [Category("Parameters")]
        public int ProcessingOption { get; set; }

        [Category("Parameters")]
        public int CacheLineSize { get; set; }

        [Category("Parameters")]
        public int ClusterPrintNumber { get; set; }

        [Category("Parameters")]
        public int PrintInterval { get; set; }

        [Category("Parameters")]
        public bool RemovalDiagnosticPrint { get; set; }

        [Category("Parameters")]
        public double[] MagicTemperatures { get; set; }

        [Category("Parameters")]
        public int MagicIndex { get; set; }

        [Category("Parameters")]
        public int MaxNcentPerNode { get; set; }

        [Category("Parameters")]
        public int MaxNcentTotal { get; set; }

        [Category("Parameters")]
        public int TargetNcentPerPoint { get; set; }

        [Category("Parameters")]
        public int TargetMinimumNcentPerPoint { get; set; }

        [Category("Parameters")]
        public int MaxNcentPerPoint { get; set; }

        [Category("Parameters")]
        public int MaxIntegerComponents { get; set; }

        [Category("Parameters")]
        public int MaxDoubleComponents { get; set; }

        [Category("Parameters")]
        public int MaxMPITransportBuffer { get; set; }

        [Category("Parameters")]
        public int MaxNumberAccumulationsPerNode { get; set; }

        [Category("Parameters")]
        public int MaxTransportedClusterStorage { get; set; }

        [Category("Parameters")]
        public double ExpArgumentCut1 { get; set; }

        [Category("Parameters")]
        public double ExpArgumentCut2 { get; set; }

        [Category("Parameters")]
        public double ExpArgumentCut3 { get; set; }

        [Category("Parameters")]
        public double Tminimum { get; set; }

        [Category("Parameters")]
        public int InitalNcent { get; set; }

        [Category("Parameters")]
        public double MinimumCountForClusterCk { get; set; }

        [Category("Parameters")]
        public double MinimumCountForClusterCkWithSponge { get; set; }

        [Category("Parameters")]
        public int MinimuCountForClusterPoints { get; set; }

        [Category("Parameters")]
        public double CountForClusterCkToBeZero { get; set; }

        [Category("Parameters")]
        public double AddSpongeScaledWidthSquared { get; set; }

        [Category("Parameters")]
        public double InitialCoolingFactor { get; set; }

        [Category("Parameters")]
        public double FineCoolingFactor { get; set; }

        [Category("Parameters")]
        public int WaitIterations { get; set; }

        [Category("Parameters")]
        public int IterationAtEnd { get; set; }

        [Category("Parameters")]
        public int ConvergenceLoopLimit { get; set; }

        [Category("Parameters")]
        public double FreezingLimit { get; set; }

        [Category("Parameters")]
        public double MalphaMaxChange { get; set; }

        [Category("Parameters")]
        public int MaxNumberSplitClusters { get; set; }

        [Category("Parameters")]
        public bool ConvergeIntermediateClusters { get; set; }

        [Category("Parameters")]
        public double TooSmallToSplit { get; set; }

        [Category("Parameters")]
        public double ScaledWidthSquaredToSplit { get; set; }

        [Category("Parameters")]
        public int ClusterLimitForDistribution { get; set; }

        [Category("Parameters")]
        public double TemperatureLimitForDistribution { get; set; }

        [Category("Parameters")]
        public string Pattern { get; set; }

        [Category("Parameters")]
        public int NodeCount { get; set; }

        [Category("Parameters")]
        public int ThreadCount { get; set; }

        [Category("Parameters")]
        public int MPIPerNodeCount { get; set; }

        #endregion

        #region Debug

        [Category("Debug")]
        public int DebugPrintOption { get; set; }

        [Category("Debug")]
        public bool ConsoleDebugOutput { get; set; }

        #endregion

        #region I/O

        [Category("I/O")]
        public string ClusterFile { get; set; }

        [Category("I/O")]
        public string DistanceMatrixFile { get; set; }

        [Category("I/O")]
        public string LabelFile { get; set; }

        [Category("I/O")]
        public string TimingFile { get; set; }

        [Category("I/O")]
        public string SummaryFile { get; set; }

        #endregion

        #endregion

        internal override void ExpandMacro(MacroReplacement macroReplacement)
        {
            // Empty;
        }
    }
}