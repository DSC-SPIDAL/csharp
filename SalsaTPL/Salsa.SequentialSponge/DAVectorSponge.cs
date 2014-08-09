using System;
using System.Collections;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using HPC.Utilities;
using MPI;
using Salsa.Core.Configuration;
using Salsa.Core.Configuration.Sections;
using Salsa.Core.Blas;
using Salsa.Core;

using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    class Program
    {
        #region Overall parameters

        //  Specify Clustering Algorithm
        public static bool DoKmeans = false;    //  If true use simple Kmeans
        public static bool DoLCMS = false; // If true solve 2D LCMS problem
        public static bool DoBasicDA = false;   // If true do basic Vector DA
        public static bool FollowUpKmeans = false;  // If true follow up BasicDA with Kmeans

        //  Specify Nature of Vector Data
        public static int ParameterVectorDimension = 2;        // Vector Dimension of Points

        //  Set Fixed errors in each dimension
        public static int SigmaMethod = 0; //   if 0, sigma of clusters unspecified. Otherwise specifies calculation method
        public static double[] SigmaVectorParameters_i_; // Control calculation of Cluster Sigma with dimension ParameterVectorDimension
        public static double FinalTargetTemperature = 3.0;  // Used in SigmaMethod =3
        public static double FinalTargetSigma0 = 0.0; // Used in SigmaMethod =3
        public static double InitialSigma0 = 0.0;   // Initial value

        //  Set Dataset properties
        public static int NumberDataPoints = -1;    // if -1 then find number of points from input
        public static int SelectedInputLabel = 6;   // Charge to select; if negative veto this label
        public static int RestartSelectedInputLabel = 6;   // Charge to select in Restart; if negative veto this label
        public static int InputFileType = 0;    // If 0 raw data; =1 Output file
        public static int RestartInputFileType = 1;    // If 0 raw data; =1 Output file for Restart
        public static bool Refinement = true;       // If False do no refinement
        public static int StartPointPositiononInputLine = 1;    // Starting position of Point Data on Input File
        public static int ClusterIndexonInputLine = -1;    // Position of Cluster Index on Input File for K means
        public static int FirstClusterValue = 1;    // Index for first cluster in Kmeans initialization
        public static double[][] PointPosition;  // Position of Point x is PointPosition[x, i] where i runs over vector space directions 0 ... ParameterVectorDimension-1
        public static int[] PointOriginalIndex; // Point Index on File
        public static int[] PointLabel; // Label on File
        public static int Replicate = 1;    // Replicate points for scaling tests
        public static double[][] FullPoint3DPosition;   // 3D Position of Points -- all stored
        public static int RW3DData = -1; // If positive read and write Plotviz files of this dimension

        //  Control Execution Model in general
        public static string ControlFileName = "";  // Control File Name
        public static int ProcessingOption = 0; // Processing Option
        public static int cachelinesize = 0;    // Specify Cache line minimum size
        public static bool CalculateCorrelationMatrix = false;  // If True calculate correlations
        public static bool CalculateIndividualWidths = false;  // If True calculate individual widths
        public static bool CalculateEigenvaluesfromMatrix = false;  // If false use iteration

        //  Describe mode for parallelism over clusters (as opposed to normal parallelism over points)
        public static int ClusterLimitforDistribution = -1;  // Go into distributed mode if Number of Clusters per node greater than this
        public static int ActualClusterNumberforDistribution = -1;  // Number of Clusters per node when go into distributed mode
        public static double TemperatureLimitforDistribution = -1.0;   // Go into distributed mode if Temperature less than this
        public static double ActualTemperatureforDistribution = -1.0;   // Temperature when go into distributed mode 
        public static double TimeatDistribution = 0.0;  // Time when distributed execution started

        public static int MaxIntegerComponents = 2; // Maximum Number of Integer Components in Transport (for MPI arrays)
        public static int MaxDoubleComponents = 3;  // Maximum Number of Double Components in Transport
        public static int MaxMPITransportBuffer = 500; // Maximum Number of Clusters in MPI move routines
        public static int ActualMaxMPITransportBuffer = 0; // Needed Maximum Number of Clusters in MPI move routines
        public static int MaxNumberAccumulationsperNode = 30000; // Maximum Number of Clusters for a single node in Distributed Mode
        public static int ActualMaxNumberAccumulationsperNode = 0; // Needed Maximum Number of Clusters for a single node in Distributed Mode
        public static int MaxTransportedClusterStorage = 500;  // Maximum Number of Remote Clusters stored on a node
        public static int ActualMaxTransportedClusterStorage = 0;  // Needed Maximum Number of Remote Clusters stored on a node

        //  Specify a restart from existing file
        public static double RestartTemperature = -1.0; // If positive, Restart based on previous run(s)
        public static string RestartClusterFile = "";   // First file for restart

        //  Specify sizes of arrays dependent on cluster sizes
        public static int maxNcentperNode = 0; // maximum number of cluster centers needed in any one node
        public static int maxNcentTOTAL = 0; // maximum number of cluster centers needed in full arrays
        public static int targetNcentperPoint = 20;    // Target Maximum number of cluster centers for each point (includes Sponge)
        public static int targetMinimumNcentperPoint = 1;   // Target Minimum number of cluster centers for each point (includes Sponge)
        public static int maxNcentperPoint = 25;    // Actual Maximum number of cluster centers for each point (includes Sponge)
        public static int maxNcentCreated = 200;    // Maximum Number of clusters that can be created

        public static int InitialNcent = 1; // Initial value of Ncent

        //  Specify Annealing Schedule
        public static double InitialCoolingFactor = 0.9; // InitialCooling Factor in Annealing
        public static double FineCoolingFactor = 0.99; // Refined Cooling Factor in Annealing
        public static double InitialCoolingFactor1 = 0.9; // InitialCooling Factor in Annealing
        public static double FineCoolingFactor1 = 0.99; // Refined Cooling Factor in Annealing
        public static double InitialCoolingFactor2 = 0.9; // InitialCooling Factor in Annealing
        public static double FineCoolingFactor2 = 0.99; // Refined Cooling Factor in Annealing
        public static double CoolingTemperatureSwitch = 12.0;  // Temperature when one should switch Cooling Factors
        public static int NumberTemperatureSteps = 0;   // Count Temperature Steps

        public static double Tminimum = -1000.0;    // If Positive this is minmum Temperature, If negative Divide maximum temperature by this to get target minimum
        public static double ActualStartTemperature;    // Actual Starting Temperature;
        public static double ActualEndTemperatureafterconverging;   // Final Temperature after Converging
        public static double ActualEndTemperature;  // Temperature where we decided to stop
        public static double TimeatSplittingStop;   // Time when splitting switched off
        public static double TargetEndTemperature;  // Temperature we aimed at

        //  Parameters to specify terms that are too small
        public static double ExpArgumentCut1 = 20.0;    // For numerical stability in EM Iteration DAVectorEMIterate
        public static double ExpArgumentCut2 = 40.0; // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut2 in Point Collection (used in public void SetClustersforaPoint)
        public static double ExpArgumentCut3 = 50.0; // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut3 in distributed broadcast (only calculation 0 component) in ClusterHostlinkage in DistributedClusteringSolution
        public static double ExpArgumentCut4 = 80.0;    // For inclusion in Triangle Inequality

        //  Convergence and Iteration parameters
        public static int Waititerations = 1;   // Wait this number of Temperature iterations before splitting
        public static int Waititerations_Converge = 4;   // Wait this number of Temperature iterations before splitting after Cluster Removal

        public static int Iterationatend = 2000;  // Finish up EM Loop with at most this number of iterations
        public static int ConvergenceLoopLimit = 20; // Limit on EM Convergence for each step

        public static double FreezingLimit = 0.002; // In finish stop when all freezing measures are < FreezingLimit
        public static double Malpha_MaxChange = 0.005;  // Change in average M per point to define convergence in a step (replaces epsilon limit used in Pairwise Case)
        public static double Malpha_MaxChange1 = 0.0005;  // Change in average M per point to define convergence in a step for final steps
        public static double YChangeSquared = 0.000001;  // Change compared to width of sum of y squared changes

        public static int PowerIterationLimit = 200;   //   Limit for Power Iterations
        public static double eigenvaluechange = 0.001;   // Limit 1 on Eigenvalue Changes
        public static double eigenvectorchange = 0.001;   // Limit 2 on Eigenvalue Changes

        //  Specify Nature of Splitting
        public static int MaxNumberSplitClusters = 3;   // System will split upto this number per node (in distributed mode) simultaneously
        public static double ToosmalltoSplit = 4.0;    // Size of Cluster measured by C_k_ that should not be split
        public static double MinimumScaledWidthsquaredtosplit = 6.0;    // Do not split Cluster whose scaled width is less than this

        //  Specify annealing clean up and way it deals with close clusters
        public static double MinimumCountforCluster_C_k = 0.5; // Remove clusters with fewer points than this (based on C_k_) (This only used for converged clusters)
        public static double MinimumCountforCluster_C_kwithSponge = 1.5;
        public static int MinimumCountforCluster_Points = 2; // Remove clusters with fewer points than this (based on assigned points)
        public static double CountforCluster_C_ktobezero = 0.001; // This limit used to define a zero size cluster while system running. Positions Y are not evolved and Widths set to zero in this case

        public static double TemperatureforClosenessTest = 6.0; // Validity Check removes Close Clusters below this
        public static double ScaledSquaredDistanceatClosenessTest = 1.0;   // Coalesce Clusters separated by this
        public static double[] MagicTemperatures = { 4.0, 3.5, 2.5, 2.0, 1.75 };  // Temperatures for clean up
        public static int magicindex = 0;   // Index magic temperatures

        //  Specify Sponge Features
        public static bool UseSponge = false;   // If True use a single Sponge Factor
        public static double SpongeFactor1 = 3.0;   //  Initial Sponge Factor
        public static double SpongeFactor2 = 3.0;   // Final Sponge Factor
        public static double SpongeFactor = SpongeFactor1;   // Sponge Factor has an exp(-SpongeFactor^2/2T)
        public static int SpongePoption = 1;    // Set method to calculate sponge normalization P(SpongeCluster) = 0 default; =1 Set to average of others
        public static double SpongePWeight = 0.1;
        public static double CreateSpongeScaledSquaredWidth = -1.0; // Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches this
        public static double SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
        public static double ActualSpongeTemperature = -1.0;     // Temperature where Sponge Introduced
        public static double ActualWidth = -1.0;    // Width when Sponge introduced
        public static double TimeatSponge = 0.0;    // Time when Sponge added
        public static double SpongeTemperature2 = -1.0;     // Temperature where Final Sponge Factor introduced 
        public static double AddSpongeScaledWidthSquared = -1.0; // If Positive add sponge cluster when averaged scaled width less than this

        //  Specify Deterministic Annealing Style
        public static bool ContinuousClustering = true;    // If true use the Ken Rose Continuous Clustering
        public static bool ConvergeIntermediateClusters = false; // NOT used

        //  Specify Diagnostic Print Out
        public static int ClusterPrintNumber = 5;   // Print the last number of these clusters
        public static int PrintInterval = 3;    // Output at this interval in iterations
        public static bool RemovalDiagnosticPrint = false;  // If true output histograms on clusters per point
        public static int MaxNumberClustersToPrint = 200;   // Maximum Number of Clusters to Print
        public static bool Printeigenvectors = false;   // If true print eigenvectors in Shouldwesplit()

        //  Set Quantities that record properties of clustering
        public static double MdiffSum = 0.0;    // Accumulate Mdiff values for average
        public static int NumberMdiffSums = 0;  // Number of MdiffAvg's in Sum
        public static int IterationsperStepSum = 0; // Accumulate number of iterations per step
        public static int NumberIterationSteps = 0; // Number of steps
        public static int NumberMajorSynchs1 = 0;    // Number of calls to Major Synchromization in Global mode
        public static int NumberMajorSynchs2 = 0;    // Number of calls to Major Synchromization in Distributed mode
        public static int CountClusters2 = 0;   //  Count clusters while distributed
        public static int NumberMinorSynchs = 0;    // Number of calls to Major Synchromization in Distributed mode
        public static int NumberPipelineSteps = 0;      // Number of Pipeline steps
        public static int NumberofPipelineClusters = 0; // Number of Cluster-Steps in Pipelines
        public static int NumberPipelineGroups = 0;      // Number of Pipeline group calls
        public static string FinalReason = "";  // Final reason to stop


        public static double SumUsefulCalcs = 0.0;  // Number of Useful distance calcs
        public static double SumUselessCalcs = 0.0;  // Number of Useless distance calcs
        public static double SumIgnoredCalcs = 0.0;  // Number of Ignored distance calcs
        public static double SumEigenSPCalcs = 0.0; // Sum eigenvector scalar product computations

        public static int TotalClustersDeleted_CSmall = 0;  // Number of Clusters deleted as C_k_ too small
        public static int TotalClustersDeleted_OccCount = 0;  // Number of Clusters deleted as Occupation Count too small
        public static int TotalClustersDeleted_Close = 0;  // Number of Clusters deleted as too close

        public static double NumberMsuccesses = 0.0;  // Number of steps where M succeeds
        public static double NumberYsuccesses = 0.0;  // Number of steps where Y succeeds
        public static double NumberMfailures = 0.0; // Number of steps where M fails
        public static double NumberYfailures = 0.0; // Number of steps where Y fails
        public static double AccumulateMvalues = 0.0; // Number of iterations where M succeeeds
        public static double AccumulateYvalues = 0.0; // Number of iterations where Y succeeds

        //  Kmeans Paramters
        public static double KmeansCenterChangeStop = 0.001;    // Stop when relative to radius center change reaches this limit
        public static int KmeansIterationLimit = 1000;  // Stop at this iteration limit

        //  Triangle Inequality for Kmeans Parameters
        public static int maxNcentTOTALforParallelism_Kmeans = 50;   // Use Center parallelism when this limit hit
        public static int OldCenterOption_Kmeans = 0;  // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations
        public static bool DoBackwardFacingTests_Kmeans = true; // If True do backward facing tests
        public static int UseTriangleInequality_Kmeans = 0;    // if 0 do NOT use triangle inequality; > 0 use it in a way specified by integer
        public static double TriangleInequality_Delta1_old_Kmeans = 0.1;  // Test for center change and old lower bounds (normalized by radius)
        public static double TriangleInequality_Delta1_current_Kmeans = 0.1;  // Test for Center change and current lower bounds (normalized by radius)
        public static int MaxClusterLBsperPoint_Kmeans = 50;   // Maximum number of Lower Bound values
        public static int MaxCentersperCenter_Kmeans = 49;   // Maximum number of Centers in Center Difference Array

        //  Triangle Inequality for DA Parameters
        public static int maxNcentTOTALforParallelism_DA = 50;   // Use Center parallelism when this limit hit
        public static int OldCenterOption_DA = 0;  // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations
        public static int UseTriangleInequality_DA = 0;    // if 0 do NOT use triangle inequality; > 0 use it in a way specified by integer
        public static double TriangleInequality_Delta1_old_DA = 0.1;  // Test for center change and old lower bounds (normalized by radius)
        public static double TriangleInequality_Delta1_current_DA = 0.1;  // Test for Center change and current lower bounds (normalized by radius)
        public static int MaxClusterLBsperPoint_DA = 50;   // Maximum number of Lower Bound values
        public static int MaxCentersperCenter_DA = 49;   // Maximum number of Centers in Center Difference Array
        public static double[] ClustersperCenterDiagnostics = new double[8];    // Diagnostics

        //  Very Specific LC-MS Analysis ("print out") of Results
        public static DAVectorSponge.ArbitraryClustering GoldenPeaks;    // Certified Peaks
        public static DAVectorSponge.ArbitraryClustering MedeaClusters; // Medusa Clusters
        public static DAVectorSponge.ArbitraryClustering MclustClusters; // Mclust Clusters
        public static DAVectorSponge.ArbitraryClustering OurClusters;   // Our Clusters

        public static int CompareSolution = -1; // Control solution comparison
        public static string ComparisonClusterFile = ""; // File for Comparison Data
        public static int ComparisonSelectedInputLabel = 6;   // Charge to select in Comparison; if negative veto this label
        public static int ComparisonInputFileType = 0;    // If 0 raw data; =1 Output file

        public static int ClusterCountOutput = 0;   // Control Label Output = -1 not at all, = 0 at end only, = 1 at each count
        public static int[] ClusterAssignments;  // This gives for all points their cluster assignments (set in OutputClusterLabels)
        public static ClusterQuality ClusterStatus;
        public static int ClusterNumberOutput = -1;
        public static int NumberNearbyClusters = 5;    // specify number of nearby clusters to output
        public static double NearbySpongePointLimit = -1.0;  // Multiplier for sigma used in counting sponge points near a cluster; if negative use spongefactor

        // Config Settings
        public static ConfigurationMgr _configurationManager;
        public static ParallelOptions _parallelOptions;
        #endregion

        //  User routine to specify Sigma as a function of Cluster Position given in Centre
        public static void CalculateSigma(double[] Centre, ref double[] Sigma)
        {
            
            if (Program.SigmaMethod == 0)
            {
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    Sigma[VectorIndex] = 1.0;
                return;
            }
            if (Program.SigmaMethod == 1)
            {
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    Sigma[VectorIndex] = Program.SigmaVectorParameters_i_[VectorIndex] * Program.SigmaVectorParameters_i_[VectorIndex];
                return;
            }
            if ((Program.SigmaMethod == 2) || (Program.SigmaMethod == 3))
            {   // 
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    Sigma[VectorIndex] = Program.SigmaVectorParameters_i_[VectorIndex];
                Sigma[0] = Sigma[0] * Centre[0];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    Sigma[VectorIndex] = Sigma[VectorIndex] * Sigma[VectorIndex];
                return;
            }
            return;

        }   // End CalculateSigmaSquared

        public static bool ChangeClusterSigmas(double TemperatureRatioChange, ClusteringSolution Solution)
        {
            if (Program.SigmaMethod != 3)
                return false;

            bool asymptotic = (Program.SigmaVectorParameters_i_[0] <= Program.FinalTargetSigma0);
            DAVectorUtility.SynchronizeMPIvariable(ref asymptotic);
            if (asymptotic)
                return false;
            bool justsetit = Solution.Temperature <= Program.FinalTargetTemperature;
            DAVectorUtility.SynchronizeMPIvariable(ref justsetit);
            if (justsetit)
                Program.SigmaVectorParameters_i_[0] = Program.FinalTargetSigma0;
            else
            {
                double LogTemperatureChange = Math.Min(Math.Log(Program.FinalTargetTemperature) - Math.Log(Solution.Temperature), Math.Log(TemperatureRatioChange));
                double ActualRatio = Math.Exp(Math.Log(TemperatureRatioChange) * (Math.Log(Program.FinalTargetSigma0) - Math.Log(Program.SigmaVectorParameters_i_[0]))
                    / LogTemperatureChange);
                Program.SigmaVectorParameters_i_[0] = Math.Max(Program.FinalTargetSigma0, Program.SigmaVectorParameters_i_[0] * ActualRatio);
            }

            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                Program.CalculateSigma(Solution.Y_k_i_[RealClusterIndex], ref Solution.Sigma_k_i_[RealClusterIndex]);
            }
            return true;
        }   // End ChangeClusterSigmas

        public static bool ChangeSpongeFactor(double TemperatureRatioChange, ClusteringSolution Solution)
        {

            if (Solution.SpongeCluster < 0)
                return false;
            if ((Program.SpongeTemperature2 < 0.0) || (Program.SpongeFactor2 < 0.0))
                return false;
            bool asymptotic = (Program.SpongeFactor <= Program.SpongeFactor2);
            DAVectorUtility.SynchronizeMPIvariable(ref asymptotic);
            if (asymptotic)
                return false;

            bool justsetit = Solution.Temperature <= Program.SpongeTemperature2;
            DAVectorUtility.SynchronizeMPIvariable(ref justsetit);
            if (justsetit)
                Program.SpongeFactor = Program.SpongeFactor2;
            else
            {
                double LogTemperatureChange = Math.Min(Math.Log(Program.SpongeTemperature2) - Math.Log(Solution.Temperature), Math.Log(TemperatureRatioChange));
                double ActualRatio = Math.Exp(Math.Log(TemperatureRatioChange) * (Math.Log(Program.SpongeFactor2) - Math.Log(Program.SpongeFactor))
                    / LogTemperatureChange);
                Program.SpongeFactor = Program.SpongeFactor * ActualRatio;
                Program.SpongeFactor = Program.SpongeFactor2 + (Program.SpongeFactor1 - Program.SpongeFactor2) * (Solution.Temperature - Program.SpongeTemperature2) / (Program.SpongeTemperature1 - Program.SpongeTemperature2);
                Program.SpongeFactor = Math.Max(Program.SpongeFactor2, Program.SpongeFactor);
            }
            return true;

        }   // End ChangeSpongeFactor

        // ***** Main *******
        // Parallel Vector clustering based on 
        // Kmeans or Deterministic Annealing algorithm 
        static void Main(string[] args)
        {
            // Load the command line args into our helper class which allows us to name arguments
            Salsa.Core.Arguments pargs = new Salsa.Core.Arguments(args);
            pargs.Usage = "Usage: Salsa.Sponged2DClustering.exe /configFile=<string> /nodeCount=<int> /threadCount=<int>";

            if (pargs.CheckRequired(new string[] { "configFile", "nodeCount", "threadCount" }) == false)
            {
                Console.WriteLine(pargs.Usage);
                return;
            }

            //  Read Metadata using this as source of other metadata
            _configurationManager = ConfigurationMgr.LoadConfiguration(pargs.GetValue<string>("configFile"), true);
            ReadControlFile(_configurationManager);
            Program.ControlFileName = pargs.GetValue<string>("configFile");
            DAVectorUtility.NodeCount = pargs.GetValue<int>("nodeCount");
            DAVectorUtility.ThreadCount = pargs.GetValue<int>("threadCount");

            //  Initialize Clusters
            SigmaVectorParameters_i_ = new double[Program.ParameterVectorDimension];

            Program.DoKmeans = false;
            Program.DoLCMS = false;
            Program.DoBasicDA = false;
            Program.FollowUpKmeans = true;

            // Program.DoBasicDA = false;
             Program.DoKmeans = true;
//            Program.DoLCMS = true;

            if (Program.DoLCMS)
                Program.DoBasicDA = false;
            if (Program.DoBasicDA)
                Program.DoKmeans = false;
            if (!Program.DoBasicDA)
                Program.FollowUpKmeans = false;

            Program.SetupLCMS();
            if(!Program.DoLCMS)
                Program.CompareSolution = -1;
            Program.SetupKmeans();
            Program.SetupDA();

            if (Program.UseTriangleInequality_DA > 0)
            {
                for (int DiagnosticLoop = 0; DiagnosticLoop < 6; DiagnosticLoop++)
                    Program.ClustersperCenterDiagnostics[DiagnosticLoop] = 0.0;
            }
            if (!Program.CalculateEigenvaluesfromMatrix)
                Program.CalculateCorrelationMatrix = false;
           
            //  Don't Change anything after this
            if (!Program.Refinement)
            {
                Program.RestartTemperature = Tminimum;
                RestartInputFileType = 1;
                RestartClusterFile = _configurationManager.DAVectorSpongeSection.DistanceMatrixFile;
                _configurationManager.DAVectorSpongeSection.LabelFile = "";
            }

            // if (Program.RestartTemperature > 0.0)
            //   Program.UseSponge = true;

            Program.InitialCoolingFactor2 = Math.Max(Program.InitialCoolingFactor2, Program.InitialCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Max(Program.FineCoolingFactor2, Program.FineCoolingFactor1);
            Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
            Program.FineCoolingFactor = Program.FineCoolingFactor1;

            // Set Initial Values of Annealed Parameters
            Program.InitialSigma0 = Program.SigmaVectorParameters_i_[0];
            Program.SpongeFactor = Program.SpongeFactor1;


            //  Set up TPL
            _parallelOptions = new ParallelOptions();
            _parallelOptions.MaxDegreeOfParallelism = DAVectorUtility.ThreadCount;

            //  Set up MPI Parallelism
            DAVectorParallelism.SetupParallelism(ref args);

            //  Restrict number of splits per node
            Program.MaxNumberSplitClusters = Math.Min(Program.MaxNumberSplitClusters, DAVectorUtility.MPI_Size * 32);
            if(!Program.DoKmeans)
                Program.maxNcentperPoint = Math.Max(Program.maxNcentperPoint, Program.targetNcentperPoint + Program.MaxNumberSplitClusters + 3);
            int SaveSplitNumber = Program.MaxNumberSplitClusters;

            //  Find number of points if requested
            if (Program.NumberDataPoints < 1)
            {
                if (Program.DoLCMS)
                {
                    if (DAVectorUtility.MPI_Rank == 0)
                        DAVectorReadData.AnalyzeDataFromFile(_configurationManager.DAVectorSpongeSection.DistanceMatrixFile);
                    DAVectorUtility.MPI_communicator.Broadcast<int>(ref DAVectorUtility.PointCount_Global, 0);
                }
                else
                {
                    Exception e = DAVectorUtility.SALSAError("Must Set Number of Input Data Points");
                    throw (e);
                }
            }
            else
                DAVectorUtility.PointCount_Global = Program.NumberDataPoints;

            DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global.ToString() + " Selection " + Program.SelectedInputLabel.ToString() + " Input File Type " + Program.InputFileType.ToString() + 
                " Comparison Selection " + Program.ComparisonSelectedInputLabel.ToString() + " Comparison Input File Type " + Program.ComparisonInputFileType.ToString() +
                " Restart Selection " + Program.RestartSelectedInputLabel.ToString() + " Restart Input File Type " + Program.RestartInputFileType.ToString() );
            DAVectorUtility.SALSAPrint(1, "Files " + _configurationManager.DAVectorSpongeSection.DistanceMatrixFile + " Comparison " + Program.ComparisonClusterFile + " Restarts " + Program.RestartClusterFile + " "
                + _configurationManager.DAVectorSpongeSection.LabelFile);
            if (Program.ComparisonInputFileType == 1)
                Program.CompareSolution = -1;

            //  Read data for LC-MS Clustering Comparisons
            if (Program.CompareSolution > 0)
            {
                if(!Program.DoLCMS)
                {
                    Exception e = DAVectorUtility.SALSAError("Invalid Compare Solution Option " + Program.CompareSolution.ToString() );
                    throw (e);
                }
                
                int save1 = InputFileType;
                int save2 = SelectedInputLabel;
                InputFileType = ComparisonInputFileType;
                SelectedInputLabel = ComparisonSelectedInputLabel;
                Program.GoldenPeaks = new DAVectorSponge.ArbitraryClustering(DAVectorUtility.PointCount_Global, "Golden Peaks");
                Program.MclustClusters = new DAVectorSponge.ArbitraryClustering(DAVectorUtility.PointCount_Global, "Mclust");
                Program.MedeaClusters = new DAVectorSponge.ArbitraryClustering(DAVectorUtility.PointCount_Global, "Medea");
                Program.OurClusters = new DAVectorSponge.ArbitraryClustering(DAVectorUtility.PointCount_Global,"DAVectorSponge");
                GoldenExamination.GoldenID = new int[DAVectorUtility.PointCount_Global];
                GoldenExamination.GoldenLabel = new string[DAVectorUtility.PointCount_Global];
                GoldenExamination.PeakPosition = new double[DAVectorUtility.PointCount_Global][];
                for(int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
                    GoldenExamination.PeakPosition[GlobalPointIndex] = new double[Program.ParameterVectorDimension];
                DAVectorReadData.ReadLabelsFromFile(ComparisonClusterFile);
                InputFileType = save1;
                SelectedInputLabel = save2;
            }

            // Set up Decomposition of USED points
            DAVectorParallelism.SetParallelDecomposition();
            if (Program.Replicate > 1)
            {
                DAVectorUtility.SALSAPrint(1, "Replicate " + Program.Replicate.ToString());
                DAVectorUtility.PointCount_Global *= Program.Replicate;
                DAVectorUtility.PointCount_Largest *= Program.Replicate;
                DAVectorUtility.PointStart_Process *= Program.Replicate;
                DAVectorUtility.PointCount_Process *= Program.Replicate;
                for (int ProcessIndex = 0; ProcessIndex < DAVectorUtility.MPI_Size; ProcessIndex++)
                {
                    DAVectorUtility.PointsperProcess[ProcessIndex] *= Program.Replicate;
                    for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                        DAVectorUtility.PointsperThreadperProcess[ProcessIndex][ThreadIndex] *= Program.Replicate;
                }
                for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                {
                    DAVectorUtility.StartPointperThread[ThreadIndex] *= Program.Replicate;
                    DAVectorUtility.PointsperThread[ThreadIndex] *= Program.Replicate;
                }
            }

            // Setup PointPosition and related arrays
            Program.PointOriginalIndex = new int[DAVectorUtility.PointCount_Process]; // Point Index on File
            Program.PointLabel = new int[DAVectorUtility.PointCount_Process];
            Program.PointPosition = new double[DAVectorUtility.PointCount_Process][];
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    Program.PointPosition[alpha] = new double[Program.ParameterVectorDimension];
                }
            }); // End loop initialing Point dependent quantities

            Program.ClusterAssignments = new int[DAVectorUtility.PointCount_Global];

            Program.ClusterLimitforDistribution = Math.Max(Program.ClusterLimitforDistribution, 4 * DAVectorUtility.MPI_Size);

            // Initial Processing Complete
            DAVectorUtility.MPI_communicator.Barrier(); // Make certain all processes have processed original data before writing updated

            //  read data into memory
            if (Program.DoKmeans)
            {
                KmeansTriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_Kmeans, Program.MaxClusterLBsperPoint_Kmeans, Program.MaxCentersperCenter_Kmeans,
                    Program.TriangleInequality_Delta1_old_Kmeans, Program.TriangleInequality_Delta1_current_Kmeans, Program.OldCenterOption_Kmeans, Program.DoBackwardFacingTests_Kmeans);
                Kmeans.InitializeKmeans(Program.PointPosition, Program._parallelOptions,_configurationManager.DAVectorSpongeSection.DistanceMatrixFile,
                    Program.ClusterIndexonInputLine, Program.FirstClusterValue, Program.StartPointPositiononInputLine, Program.InitialNcent, Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_Kmeans,
                    Program.ParameterVectorDimension, Program.KmeansCenterChangeStop, Program.KmeansIterationLimit);
            }
            else if(Program.DoLCMS)
                DAVectorReadData.ReadDataFromFile(_configurationManager.DAVectorSpongeSection.DistanceMatrixFile, 0);
            else
                DAVectorReadData.ReadDataFromFile(_configurationManager.DAVectorSpongeSection.DistanceMatrixFile, Program.ClusterIndexonInputLine, null, Program.StartPointPositiononInputLine);
            // Program.ClusterIndexonInputLine MUST be NULL

            //  Read 3D Data for Plotviz
            if (_configurationManager.DAVectorSpongeSection.LabelFile.Length <= 0)
                Program.RW3DData = -1;
            if (Program.RW3DData > 0)
            {
                DAVectorUtility.SALSAPrint(0, "3D Data " + _configurationManager.DAVectorSpongeSection.LabelFile);
                Program.FullPoint3DPosition = new double[DAVectorUtility.PointCount_Global][];
                for (int alpha = 0; alpha < DAVectorUtility.PointCount_Global; alpha++)
                    Program.FullPoint3DPosition[alpha] = new double[Program.RW3DData];
                if(DAVectorUtility.MPI_Rank == 0 )
                    DAVectorReadData.Read3DDataFromFile(_configurationManager.DAVectorSpongeSection.LabelFile, Program.RW3DData, 1);
            }

            //  Set up Timing
            int nonMPITimings = 14;
            DAVectorUtility.InitializeTiming(13 + nonMPITimings);
            DAVectorUtility.SetUpMPISubTimers(nonMPITimings, "");
            DAVectorUtility.SetUpSubTimer(0, "Decide Splitting");
            DAVectorUtility.SetUpSubTimer(1, "Validity");
            DAVectorUtility.SetUpSubTimer(2, "EMIterate-1");
            DAVectorUtility.SetUpSubTimer(3, "EMIterate-2");
            DAVectorUtility.SetUpSubTimer(4, "ClusterAvgs");
            DAVectorUtility.SetUpSubTimer(5, "Do Splitting");
            DAVectorUtility.SetUpSubTimer(6, "Major Synch - 1");
            DAVectorUtility.SetUpSubTimer(7, "Minor Synch");
            DAVectorUtility.SetUpSubTimer(8, "Final Status");
            DAVectorUtility.SetUpSubTimer(9, "EMIterate-3");
            DAVectorUtility.SetUpSubTimer(10, "EMIterate-4");
            DAVectorUtility.SetUpSubTimer(11, "Major Synch - 2");
            DAVectorUtility.SetUpSubTimer(12, "Major Synch - 3");
            DAVectorUtility.SetUpSubTimer(13, "Redo Clusters per Point");
            DAVectorUtility.TimingOutputOrder[0] = 0; 
            DAVectorUtility.TimingOutputOrder[1] = 5;
            DAVectorUtility.TimingOutputOrder[2] = 1;
            DAVectorUtility.TimingOutputOrder[3] = 13;
            DAVectorUtility.TimingOutputOrder[4] = 4;
            DAVectorUtility.TimingOutputOrder[5] = 8;
            DAVectorUtility.TimingOutputOrder[6] = 2;
            DAVectorUtility.TimingOutputOrder[7] = 3;
            DAVectorUtility.TimingOutputOrder[8] = 9;
            DAVectorUtility.TimingOutputOrder[9] = 10;
            DAVectorUtility.TimingOutputOrder[10] = 7;
            DAVectorUtility.TimingOutputOrder[11] = 6;
            DAVectorUtility.TimingOutputOrder[12] = 11;
            DAVectorUtility.TimingOutputOrder[13] = 12;
            if ((!Program.DoLCMS) && (!Program.DoKmeans))
            {
                DAVectorUtility.SetUpSubTimer(6, "Full Triangle Ineq");
                DAVectorUtility.SetUpSubTimer(11, "CenterFacing");
                DAVectorUtility.SetUpSubTimer(12, "Point LB");
            }
            if (Program.DoKmeans)
            {
                DAVectorUtility.SetUpSubTimer(0, "Pure Kmeans");
                DAVectorUtility.SetUpSubTimer(1, "End Iteration");
                DAVectorUtility.SetUpSubTimer(2, "Center Facing");
                DAVectorUtility.SetUpSubTimer(3, "Point LB");
                DAVectorUtility.SetUpSubTimer(4, "Cluster Centers");
                DAVectorUtility.SetUpSubTimer(5, "Individual Centers");
                DAVectorUtility.SetUpSubTimer(6, "Center Centers");
                DAVectorUtility.SetUpSubTimer(7, "Update Centers");
                for (int timingloop = 9; timingloop < nonMPITimings; timingloop++)
                    DAVectorUtility.SubTimingEnable[timingloop] = false;
                for(int timingloop =0; timingloop < nonMPITimings; timingloop++)
                    DAVectorUtility.TimingOutputOrder[timingloop] = timingloop;
            }

            //  Set up basic clusters
            ClusteringSolution.SetParameters(DAVectorUtility.PointCount_Process, Program.maxNcentCreated, Program.maxNcentTOTAL, Program.maxNcentperNode,
                Program.cachelinesize, Program.targetNcentperPoint, Program.targetMinimumNcentperPoint, Program.maxNcentperPoint, Program.ExpArgumentCut2);

            //  Set up Distributed Clusters
            DistributedClusteringSolution DistributedSetup;
            if(Program.DoLCMS)
                DistributedSetup = new DistributedClusteringSolution(MaxMPITransportBuffer, MaxTransportedClusterStorage, MaxNumberAccumulationsperNode,
                MaxDoubleComponents, MaxIntegerComponents);
            DAVectorUtility.SALSAPrint(0, "Setup Finished");

            //  Do Clustering
            VectorAnnealIterate RunVectorSpongeDA;
            ControlKmeans RunKmeansClustering;
            if (DoKmeans)
            {
                RunKmeansClustering = new ControlKmeans();
            }
            else
            {
                RunVectorSpongeDA = new VectorAnnealIterate();
                RunVectorSpongeDA.ControlVectorSpongeDA();
            }
            Program.ActualEndTemperatureafterconverging = ParallelClustering.RunningSolution.Temperature;

            //  End Timing
            DAVectorUtility.EndTiming();
            
            // Calculate Cluster Statistics
            DAVectorUtility.StartSubTimer(8);

            // Calculate Occupation Counts and other statistics for LCMS
            //  This can alter "Sponge Confused Points"
            if(!Program.DoKmeans)
                ParallelClustering.RunningSolution.FindOccupationCounts();

            if (Program.ClusterCountOutput >= 0 && (!Program.DoKmeans) )
                VectorAnnealIterate.OutputClusteringResults("Final");

            if(Program.DoBasicDA && (Program.RW3DData > 0))
                VectorAnnealIterate.Output3DClusterLabels("DA");
            if (Program.DoKmeans)
            {
                double[][] ClusterCenters;
                
                ControlKmeans.CaptureKmeans(Program.ClusterAssignments, out ClusterCenters);
                if (Program.ClusterCountOutput >= 0)
                {
                    VectorAnnealIterate.SimpleOutputClusteringResults("Kmeans", ClusterCenters);
                    DAVectorUtility.SALSAPrint(0, "End Output of Full Clusters");
                }
                if (Program.RW3DData > 0)
                {
                    VectorAnnealIterate.Output3DClusterLabels("Kmeans");
                    DAVectorUtility.SALSAPrint(0, "End Output of 3D Clusters");
                }
            }
            if (Program.DoLCMS)
                VectorAnnealIterate.CalculateClusterStatus();

            // Output Results
            string nextline = "\nNode 0 Center Averages (Counts) [Distce] ";
            int Totnumber = ClusteringSolution.NumberLocalActiveClusters;
            if (Totnumber > Program.MaxNumberClustersToPrint)
            {
                nextline += Totnumber.ToString() + " Truncated at " + Program.MaxNumberClustersToPrint.ToString() + " ";
                Totnumber = Program.MaxNumberClustersToPrint;
            }
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < Totnumber; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                double meandist = ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex];
                string cees = "";
                if (!Program.DoKmeans)
                    cees = ParallelClustering.RunningSolution.C_k_[RealClusterIndex].ToString("F2") + " - ";
                nextline += RealClusterIndex.ToString() + " " + cees 
                     + ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex].ToString()
                    + " [Wid " + meandist.ToString("F1") + "]";
                if (!Program.DoKmeans)
                    nextline += "[Frz " + ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex].ToString("E2") + "]";
                nextline += " * ";
            }
            bool save = DAVectorUtility.ConsoleDebugOutput;
            if( Totnumber > 100)
                DAVectorUtility.ConsoleDebugOutput = false;
            DAVectorUtility.SALSAPrint(0, nextline);
            DAVectorUtility.ConsoleDebugOutput = save;

            if (!Program.DoKmeans)
            {
                double tmp1 = 0.0;
                if (Program.NumberMdiffSums > 0)
                    Program.MdiffSum = Program.MdiffSum / Program.NumberMdiffSums;
                if( Program.NumberIterationSteps > 0 )
                    tmp1 = (double)Program.IterationsperStepSum / (double)Program.NumberIterationSteps;

                DAVectorUtility.SALSAPrint(0, "\nT " + ParallelClustering.RunningSolution.Temperature.ToString("F5") + " Cluster " + ParallelClustering.RunningSolution.Ncent_Global.ToString()
                    + " Iter " + VectorAnnealIterate.EMIterationCount.ToString() + " Extra Iter " + VectorAnnealIterate.Extra_EMIterationCount.ToString()
                    + " Average Mdiff " + Program.MdiffSum.ToString("F5") + " Number of Steps "
                    + Program.NumberIterationSteps.ToString() + " Iterations per Step " + tmp1.ToString("F2"));
            }

            DAVectorUtility.SALSAPrint(0, "\n" + DAVectorUtility.PatternLabel);
            DAVectorUtility.SALSAPrint(0, "Labels File: " + _configurationManager.DAVectorSpongeSection.ClusterFile);
            DAVectorUtility.SALSAPrint(0, "Timing Output: " + _configurationManager.DAVectorSpongeSection.TimingFile);
            string message = " Charge Value ";
            if (InputFileType == 1)
                message = " Selected with Cluster ";
            if (Program.DoKmeans)
                message = " Selection ";
            DAVectorUtility.SALSAPrint(0, "Data Points " + DAVectorUtility.PointCount_Global.ToString() + message + Program.SelectedInputLabel.ToString());
            DAVectorUtility.SALSAPrint(0, "Vector Dimension: " + Program.ParameterVectorDimension.ToString());
            DAVectorUtility.SALSAPrint(0, "Continuous Clustering: " + Program.ContinuousClustering.ToString());
           
            string sigmamethodString = "Sigma Method: " + Program.SigmaMethod.ToString();
            if (Program.SigmaMethod > 0)
            {
                sigmamethodString += " Parameters:";
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    sigmamethodString += " " + Program.SigmaVectorParameters_i_[VectorIndex].ToString("E4");
            }
            DAVectorUtility.SALSAPrint(0, sigmamethodString);
            if (Program.SigmaMethod == 3)
                DAVectorUtility.SALSAPrint(0, "Initial Sigma[0] " + Program.InitialSigma0.ToString("E4") + " Final Sigma[0] = " + Program.FinalTargetSigma0.ToString("E4") +
                    " /Center[0] at Temperature " + Program.FinalTargetTemperature.ToString("F4"));

            DAVectorUtility.SALSAPrint(0, "Initial Number of Centers: " + Program.InitialNcent.ToString());
            DAVectorUtility.SALSAPrint(0, "Final Number of Centers: " + ParallelClustering.RunningSolution.Ncent_Global.ToString());
            DAVectorUtility.SALSAPrint(0, "Maximum Number of Centers: " + Program.maxNcentperNode.ToString());
            DAVectorUtility.SALSAPrint(0, "Target Maximum number of cluster centers for each point (includes Sponge): " + Program.targetNcentperPoint.ToString());
            DAVectorUtility.SALSAPrint(0, "Actual Maximum number of cluster centers for each point (includes Sponge): " + Program.maxNcentperPoint.ToString());

            if (!Program.DoKmeans)
            {
                ActualMaxMPITransportBuffer = DAVectorUtility.MPI_communicator.Allreduce<int>(ActualMaxMPITransportBuffer, Operation<int>.Max);
                ActualMaxNumberAccumulationsperNode = DAVectorUtility.MPI_communicator.Allreduce<int>(ActualMaxNumberAccumulationsperNode, Operation<int>.Max);
                ActualMaxTransportedClusterStorage = DAVectorUtility.MPI_communicator.Allreduce<int>(ActualMaxTransportedClusterStorage, Operation<int>.Max);
                DAVectorUtility.SALSAPrint(0, "Needed Maximum: MaxMPITransportBuffer " + ActualMaxMPITransportBuffer.ToString() + " MaxNumberAccumulationsperNode " + ActualMaxNumberAccumulationsperNode.ToString()
                    + " MaxTransportedClusterStorage "  + ActualMaxTransportedClusterStorage.ToString());
                int MaxCreatedIndex = ClusteringSolution.CenterMaxforCreatedIndex;
                MaxCreatedIndex = DAVectorUtility.MPI_communicator.Allreduce<int>(MaxCreatedIndex, Operation<int>.Max);
                DAVectorUtility.SALSAPrint(0, "Created Index Space " + MaxCreatedIndex.ToString() + " times " + ClusteringSolution.PACKINGMULTIPLIER);

                DAVectorUtility.SALSAPrint(0, "Iterations at the End: " + Program.Iterationatend.ToString());
                DAVectorUtility.SALSAPrint(0, "Limit on EM Convergence for each step: " + Program.ConvergenceLoopLimit.ToString());
                DAVectorUtility.SALSAPrint(0, "Change in Average M per point: " + Program.Malpha_MaxChange.ToString("E4") + " In final Loop " + Program.Malpha_MaxChange1.ToString("E4")
                    + " Y Change (final only) " + YChangeSquared.ToString("E4") );
                DAVectorUtility.SALSAPrint(0, "Freezing Limit for convergence: " + Program.FreezingLimit.ToString("E4"));
                DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Iteration Code " + Program.ExpArgumentCut1.ToString("F4"));
                DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Associating Clusters with Points " + Program.ExpArgumentCut2.ToString("F4") + " (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut");
                DAVectorUtility.SALSAPrint(0, "Mean Cluster Count per point " + VectorAnnealIterate.MeanClusterCount.ToString("F2") + " Pts with Just 1 Cluster " + VectorAnnealIterate.PointswithClusterCount1.ToString("F0"));
                
                DAVectorUtility.SALSAPrint(0, "Tmininimum " + Program.Tminimum.ToString("F4") + " If Positive this is minmum Temperature, If negative Divide maximum temperature by this to get target minimum");

                DAVectorUtility.SALSAPrint(0, "\nNumber of Clusters Split Simultaneously: " + Program.MaxNumberSplitClusters.ToString() + " Initially " + SaveSplitNumber.ToString());

                DAVectorUtility.SALSAPrint(0, "Do not split Clusters smaller than this: " + Program.ToosmalltoSplit.ToString("F3"));
                DAVectorUtility.SALSAPrint(0, "Do not split Clusters with Scaled Squared Width Less than this " + Program.MinimumScaledWidthsquaredtosplit.ToString("F3"));
                DAVectorUtility.SALSAPrint(0, "Wait stages between splits: " + Program.Waititerations.ToString());
                DAVectorUtility.SALSAPrint(0, "Initial Cooling Factor in Annealing: " + Program.InitialCoolingFactor1.ToString("F7") + " " + Program.InitialCoolingFactor2.ToString("F7"));
                DAVectorUtility.SALSAPrint(0, "Refined Cooling Factor in Annealing: " + Program.FineCoolingFactor1.ToString("F7") + " " + Program.FineCoolingFactor2.ToString("F7")
                    + " Switching at " + Program.CoolingTemperatureSwitch.ToString("F3"));

                int cleanuplength = Program.MagicTemperatures.Length;
                message = "Clean up Temperatures ";
                for (int loop = 0; loop < cleanuplength; loop++)
                    message += Program.MagicTemperatures[loop].ToString("F2") + " ";
                DAVectorUtility.SALSAPrint(0, message);

                if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                {
                    DAVectorUtility.SALSAPrint(0, "\nSponge Cluster: " + ParallelClustering.RunningSolution.SpongeCluster.ToString());
                    if (Program.UseSponge)
                        DAVectorUtility.SALSAPrint(0, "Sponge Created Initially");
                    DAVectorUtility.SALSAPrint(0, "Final Sponge Factor: " + Program.SpongeFactor.ToString("F2") + " Initial " + Program.SpongeFactor1.ToString("F2") + " Started at " + Program.SpongeTemperature1.ToString("F3")
                        + " Aiming at " + Program.SpongeFactor2.ToString("F2") + " at Temperature " + Program.SpongeTemperature2.ToString("F3"));
                    DAVectorUtility.SALSAPrint(0, "Sponge p(k) Option: " + Program.SpongePoption.ToString() + " Weighting " + Program.SpongePWeight.ToString("E3"));
                    DAVectorUtility.SALSAPrint(0, "Nearby Sponge Point Limit: " + Program.NearbySpongePointLimit.ToString());
                    if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                        DAVectorUtility.SALSAPrint(0, " Points in Sponge Cluster " + ParallelClustering.RunningSolution.C_k_[ParallelClustering.RunningSolution.SpongeCluster].ToString("F2") + " Occupation Count " +
                            ParallelClustering.RunningSolution.OccupationCounts_k_[ParallelClustering.RunningSolution.SpongeCluster].ToString());
                    DAVectorUtility.SALSAPrint(0, "Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches " + Program.CreateSpongeScaledSquaredWidth.ToString("F2"));
                    if (Program.ActualSpongeTemperature > 0.0)
                    {
                        double FractionTimeatSponge = Program.TimeatSponge / DAVectorUtility.HPDuration;
                        DAVectorUtility.SALSAPrint(0, "Temperature where Sponge Introduced " + Program.ActualSpongeTemperature.ToString("F4")
                            + " Width Then " + Program.ActualWidth.ToString("F3") + " Time fraction " + FractionTimeatSponge.ToString("F4"));
                    }
                }
            }

            DAVectorUtility.SALSAPrint(0, "\nCluster Deletion Statistics: C_k_ Small " + Program.TotalClustersDeleted_CSmall.ToString() + " Occ. Count Small " + Program.TotalClustersDeleted_OccCount.ToString()
                + " Clusters Close " + Program.TotalClustersDeleted_Close.ToString());
            DAVectorUtility.SALSAPrint(0, "Minimum Cluster Count Average Point for Absorbing: " + Program.MinimumCountforCluster_C_k.ToString("F2")
                                            + " With Sponge " + Program.MinimumCountforCluster_C_kwithSponge.ToString("F2"));
            DAVectorUtility.SALSAPrint(0, "Minimum Cluster Identified Point Count for Absorbing: " + Program.MinimumCountforCluster_Points.ToString());
            DAVectorUtility.SALSAPrint(0, "Minimum Cluster Average Count for considering as nonzero " + Program.CountforCluster_C_ktobezero.ToString("F5"));
            DAVectorUtility.SALSAPrint(0, "Scaled Squared Distance Cut for Closeness " + Program.ScaledSquaredDistanceatClosenessTest.ToString("F2") + " Temperature " + Program.TemperatureforClosenessTest.ToString("F2"));

            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                double FractionTimeatDist = Program.TimeatDistribution / DAVectorUtility.HPDuration;
                DAVectorUtility.SALSAPrint(0, "\nDistributed Execution entered at Temperature " + Program.ActualTemperatureforDistribution.ToString("F4") + " Time(fraction) " + FractionTimeatDist.ToString("F4") +
                    " Cut " + Program.TemperatureLimitforDistribution.ToString("F4") + " Clusters "
                    + Program.ActualClusterNumberforDistribution.ToString() + " Cut " + Program.ClusterLimitforDistribution.ToString());
                DAVectorUtility.SALSAPrint(0, "Exponential Cut Used in Distributed Broadcast " + Program.ExpArgumentCut3.ToString("F4") + " (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut -- 0th component only");
                DAVectorUtility.SALSAPrint(0, "Number of Minor Synchronizations " + Program.NumberMinorSynchs.ToString() + " Number of Major Synchs while Global " + Program.NumberMajorSynchs1.ToString()
                    + " Number of Major Synchs while Distributed " + Program.NumberMajorSynchs2.ToString());
                double ClustersperPipeline = 0.0;
                if (Program.NumberPipelineSteps > 0)
                    ClustersperPipeline = (double) Program.NumberofPipelineClusters / (double) Program.NumberPipelineSteps;
                DAVectorUtility.SALSAPrint(0, " Number of Pipeline steps " + Program.NumberPipelineSteps.ToString() + " with average cluster load " + ClustersperPipeline.ToString("F3") + " In " + Program.NumberPipelineGroups.ToString() + " sets");
                double AverageClusterCount = (double)Program.CountClusters2 / (double)Program.NumberMajorSynchs2;
                double[] AverageClusterCountvalues = new double[DAVectorUtility.MPI_Size];
                AverageClusterCountvalues = DAVectorUtility.MPI_communicator.Allgather<double>(AverageClusterCount);
                double maxClusterCount = 0.0;
                double minClusterCount = 1.0E+08;
                double avgClusterCount = 0.0;
                for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
                {
                    avgClusterCount += AverageClusterCountvalues[nodes];
                    if( AverageClusterCountvalues[nodes] > maxClusterCount)
                        maxClusterCount = AverageClusterCountvalues[nodes];
                    if (AverageClusterCountvalues[nodes] < minClusterCount)
                        minClusterCount = AverageClusterCountvalues[nodes];
                }
                avgClusterCount = avgClusterCount / DAVectorUtility.MPI_Size;
                AverageClusterCount = (double) ParallelClustering.RunningSolution.Ncent_ThisNode;
                AverageClusterCountvalues = DAVectorUtility.MPI_communicator.Allgather<double>(AverageClusterCount);
                double maxClusterCountFinal = 0.0;
                double minClusterCountFinal = 1.0E+08;
                double avgClusterCountFinal = 0.0;
                for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
                {
                    avgClusterCountFinal += AverageClusterCountvalues[nodes];
                    if (AverageClusterCountvalues[nodes] > maxClusterCountFinal)
                        maxClusterCountFinal = AverageClusterCountvalues[nodes];
                    if (AverageClusterCountvalues[nodes] < minClusterCountFinal)
                        minClusterCountFinal = AverageClusterCountvalues[nodes];
                }
                avgClusterCountFinal = avgClusterCountFinal / DAVectorUtility.MPI_Size;
                DAVectorUtility.SALSAPrint(0, "Clusters per node Average " + avgClusterCount.ToString("F2") + " Min " + minClusterCount.ToString("F2") + " Max " + maxClusterCount.ToString("F2"));
                DAVectorUtility.SALSAPrint(0, "Clusters per node Final " + avgClusterCountFinal.ToString("F2") + " Min " + minClusterCountFinal.ToString("F2") + " Max " + maxClusterCountFinal.ToString("F2"));
            }

            DAVectorUtility.SALSAPrint(0, "\nActual Starting Temperature " + Program.ActualStartTemperature.ToString("F4"));
            DAVectorUtility.SALSAPrint(0, "Final Temperature after Converging " + Program.ActualEndTemperatureafterconverging.ToString("F4"));
            double FractionTimeatSplitStop = Program.TimeatSplittingStop / DAVectorUtility.HPDuration;
            DAVectorUtility.SALSAPrint(0, "Temperature where we decided to stop " + Program.ActualEndTemperature.ToString("F4") + " Fractional Time " + FractionTimeatSplitStop.ToString("F4"));
            DAVectorUtility.SALSAPrint(0, "Temperature where we aimed to stop at " + Program.TargetEndTemperature.ToString("F4"));
            if (!Program.DoKmeans)
            {
                double TotalCalcs = Program.SumUsefulCalcs + Program.SumUselessCalcs + Program.SumIgnoredCalcs;
                DAVectorUtility.SALSAPrint(0, "\nNumber of distance calcs in EMIterate " + TotalCalcs.ToString("E4") + " Useful " + (Program.SumUsefulCalcs / TotalCalcs).ToString("F4")
                    + " Useless " + (Program.SumUselessCalcs / TotalCalcs).ToString("F4") + " Ignored " + (Program.SumIgnoredCalcs / TotalCalcs).ToString("F4") );
                DAVectorUtility.SALSAPrint(0, "Scalar Products in Eigenvector part " + Program.SumEigenSPCalcs.ToString("E4") + " Sum Useful + Useless " + (Program.SumUsefulCalcs + Program.SumUselessCalcs).ToString("E4") );

                double MIterationAverage = Program.AccumulateMvalues / Program.NumberMsuccesses;
                DAVectorUtility.SALSAPrint(0, "Mchange Test Failures " + Program.NumberMfailures.ToString("F0") + " Successes " + Program.NumberMsuccesses.ToString("F0") +
                    " Average Iterations " + MIterationAverage.ToString("F2") + " Max " + Program.ConvergenceLoopLimit.ToString());
                double YIterationAverage = Program.AccumulateYvalues / Program.NumberYsuccesses;
                DAVectorUtility.SALSAPrint(0, "Ychange Test Failures " + Program.NumberYfailures.ToString("F0") + " Successes " + Program.NumberYsuccesses.ToString("F0") +
                    " Average Iterations " + YIterationAverage.ToString("F2") + " Max " + Program.ConvergenceLoopLimit.ToString());
            }
            if (Program.UseTriangleInequality_DA > 0)
            {
                double weighting = 1.0 / Program.ClustersperCenterDiagnostics[0];
                for (int DiagnosticLoop = 1; DiagnosticLoop < 8; DiagnosticLoop++)
                    Program.ClustersperCenterDiagnostics[DiagnosticLoop] *= weighting;
                DAVectorUtility.SALSAPrint(0, "Dynamic Clusters in Triangle Inequality Failures " + Program.ClustersperCenterDiagnostics[1].ToString("F4") + " Succeses " + Program.ClustersperCenterDiagnostics[2].ToString("F4")
                    + " Very Small " + Program.ClustersperCenterDiagnostics[3].ToString("F4") + " Small " + Program.ClustersperCenterDiagnostics[4].ToString("F4")
                    + " Too Big " + Program.ClustersperCenterDiagnostics[5].ToString("F4") + " Average reset M " + Program.ClustersperCenterDiagnostics[6].ToString("F4") +
                    " Average Changes in Clusters per Point " + Program.ClustersperCenterDiagnostics[7].ToString("F4") );
            }

            DAVectorUtility.SALSAPrint(0, "\nPoints " + DAVectorUtility.PointCount_Global.ToString() + " Selection " + Program.SelectedInputLabel.ToString() + " Number of Clusters "
                + ParallelClustering.RunningSolution.Ncent_Global.ToString() + " Temperature Steps " + Program.NumberTemperatureSteps.ToString() + " " + Program.FinalReason);

            //  Cluster Widths
            string widthmessage = "";
            if (!Program.CalculateIndividualWidths)
                widthmessage = " Average Width " + ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth.ToString("E3");
            else
            {
                widthmessage = " Average Widths excluding Sponge Points (These are in Hamiltonian value) ";
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    widthmessage += ParallelClustering.RunningSolution.AverageWidth[VectorIndex].ToString("E3") + " ";
            } 
            DAVectorUtility.SALSAPrint(0, widthmessage + " These are squared and divided by squares of appropriate sigmas with Hamiltonian " + ParallelClustering.RunningSolution.PairwiseHammy.ToString("E4"));

            if (Program.DoBasicDA)
            {   //  General Statistics
                ClusterQuality.CaculateTemperatureClusterCountPlot();
                ClusteringSolution.TotalClusterSummary.HistogramClusterProperties();
                if(Program.UseTriangleInequality_DA > 0 )
                    DATriangleInequality.PrintDiagnostics();
            }

            /* Compute the duration between the initial and the end time ignoring print out. */
            DAVectorUtility.StopSubTimer(8);
            TimeSpan duration = DAVectorUtility.endTime - DAVectorUtility.startTime;
            DAVectorUtility.SALSAPrint(0, "\nTotal Time excluding I/O  " + duration.ToString() + " " + (DAVectorUtility.HPDuration * .001).ToString("F0"));
            nextline = "Partial Times ";
            for (int jtimer = 0; jtimer < DAVectorUtility.NumberofSubTimings; jtimer++)
            {
                int itimer = DAVectorUtility.TimingOutputOrder[jtimer];
                if (itimer == nonMPITimings)
                    nextline += "\n\n";
                if (DAVectorUtility.SubTimingEnable[itimer])
                {
                    double tmp = DAVectorUtility.SubDurations[itimer] / DAVectorUtility.HPDuration;
                    nextline += DAVectorUtility.SubTimingNames[itimer] + " Time=" + (DAVectorUtility.SubDurations[itimer] * .001).ToString("F0") + " #=" + DAVectorUtility.SubTimingCalls[itimer].ToString() + " " + tmp.ToString("F4") + " ** ";
                }
            }
            DAVectorUtility.SALSAPrint(0, nextline);

            //  Add a Fast K means step using triangle inequality speed up
            if (Program.DoBasicDA && Program.FollowUpKmeans && (ParallelClustering.RunningSolution.SpongeCluster < 0)
                && !ParallelClustering.RunningSolution.DistributedExecutionMode )
            {
                Program.DoKmeans = true;

                Program.TemperatureLimitforDistribution = -1.0;
                Program.ClusterLimitforDistribution = -1;

                Program.CalculateEigenvaluesfromMatrix = false;
                Program.UseSponge = false;
                Program.ContinuousClustering = false;
                Program.SigmaMethod = 0;
                Program.targetNcentperPoint = ParallelClustering.RunningSolution.Ncent_Global;
                Program.maxNcentperPoint = ParallelClustering.RunningSolution.Ncent_Global;
                Program.targetMinimumNcentperPoint = 2;

                Program.InitialNcent = ParallelClustering.RunningSolution.Ncent_Global;
                Program.maxNcentperNode = ParallelClustering.RunningSolution.Ncent_Global;
                Program.maxNcentTOTAL = ParallelClustering.RunningSolution.Ncent_Global;
                Program.maxNcentTOTALforParallelism_Kmeans = 50;
                Program.maxNcentCreated = 200;
                Program.NumberDataPoints = 200000;

                Program.UseTriangleInequality_Kmeans = 0;  // 0 is Pure K means
                Program.MaxClusterLBsperPoint_Kmeans = ParallelClustering.RunningSolution.Ncent_Global;
                Program.MaxCentersperCenter_Kmeans = ParallelClustering.RunningSolution.Ncent_Global;

                Program.TriangleInequality_Delta1_old_Kmeans = 0.2;
                Program.TriangleInequality_Delta1_current_Kmeans = 0.2;
                Program.KmeansCenterChangeStop = 0.00001;
                Program.KmeansIterationLimit = 1000;

                KmeansTriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_Kmeans, Program.MaxClusterLBsperPoint_Kmeans, Program.MaxCentersperCenter_Kmeans,
                Program.TriangleInequality_Delta1_old_Kmeans, Program.TriangleInequality_Delta1_current_Kmeans, Program.OldCenterOption_Kmeans, Program.DoBackwardFacingTests_Kmeans);
                Kmeans.InitializeKmeans(Program.PointPosition, Program._parallelOptions, "",
                    Program.ClusterIndexonInputLine, Program.FirstClusterValue, Program.StartPointPositiononInputLine, Program.InitialNcent, Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_Kmeans,
                    Program.ParameterVectorDimension, Program.KmeansCenterChangeStop, Program.KmeansIterationLimit);
                Kmeans.SetupKmeans(Program.ClusterAssignments);
                RunKmeansClustering = new ControlKmeans();

                double[][] ClusterCenters;
                ControlKmeans.CaptureKmeans(Program.ClusterAssignments, out ClusterCenters);
                if (Program.ClusterCountOutput >= 0)
                    VectorAnnealIterate.SimpleOutputClusteringResults("Kmeans", ClusterCenters);
                if(Program.RW3DData > 0)
                    VectorAnnealIterate.Output3DClusterLabels("Kmeans");

            }

            if (DAVectorUtility.MPI_Rank == 0)
            {
                DAVectorUtility.WriteResults_Cluster(_configurationManager.DAVectorSpongeSection.SummaryFile, DAVectorUtility.CosmicOutput);

                WriteTimingFile(_configurationManager.DAVectorSpongeSection.TimingFile,
                    duration,
                    DAVectorUtility.HPDuration,
                    DAVectorUtility.ThreadCount,
                    DAVectorUtility.MPIperNodeCount,
                    DAVectorUtility.NodeCount,
                    DAVectorUtility.PointCount_Process,
                    DAVectorUtility.PointCount_Global,
                    Program.maxNcentperNode,
                    DAVectorUtility.SubDurations[0] * 0.001,  // convert to milliseconds
                    DAVectorUtility.SubDurations[0] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[1] * 0.001,  // convert to milliseconds
                    DAVectorUtility.SubDurations[1] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[2] * 0.001,  // convert to milliseconds
                    DAVectorUtility.SubDurations[2] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[3] * 0.001,  // convert to milliseconds
                    DAVectorUtility.SubDurations[3] / DAVectorUtility.HPDuration,
                    DAVectorUtility.SubDurations[4] * 0.001,  // convert to milliseconds
                    DAVectorUtility.SubDurations[4] / DAVectorUtility.HPDuration,
                    _configurationManager.DAVectorSpongeSection.DistanceMatrixFile,
                    DateTime.Now,
                    MPI.Environment.ProcessorName);
            }

            DAVectorParallelism.TearDownParallelism(); //  Finalize MPI
            return;

        }
        // End Main



        public static void ReadControlFile(ConfigurationMgr mgr)
        {
            DAVectorSpongeSection section = mgr.DAVectorSpongeSection;

            Program.NumberDataPoints = section.NumberDataPoints;
            Program.ProcessingOption = section.ProcessingOption;
            DAVectorUtility.ParallelPattern = section.Pattern;
            Program.maxNcentperNode = section.MaxNcentPerNode;
            DAVectorUtility.ThreadCount = section.ThreadCount;     // Number of Threads
            DAVectorUtility.NodeCount = section.NodeCount;       // Number of Nodes
            DAVectorUtility.MPIperNodeCount = section.MPIPerNodeCount; // Number of MPI processes per Node
            Program.ToosmalltoSplit = section.TooSmallToSplit;
            Program.Waititerations = section.WaitIterations;
            Program.InitialCoolingFactor = section.InitialCoolingFactor;  // InitialCooling Factor in Annealing
            Program.FineCoolingFactor = section.FineCoolingFactor;    // Refined Cooling Factor in Annealing
            Program.Malpha_MaxChange = section.MalphaMaxChange;     //converge test condition for change in epsi
            Program.Iterationatend = section.IterationAtEnd;       // Finish up EM Loop with this number of iterations
            Program.ConvergenceLoopLimit = section.ConvergenceLoopLimit; // Limit on EM Convergence
            Program.FreezingLimit = section.FreezingLimit;       // In finish stop when all freezing measures are < FreezingLimit
            Program.ConvergeIntermediateClusters = section.ConvergeIntermediateClusters;
            DAVectorUtility.DebugPrintOption = section.DebugPrintOption;
            DAVectorUtility.ConsoleDebugOutput = section.ConsoleDebugOutput;

            Program.MaxMPITransportBuffer = section.MaxMPITransportBuffer;
            Program.MaxNumberAccumulationsperNode = section.MaxNumberAccumulationsPerNode;
            Program.MaxTransportedClusterStorage = section.MaxTransportedClusterStorage;
            Program.SpongeFactor1 = section.SpongeFactor1;
            Program.SpongeFactor2 = section.SpongeFactor2;
            Program.SpongeTemperature1 = section.SpongeTemperature1;
            Program.SpongeTemperature2 = section.SpongeTemperature2;
            Program.SpongePoption = section.SpongePOption;
            Program.SpongePWeight = section.SpongePWeight;
            Program.UseSponge = section.UseSponge;
            Program.MaxNumberSplitClusters = section.MaxNumberSplitClusters;
            Program.maxNcentTOTAL = section.MaxNcentTotal;
            Program.maxNcentperNode = section.MaxNcentPerNode;
            Program.maxNcentperPoint = section.MaxNcentPerPoint;
            Program.MaxDoubleComponents = section.MaxDoubleComponents;
            Program.MaxIntegerComponents = section.MaxIntegerComponents;
            Program.MinimumScaledWidthsquaredtosplit = section.ScaledWidthSquaredToSplit;
            Program.RestartTemperature = section.RestartTemperature;
        }

        public static void WriteControlFile(ConfigurationMgr mgr)
        {
            DAVectorSpongeSection section = mgr.DAVectorSpongeSection;

            section.SpongeTemperature1 = Program.SpongeTemperature1;
            section.SpongeTemperature2 = Program.SpongeTemperature2;
            section.RestartTemperature = Program.RestartTemperature;

            section.NumberDataPoints = Program.NumberDataPoints;
            section.ProcessingOption = Program.ProcessingOption;
            section.Pattern = DAVectorUtility.ParallelPattern;
            section.MaxNcentPerNode = Program.maxNcentperNode;
            section.ThreadCount = DAVectorUtility.ThreadCount;
            section.NodeCount = DAVectorUtility.NodeCount;
            section.MPIPerNodeCount = DAVectorUtility.MPIperNodeCount;
            section.TooSmallToSplit = Program.ToosmalltoSplit;
            section.WaitIterations = Program.Waititerations;
            section.InitialCoolingFactor = Program.InitialCoolingFactor;
            section.FineCoolingFactor = Program.FineCoolingFactor;
            section.MalphaMaxChange = Program.Malpha_MaxChange;
            section.IterationAtEnd = Program.Iterationatend;
            section.ConvergenceLoopLimit = Program.ConvergenceLoopLimit;
            section.FreezingLimit = Program.FreezingLimit;
            section.ConvergeIntermediateClusters = Program.ConvergeIntermediateClusters;
            section.DebugPrintOption = DAVectorUtility.DebugPrintOption;
            section.ConsoleDebugOutput = DAVectorUtility.ConsoleDebugOutput;
        }

        public static void WriteTimingFile(string fileName, TimeSpan duration, double HPDuration, int ThreadCount, int MPIperNodeCount, int NodeCount, int PointCount_Process, int PointCount_Global, int maxNcent, double PTsplitting, double Ratio1, double MPIReduce, double Ratio2, double MPISRBasic, double Ratio3, double MPISREigen, double Ratio4, double MPIBdcast, double Ratio5, string DataFileName, DateTime CurrentTime, string ProcessorName)
        {
            string header = String.Format("Duration(ms)\tHPDuration(us)\tThread#\tMPIperNode\tNode#\tPattern\tParallelism\tPoint#/Process\tPoint#/Global\tmaxNcent\tPTsplitting\tRatio1\tMPIReduce\tRatio2\tMPISRBasic\tRatio3\tMPISREigen\tRatio4\tMPIBdcast\tRatio5\tDataFileName\tCurrentTime\tProcessorName");
            string format = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}";
            string line = string.Format(format, Math.Round(duration.TotalMilliseconds), // Total time
            Math.Round(HPDuration, 0),                                                  // High performance timer
            ThreadCount,                                                                // Thread# per process
            MPIperNodeCount,                                                            // Process# per node
            NodeCount,                                                                  // Node# aquired
            string.Format("{0}x{1}x{2}", ThreadCount, MPIperNodeCount, NodeCount),      // Pattern
            (ThreadCount * MPIperNodeCount * NodeCount).ToString(), PointCount_Process, // Local points
            PointCount_Global,                                                          // Global points
            maxNcent,                                                                   // cluster#
            Math.Round(PTsplitting, 0),                                                 // Partial times splitting
            Math.Round(Ratio1, 4),                                                      // PTsplitting vs HPDuration
            Math.Round(MPIReduce, 0),                                                   // MPI reduce
            Math.Round(Ratio2, 4),                                                      // MPIReduce vs HPDuration
            Math.Round(MPISRBasic, 0),                                                  // MPI SR Basic
            Math.Round(Ratio3, 4),                                                      // MPISRBasic vs HPDuration
            Math.Round(MPISREigen, 0),                                                  // MPI SR Eigen
            Math.Round(Ratio4, 4),                                                      // MPISREigen vs HPDuration
            Math.Round(MPIBdcast, 0),                                                   // MPI broadcast
            Math.Round(Ratio5, 4),                                                      // MPIBdcast vs HPDuration
            DataFileName,                                                               // Name of input data file
            CurrentTime, ProcessorName);

            if (!File.Exists(fileName))
            {
                File.AppendAllText(fileName, header);
                File.AppendAllText(fileName, System.Environment.NewLine + line);
            }
            else
            {
                File.AppendAllText(fileName, System.Environment.NewLine + line);
            }
        }

        public static ParallelOptions ParallelOptions
        {
            get
            {
                return _parallelOptions;
            }
        }

        public static void SetupDA()
        {
            if (!Program.DoBasicDA)
                return;
            
            //  Avoid Distributed Execution
            Program.TemperatureLimitforDistribution = -1.0;
            Program.ClusterLimitforDistribution = -1;

            Program.RW3DData = 3;   // Read and Write Plotviz

            Program.CalculateEigenvaluesfromMatrix = false;
            Program.UseSponge = false;
            Program.ContinuousClustering = true;
            Program.SigmaMethod = 0;

            Program.maxNcentperNode = 140;
            Program.maxNcentperNode = 10;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 8;
            Program.InitialNcent = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
            Program.maxNcentCreated = 500 * Program.maxNcentperNode;
            
            Program.NumberDataPoints = 200000;
            Program.StartPointPositiononInputLine = 0;
            Program.ParameterVectorDimension = 74;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;

            Program.MagicTemperatures[0] = -1.0;
            Program.MagicTemperatures[1] = -1.0;
            Program.MagicTemperatures[2] = -1.0;
            Program.MagicTemperatures[3] = -1.0;
            Program.MagicTemperatures[4] = -1.0;

            Program.Iterationatend = 800;
            Program.ConvergenceLoopLimit = 400;
            Program.Malpha_MaxChange = 0.005;
            Program.Malpha_MaxChange1 = 0.0005;
            Program.YChangeSquared = 0.00000001;
            Program.ExpArgumentCut1 = 20.0;
            Program.ExpArgumentCut2 = 40.0;
            Program.ExpArgumentCut2 = 20.0;
            Program.Waititerations = 1;
            Program.Waititerations_Converge = 4;
            Program.MinimumCountforCluster_Points = -1;
            Program.InitialCoolingFactor1 = 0.95;
            Program.FineCoolingFactor1 = 0.995;
            Program.InitialCoolingFactor2 = 0.95;
            Program.FineCoolingFactor2 = 0.995;
            Program.MinimumScaledWidthsquaredtosplit = -1.0;
            Program.ScaledSquaredDistanceatClosenessTest = 0.01;
            Program.TemperatureforClosenessTest = 0.01;
            Program.MinimumCountforCluster_C_k = 8.0; // Remove clusters with fewer points than this (based on C_k_) (This only used for converged clusters)
            Program.MinimumCountforCluster_Points = 10; // Remove clusters with fewer points than this (based on assigned points)

            Program.MaxNumberSplitClusters = 16;
            Program.ToosmalltoSplit = 200.0;
            Program.eigenvaluechange = 0.001;   // Limit 1 on Eigenvalue Changes
            Program.eigenvectorchange = 0.001;   // Limit 2 on Eigenvalue Changes

            Program.UseTriangleInequality_DA = 1;
            Program.OldCenterOption_DA = -1;
            Program.maxNcentTOTALforParallelism_DA = 20;   // Use Center parallelism when this limit hit
            Program.TriangleInequality_Delta1_old_DA = 0.1;  // Test for center change and old lower bounds (normalized by radius)
            Program.TriangleInequality_Delta1_current_DA = 0.1;  // Test for Center change and current lower bounds (normalized by radius)
            Program.MaxClusterLBsperPoint_DA = 140;   // Maximum number of Lower Bound values
            Program.MaxCentersperCenter_DA = 200;   // Maximum number of Centers in Center Difference Array

            Program.Tminimum = -20000.0;

            Program.PrintInterval = 5;
            Program.ClusterPrintNumber = 10;
            DAVectorUtility.DebugPrintOption = 2;

            /* 85399 54D
            Program.maxNcentperNode = 100;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 4;
            Program.InitialNcent = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
            Program.maxNcentCreated = 500 * Program.maxNcentperNode;

            Program.NumberDataPoints = 85399;
            Program.StartPointPositiononInputLine = 0;
            Program.ParameterVectorDimension = 54;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;
            Program.MaxNumberSplitClusters = 4;
            85399 54D */


            /*  Cmeans FLAME Test
            Program.MaxNumberSplitClusters = 1;
            Program.UseTriangleInequality_DA = 1;

            Program.maxNcentperNode = 5;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 5;
            Program.InitialNcent = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
            Program.maxNcentCreated = 500 * Program.maxNcentperNode;

            Program.NumberDataPoints = 22014;
            Program.StartPointPositiononInputLine = 1;
            Program.ParameterVectorDimension = 4;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;

            Program.UseTriangleInequality_DA = 0;
            */
            /* Lung
            Program.MaxNumberSplitClusters = 1;
            Program.UseTriangleInequality_DA = 1;

            Program.maxNcentperNode = 20;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 5;
            Program.InitialNcent = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
            Program.maxNcentCreated = 500 * Program.maxNcentperNode;

            Program.NumberDataPoints = 20054;
            Program.StartPointPositiononInputLine = 0;
            Program.ParameterVectorDimension = 4;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;

            Program.UseTriangleInequality_DA = 1;
            Lung */ 

            //  1000 point Dating Run
            Program.MaxNumberSplitClusters = 4;
            Program.Waititerations = 4;
            Program.InitialCoolingFactor1 = 0.995;
            Program.FineCoolingFactor1 = 0.9995;
            Program.InitialCoolingFactor2 = 0.995;
            Program.FineCoolingFactor2 = 0.9995;
            Program.ToosmalltoSplit = 15.0;

            Program.maxNcentperNode = 30;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 30;
            Program.InitialNcent = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_DA = Program.maxNcentTOTAL;
            Program.maxNcentCreated = 500 * Program.maxNcentperNode;

            Program.NumberDataPoints = 1000;
            Program.StartPointPositiononInputLine = 1;
            Program.ParameterVectorDimension = 3;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;

            Program.UseTriangleInequality_DA = 0;

        }   // End SetupDA

        public static void SetupKmeans()
        {
            if (!Program.DoKmeans)
                return;

            Program.RW3DData = 3;   // Read and Write Plotviz

            //  Avoid Distributed Execution
            Program.TemperatureLimitforDistribution = -1.0;
            Program.ClusterLimitforDistribution = -1;
            Program.SelectedInputLabel = 0; 

            Program.CalculateEigenvaluesfromMatrix = false;
            Program.UseSponge = false;
            Program.ContinuousClustering = false;
            Program.SigmaMethod = 0;
            Program.targetNcentperPoint = 124;

            Program.maxNcentperNode = 140;
            Program.maxNcentperNode = 10;
            Program.InitialNcent = Program.maxNcentperNode;
            Program.targetNcentperPoint = Program.maxNcentperNode;
            Program.maxNcentperPoint = Program.maxNcentperNode;
            Program.targetMinimumNcentperPoint = 8;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentCreated = 200;
            Program.maxNcentTOTALforParallelism_Kmeans = 20;

            Program.UseTriangleInequality_Kmeans = 0;  // 0 is Pure K means
            Program.UseTriangleInequality_Kmeans = 1;  // 0 is Pure K means

            Program.MaxClusterLBsperPoint_Kmeans = Program.maxNcentperNode;
            Program.MaxCentersperCenter_Kmeans = Program.maxNcentperNode;

            Program.NumberDataPoints = 200000;
            Program.ClusterIndexonInputLine = 75;
            Program.ClusterIndexonInputLine = -1;
            Program.ParameterVectorDimension = 74;
            Program.StartPointPositiononInputLine = 0;
            Program.FirstClusterValue = 0;

            Program.KmeansCenterChangeStop = 0.00001;
            Program.KmeansIterationLimit = 1000;
            Program.TriangleInequality_Delta1_old_Kmeans = 0.2;
            Program.TriangleInequality_Delta1_current_Kmeans = 0.2;

            // Dating-1000-1
            Program.ContinuousClustering = true;
            InitialNcent = 1;
            maxNcentperNode = 3;
            targetNcentperPoint = maxNcentperNode;
            maxNcentperPoint = maxNcentperNode;
            UseTriangleInequality_Kmeans = 0; // 0 is Pure K means
            NumberDataPoints = 1000;
            ParameterVectorDimension = 3;

            // Small Crandall Data
            /*Program.RW3DData = -1;
            Program.ParameterVectorDimension = 2048;
            Program.StartPointPositiononInputLine = 4;
            Program.ClusterCountOutput = -1;    // No output

            Program.SelectedInputLabel = 0;
            Program.InputFileType = 0;
            Program.Replicate = 1;
            Program.FirstClusterValue = -1;

            Program.NumberDataPoints = 76800;

            Program.maxNcentperNode = 3200;
            Program.InitialNcent = Program.maxNcentperNode;
            Program.targetNcentperPoint = 1;
            Program.maxNcentperPoint = 1;
            Program.targetMinimumNcentperPoint = 1;
            Program.maxNcentTOTAL = Program.maxNcentperNode;
            Program.maxNcentCreated = Program.maxNcentperNode;
            Program.maxNcentTOTALforParallelism_Kmeans = 20;

            Program.UseTriangleInequality_Kmeans = 0;  // 0 is Pure K means
            Program.ClusterIndexonInputLine = -3;
            Program.OldCenterOption_Kmeans = 100;
            Program.DoBackwardFacingTests_Kmeans = false;*/

            int CutBounds1 = 8;
            int CutBounds2 = 8;
            Program.MaxClusterLBsperPoint_Kmeans = Program.maxNcentperNode/CutBounds1;
            Program.MaxCentersperCenter_Kmeans = Program.maxNcentperNode/CutBounds2;

            /*  Cmeans/Flame
            Program.NumberDataPoints = 22014;
            Program.StartPointPositiononInputLine = 1;
            Program.ParameterVectorDimension = 4;
            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = 5;
            Program.InputFileType = 0;
            Program.Replicate = 1;
            Program.FirstClusterValue = 1;


            Program.InitialNcent = 5;
            Program.targetNcentperPoint = Program.InitialNcent;
            Program.maxNcentperNode = Program.InitialNcent;
            Program.maxNcentTOTAL = Program.InitialNcent;
            Program.maxNcentTOTALforParallelism_Kmeans = Program.InitialNcent + 1;
            Program.maxNcentCreated = Program.InitialNcent;

            Program.UseTriangleInequality_Kmeans = 0;
            */
            /* Lung Data

            Program.InitialNcent = 100;
            Program.targetNcentperPoint = Program.InitialNcent;
            Program.maxNcentperNode = Program.InitialNcent;
            Program.maxNcentTOTAL = Program.InitialNcent;
            Program.maxNcentTOTALforParallelism_Kmeans = Program.InitialNcent + 1;
            Program.maxNcentCreated = Program.InitialNcent;

            Program.NumberDataPoints = 85399;
            Program.StartPointPositiononInputLine = 0;
            Program.ParameterVectorDimension = 54;

            Program.SelectedInputLabel = 0;
            Program.ClusterIndexonInputLine = -1;
            Program.InputFileType = 0;
            Program.Replicate = 1;
            Program.MaxNumberSplitClusters = 4;
            Program.RW3DData = 3;   // Read and Write Plotviz
            */

            /* Read Plot file
            Program.ParameterVectorDimension = 3;
            Program.StartPointPositiononInputLine = 1;
            */

        }   // End SetupKmeans()

        public static void SetupLCMS()
        {   // Set up LC MS 2D analysis

            if (!Program.DoLCMS)
                return;

            int LCMSmode = 1;
  
            //  Generic Set up
            LCMSAnalyze.InitializeLCMS();

            //  Do a Harvard Format Full Analysis
            if (LCMSmode == 1)
            {
                LCMSAnalyze.FreshFullAnalysis();
                return;
            }

            //  Analyse Sponge Data from a Previous Run
            if (LCMSmode == 2)
            {
                LCMSAnalyze.SpongePointAnalysis();
                return;
            }

            //  Join and Equilibriate 2 Files together (original Full Analysis and Analaysis of its Sponge)
            //  Or just equilibriate a single file
            if (LCMSmode == 3)
            {
                LCMSAnalyze.JoinEquilibriate();
                return;
            }

            //  Change Sponge Factor
            if (LCMSmode == 4)
            {
                LCMSAnalyze.ChangeSpongeFactor();
                return;
            }

            //  Just Analyze a file
            if (LCMSmode == 5)
            {
                LCMSAnalyze.JustAnalyze();
                return;
            }


        }
    }
}
