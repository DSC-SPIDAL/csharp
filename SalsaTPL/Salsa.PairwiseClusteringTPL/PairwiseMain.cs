using System;
using System.Collections;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using HPC.Utilities;
using MPI;
using Salsa.Core.Configuration;
using Salsa.Core.Configuration.Sections;
using Salsa.Core.Blas;
using Salsa.Core;

using SALSALibrary;

#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using TDistance = System.Int16;
#else
using TDistance = System.Double;
#endif

//	Class Program  ********************************************************

//	ConvergenceTest decides if to stop EM Computation based on change of epsilon
//	calculated in Dist.calculateEpsi

//	shouldSplit invokes vectorclass.getMinimumEigenvalue to see if to split clusters
//	and which one to split

//	split generates an extra cluster and sets initial Malpha(k) -- not epsilonalpha(k)

//	initializeDist.Temperature calculates the maximum temperature from average distance

//	ReadDataFromFile reads a symmetric distance matrix into memory

//	getArrayindexFromMatrix is a general utility to find location in memory of a particular distance

//	Main reads arguments, sets initial conditions, calls ReadDataFromFile and then getDist
//	Finally it outputs results

//	Below t runs over threads, alpha runs over points, k runs over clusters
//	Dist.oldepsi[alpha][k] 		Old value of epsilon from previous iteration
//	Dist.epsi[alpha][k]		Value of epsilon calculated on current iteration
//	Dist.epsidiff[k] 	Change of epsilon on this iteration used for convergence
//	Dist.Malpha_k_[alpha][k]	Value of Malpha[k] at start of iteration
//	Dist.partialsum_C_k_[t][k] 	Contribution to C[k] from thread t
//	Dist.partialsum_A_k_[t][k] 	Contribution to A[k] from thread t
//	Dist.A_k_[k]			Value of A[k]
//	Dist.Balpha_k_[alpha][k] 	Value of Balpha[k]
//	Dist.C_k_[k] 			Value of C[k]

//	vectorclass.initialvector[alpha]		Initial vector in data point space for Power method (same vector used for each cluster)
//	vectorclass.Ax[k][alpha]			Full iteration vector (in data point space) for cluster k
//	vectorclass.oldAx[k][alpha]			Previous full iteration vector (in data point space) for cluster k
//	vectorclass.eigens[k]				Current eigenvalue per cluster k for powermethod
//	vectorclass.maxeigens[k]			Maximum eigenvalue per cluster k for powermethod
//	vectorclass.eigenconverged[k]			Indicator of status of eigenvalue determination per cluster k

//	Dist.Ncent			The current number of clusters initialized to 2
//	Program.ClustertoSplit			The cluster number that will be picked up to be split
//	Program.T			The current temperature
//	Program.iter			Iteration number counting EM convergence steps
//  Program.oldepsiset = 0              =0 says epsi has no previous value		       

//	Program.lenindex[t]		Number of data points each thread t takes care of
//	Program.startindex[t]		Index of the first datapoint that thread t takes care
//	DistanceMatrix[Program.getArrayindexFromMatrix[ClusterCenter,ClusterIndex]] is distance between ClusterCenter and ClusterIndex (DistanceMatrix passed as argument)

//	Program.MAX_VALUE		maximum number of parallel threads in a process (read as an argument)
//	PWCUtility.PointCount_Global		Number of Points (read as an argument) summed over threads and processes
//	Program.dataFile		Input data file name (read as an argument)
//	Program.maxNcent		Maximum number of cluster centers (read as an argument)
//	PWCUtility.pattern			String Pattern of parallel execution (e.g. 8x1x8 indicates [threads/process][processes/node][nodes])
//					(read as an argument)
//	Program.splitorexpandit		Execution parameter, scale factor of input data seta (read as an argument)

//	Program.outputFile		Output (Excel) file name (preset)


//	Program.cachelinesize = 32	Fixed definition of Cache line size in units of double
//	Program.max_change = 0.001	Convergence test condition for change in epsilon

//	Program.Tmin			Minimal temperature (preset)
//	Program.Tinitial		Initial (largest) temperature calculated in initializeDist.Temperature
//	Program.CoolingFactor = 0.95	Cooling Factor in Annealing (preset)
//	Program.eigenvaluechange = .001		Convergence test in power method for eigenvalues
//	Program.eigenvectorchange = .001	Convergence test in power method for eigenvector
//	Program.poweriterations = 50		Maximum iteration count in power method

//	MPI version
//	We arrange points in any order and assign them to processes/threads by simple one dimension decomposition
//	Points assigned to each process p run from PointIndex_Global1[p] to PointIndex_Global2[p] with PointTotal[p] points in all
//	PointTotal[p] = PointIndex_Global2[p] - PointIndex_Global1[p] + 1


//	A) First Calculate Malpha(k) in parallel over alpha (unless we have this preset either for initial iteration or for split cluster)
//	B) Then Calculate C(k) as sums over alpha of Malpha(k). This is straightforwardly parallel over alpha with final reduction

//	C) Then we need to calculate Balpha(k) as sums over beta (all Data points). This is dominant computation -- parallel over beta
//	It needs Malpha(k) Mbeta(k) Distance between alpha and beta, it calculates Balpha(k)and Bbeta(k)
//
//
//	D) Then we calculate A(k) by summing weighted Balpha(k) in parallel over alpha
//	E) This allows epsilon-alpha(k) and oldepsilon-alpha(k) to be calculated
//
//	F) Now we check for convergence. This is a simple parallel reduction operation
//	If not converged, repeat from A)
//	G) If converged, check if possible to split (we have Tmin and maximum cluster counts)
//	If possible to split, calculate eigenvalues by power method
//	J) If should split, reset cluster numbers and Malpha(k) and go to A)
//	If split not needed, then reduce temperaturte and go to A)

//	Calculation of critical eigenvalues is done in two almost identical steps considering Lmax maximum of C and then maximum of Lmax-C
//	These are done sepately for each cluster

//	Power Method algorithm 
//	This needs to iterate a computation Matrix * powervector where 
//	powervector-alpha(k) is calculated in parallel over alpha similarly to calculation of Balpha(k)
//	H) It needs A(k), Malpha(k) Mbeta(k) Distance between alpha and beta, Balpha(k), Bbeta(k)
//	I) This is followed by a simple convergence test for power vector

//	We could chose approach below which minimizes memory use but probably simpler approach that stores full matrix and not triangular part of it
//	allows better use of cache


//	Storage of distances reflects order of processing terms that depend on two data points ClusterCenter,ClusterIndex. ClusterCenter,ClusterIndex components are calculated in home of ClusterCenter
//	For each data point ClusterCenter, there are indices index-nn defined into CheckeredDistance[index-nn]
//	CD-FirstPoint1[ClusterCenter] is first data point label ClusterIndex of type 1 associated with ClusterCenter and it is stored in CD-index1[ClusterCenter]'th Position in CheckeredDistance
//	CD-FirstPoint2[ClusterCenter] is first data point label ClusterIndex of type 2 associated with ClusterCenter and it is stored in CD-index2[ClusterCenter]'th Position in CheckeredDistance
//	CD-LastPoint2[ClusterCenter] is last data point label ClusterIndex of type 2 associated with ClusterCenter
//	CD-NumPoints[ClusterCenter] is number (type 1 and 2) of points ClusterIndex associated with ClusterCenter in CheckeredDistance
//	values of ClusterIndex start at CD-FirstPoint1[ClusterCenter] and increment by 2 as long as ClusterIndex < CD-FirstPoint2[ClusterCenter]
//	Then they start again at CD-FirstPoint2[ClusterCenter] and increment by 2 as long as ClusterIndex <= CD-LastPoint2[ClusterCenter]


//	Block Computation is done for steps C) and H). This does parallelism minimizing memory access to Distance Matrix and b
//	Block is arranged to enhance communication and cache use
//	Other steps are straightforwardly parallel
//	There are two blocks -- one associated with Home index alpha and one with Away index beta
//	In C) and H), each Block has Malpha(k) and Balpha(k)

#region Deprecated
// Run Metadata
//public class RunMetadata
//{
//    public string DataFileName = " "; // Input data file name
//    public string TimingOutputFileName = " "; // Timing Output file name
//    public string BaseResultDirectoryName = "";  // Base Directory holding all Runs
//    public string ResultDirectoryExtension = ""; // 
//    public string ControlDirectoryName = ""; // Directory holding Control Information for this RunSet
//    public string DataLabelsFileName = "";   // File with application specific labels for this data
//    public string SMACOFInputFileName = "";  // File in ResultDirectoryName holding SMACOF input
//    public string RunSetLabel = "";  //  Label for RunSet
//    public int RunNumber = 0;    // Unique number of run
//    public int DataPoints = 0;  // Number of Data Points
//    public int ProcessingOption = 0; // Control Requested Processing
//    public string pattern = ""; // Parallel Pattern
//    public int splitorexpandit = 1; // Split or Expand Option
//    public int maxNcent = 1;    // Maximum Number of Centers
//    public int ThreadCount = 1;    // Number of Threads
//    public int NodeCount = 1;    // Number of Nodes
//    public int MPIperNodeCount = 1;    // Number of MPI Processes per Node
//    public int MPIIOStrategy = 0;   // Strategy for I/O in MPI
//    public int HistogramBinCount = 100; // Bin Count for Histograms
//    public string OutputDataFileName = " "; // Output data file name for Processing Option > 0
//    public string InputDataFileName = " "; // Input data file name for Processing Option > 0 (except 4)
//    public string InputClusteredFileName = " "; // Input data file name for Processing Option = 4
//    public string PointPropertiesFileName = " ";
//    public string FamilySMACOFFileName = " ";
//    public string ClusterSMACOFFileName = " ";
//    public string FamilySpecification = "ALULabel+consensus";
//    public string InputFormatStyle = "JustDistances";
//    public double RandomInputCut = 1.1;
//    public int GroupCenterSize = 10;
//    public string SelectedClusters = ""; // Comma separated list of Clusters
//    public double ToosmalltoSplit = 50.0;   // Size of Cluster Too small to split
//    public int Waititerations = 10; // Iterations to wait between splits
//    public double MinEigtest = -0.01;   // Fractional Eigenvalue Test
//    public bool ConvergeIntermediateClusters = false;   // If true converge all intermediate cluster
//    public double MinimumFeatureOverlap = 50.0;   // Cutoff for applying fancy formula
//    public int DistanceUndefinedOption = 0; // Control what to do when distance undefined
//    public int DistanceInfiniteOption = 0; // Control what to do when distance infinite
//    public int DistanceChoice = 0; // Control Distance calculation
//    public string Comment = ""; // User Comment
//    public int TransformDimension = 4;  // In reduction use this dimension (only 2 or 4 implemented)
//    public string Extradata1 = "";  // Use in special ways
//    public string Extradata2 = "";  // Use in special ways
//    public string Extradata3 = "";  // Use in special ways
//    public string Extradata4 = "";  // Use in special ways
//    public int ExtraOption1 = 0;    // use in various ways
//    public int DebugPrintOption = 1; // Control Debug Printing (= 0 None, = 1 Full, ==2 Summary)
//    public bool ConsoleDebugOutput = false; // If true send debug output to Console


//    public static void Membercopy(ref RunMetadata One, ref RunMetadata Two)
//    {
//        Two = (Salsa.PairwiseClustering.RunMetadata)One.MemberwiseClone();
//        return;
//    }

//}
#endregion

namespace Salsa.PairwiseClusteringTPL
{
    class Program
    {
        #region Members

        public static int ClusterCountOutput = 1;   // Control Label Output = 0 not at all, = 1 at each count, = 2 also find centers
        public static int[] ClusterAssignments = null;  // This gives for all points their cluster assignments (set in OutputClusterLabels)

        public static int cachelinesize = 32;

        public static bool ReadPartialMatrix = false; // Distance calculations are stored as partial matrices.

        public static string ControlFileName = "";  // Control File Name
        public static int ProcessingOption = 0; // Processing Option
        public static int PrintInterval = 5;   // Output at this interval

        public static int FirstClusterNumber = 0;   // Number (0 or 1) of first cluster
        public static bool ContinuousClustering = false;    // If true use the Ken Rose Continuous Clustering
        public static int maxNcent = 0; // maximum number of cluster center
        public static int InitialNcent = 1; // Initial value of Ncent
        public static int minimumclustercount = 25; // absorb clusters with fewer points than this

        public static double InitialCoolingFactor = 0.95; // InitialCooling Factor in Annealing
        public static double FineCoolingFactor = 0.995; // Refined Cooling Factor in Annealing
        public static double ConvergingCoolingFactor = 0.98; // Refined Cooling Factor in final iterations
        public static int Waititerations = 10;   // Wait this number of Temperature iterations before splitting

        public static int Iterationatend = 1000;  // Finish up EM Loop with at most this number of iterations
        public static int ConvergenceLoopLimit = 50; // Limit on EM Convergence for each step at a given temperature
        public static int CountPMExtraIterations = 0;   // Count in Pairwise converging P-M for Continuous Clustering
        public static int EMIterationStepCountMax = -1; // Count Maximum Iterations per loop

        public static double Epsi_max_change = 0.001;    //converge test condition for change in epsi
        public static double FreezingLimit = 0.002; // In finish stop when all freezing measures are < FreezingLimit

        public static bool ConvergeIntermediateClusters = true; // If true converge intermediate cluster cases
        public static double ToosmalltoSplit = 50.0;    // Size of Cluster measured by C_k_ that should not be split

        public static int PowerIterationLimit = 200; // Maximum iteration count in power method 
        public static double MinEigtest = -0.01; // Factor of Pass 0 eigenvalue used to test negative full eigenvalue (can be negative).  Test easier to satisfy as MinEigtest gets bigger
        public static int EMlimit_Duplication = 200;  // Limit on EM loops in Duplication mode
        public static int Eigenvalue_Methodology = 3;   // Specify matrix whose second derivative is to be calculated
        // =0 Do not look at second derivative for splitting
        // =1 original simple method using unaltered second derivative matrix
        // =2  estimate of Duplicated second derivative matrix valid for continuousclustering true
        // =3 EM estimate of Duplicated second derivative matrix
        // =4 Sophisticated estimate of Duplicated second derivative matrix (used in test only as 3 faster and same)
        public static double eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
        public static double eigenvectorchange = 0.001; // convergence test in power method for eigenvector
        public static int MaxNumberSplitClusters = 3;   // System will split upto this number simultaneously
        public static int CountTotalPowerIterations = 0; // Count Power Iterations
        public static int CountTotalPowerErrors = 0; // Count Power Iterations

        public static double ExpectedChange = 0.0;  // Expected Change in Cluster Count in Perturbation at Split
        public static int TestExpectedChange = -1;  // Cluster Index for Expected Change in Cluster Count in Perturbation at Split; if = -1 Nothing Set
        public static int ExpectedChangeMethod = -1;    // Method for Expected Change in Cluster Count in Perturbation at Split; 0 Change in C; 1 Special Ncent=1 case; =2 size of split vector
        public static double PreviousC = 0.0;   // Original Cluster Count in Perturbation at Split
        public static double[] deltaCoverC = new double[5]; //  Array to average splitting changes
        public static int[] CountdeltaCoverC = new int[5]; //  Count to average splitting changes
        public static int Countmethodzero = 0;  // Count Method 0 Perturbations
        public static int Countmethodtwo = 0;   // Count Method 2 Perturbations

        public static bool PerformEigenTest = false; // Set special Eigenvalue Testing
        public static double Epsi_max_change_Duplication = 0.001;    //converge test condition for change in epsi for Duplication mode

        public static int PerturbationVehicle = 0; // = 0 Perturb Epsilon ; if 1 Perturb Malpha_k_
        public static double SplitPerturbationFactor = 1.0; // Normalization of Perturbation for Split

        public static int JiggleOption = 0;   //    Control Perturbation of System every JiggleOption iterations
        public static int JigglePerturbation = 1;   //    Control How Perturbation implemented
        //  JigglePerturbation = 1 do randomly
        //  JigglePerturbation = 2 do using lowest second derivative eigenvalue
        public static double JigglePerturbationFactor = 10.0; // Normalization of Perturbation for Jiggle

        // Config Settings
        private static ConfigurationMgr _configurationManager;
        private static ParallelOptions _parallelOptions;
        #endregion

        // ***** Main *******
        // Pairwise (non-vector) Parallel clustering based on 
        // Deterministic Annealing algorithm 
        static void Main(string[] args)
        {
            // Load the command line args into our helper class which allows us to name arguments
            Salsa.Core.Arguments pargs = new Salsa.Core.Arguments(args);
            pargs.Usage = "Usage: Salsa.PairwiseClustering.exe /configFile=<string> /nodeCount=<int> /threadCount=<int>";

            if (pargs.CheckRequired(new string[] { "configFile", "nodeCount", "threadCount" }) == false)
            {
                Console.WriteLine(pargs.Usage);
                return;
            }

            //  Read Metadata using this as source of other metadata
            _configurationManager = ConfigurationMgr.LoadConfiguration(pargs.GetValue<string>("configFile"), true);
            ReadControlFile(_configurationManager);
            Program.ControlFileName = pargs.GetValue<string>("configFile");
            PWCUtility.NodeCount = pargs.GetValue<int>("nodeCount");
            PWCUtility.ThreadCount = pargs.GetValue<int>("threadCount");

            //  Set up TPL
            _parallelOptions = new ParallelOptions();
            _parallelOptions.MaxDegreeOfParallelism = PWCUtility.ThreadCount;

            //  Set up Parallelism
            PWCParallelism.SetupParallelism(ref args);

            // Set up Decomposition of USED points
            PWCParallelism.SetParallelDecomposition();

            PWCUtility.PatternLabel = String.Format("==== mpi-cluster({0}) ==== Threads:{1} Clusters:{2} PointCount_Global:{3} ==== {4} ({5}) ==== ",
                    PWCUtility.ParallelPattern,
                    PWCUtility.ThreadCount.ToString(),
                    Program.maxNcent.ToString(),
                    PWCUtility.PointCount_Global,
                    _configurationManager.PairwiseSection.DistanceMatrixFile,
                    DateTime.Now);

            PWCUtility.SALSAPrint(0, "\n" + PWCUtility.PatternLabel);

            // Initial Processing Complete
            PWCUtility.MPI_communicator.Barrier(); // Make certain all processes have processed original data before writing updated

            //  read data into memory
            if (Program.ReadPartialMatrix)
            {
                string singleDistanceFile = _configurationManager.PairwiseSection.DistanceMatrixFile;
                string path = Path.GetDirectoryName(singleDistanceFile);
                string prefix = Path.GetFileNameWithoutExtension(singleDistanceFile);
                string ext = Path.GetExtension(singleDistanceFile);                    
                PWCParallelism.ReadDataFromFiles(path, prefix, ext);
            }
            else
            {
                PWCParallelism.ReadDataFromFile(_configurationManager.PairwiseSection.DistanceMatrixFile);
            }

            // Start real work
            Program.ClusterAssignments = new int[PWCUtility.PointCount_Global];    // Set whenever clusters output

            if (Program.ProcessingOption == 101)
            {   // Find centers of Previously determined Clusters

                int NumberofClusters = 0;
                Program.ReadClusterNumbers(PWCUtility.ClusterNumberfile, out NumberofClusters, Program.ClusterAssignments, 0, PWCUtility.PointCount_Global);
                
                PWCUtility.SALSAPrint(0, PWCUtility.PatternLabel);
                PWCUtility.SALSAPrint(0, "Labels File: " + _configurationManager.PairwiseSection.ClusterFile);
                string file = "CenterFile-M" + Program.maxNcent.ToString() + "-C" + NumberofClusters + ".txt";
                string directory1 = Path.GetDirectoryName(_configurationManager.PairwiseSection.ClusterFile);
                string CenterFileName = Path.Combine(directory1, file);
                string place = Path.GetDirectoryName(_configurationManager.PairwiseSection.DistanceMatrixFile);
                string MDSFileName = Path.Combine(place, PWCUtility.addMDSfile);
                string FullLabelFileName = PWCUtility.Labelfile;

                PWCUtility.LengthCut1 = 250;    // Not used as yet
                PWCUtility.LengthCut2 = 250;    // Cut on Second Sequence Length
                PWCUtility.MinimumDistanceCut = 0.001;
                PWCUtility.LinkCountinCenterFinding = 30;

                PWCUtility.SALSAPrint(0, "Find Centers for " + NumberofClusters + " Cluster Numbers read from " + PWCUtility.ClusterNumberfile
                    + " Labels from " + FullLabelFileName + " MDS Read from " + MDSFileName + "\nLength Cuts " + PWCUtility.LengthCut1.ToString()
                    + " " + PWCUtility.LengthCut2.ToString() + " Minimum Sequence Distance " + PWCUtility.MinimumDistanceCut.ToString("F3")
                    + " Minimum Link Count " +PWCUtility.LinkCountinCenterFinding.ToString() );

                //  Create Centers
                FindCenters.FindGroupCenters(Program.ClusterAssignments, NumberofClusters, 0, CenterFileName, MDSFileName, FullLabelFileName);
                // Produce pviz file
                PWCUtility.MPI_communicator.Barrier();
                if (PWCUtility.MPI_Rank == 0)
                {
                    string plotDescription = PWCUtility.PointCount_Global + "_points_into_" + NumberofClusters +
                                             "_clusters_with_centers";
                    string clusterNumberFile = _configurationManager.PairwiseSection.ClusterNumberFile;
                    PlotTools.CreatePlotWithCenters(CenterFileName, MDSFileName, clusterNumberFile,
                                                    PWCUtility.CenterPointsPerCenterTypeInOutput,
                                                    PWCUtility.CenterPlotFile, plotDescription);
                }
                PWCUtility.MPI_communicator.Barrier(); 

                if ((PWCUtility.MPIIOStrategy > 0) || (PWCUtility.MPI_Rank == 0))
                    PWCUtility.WriteResults_Cluster(_configurationManager.PairwiseSection.SummaryFile, PWCUtility.CosmicOutput);
                PWCParallelism.TearDownParallelism(); //  Finalize MPI
                return;
            }

            for (int CountDCoverC = 0; CountDCoverC < 5; CountDCoverC++)
            {
                Program.deltaCoverC[CountDCoverC] = 0.0;
                Program.CountdeltaCoverC[CountDCoverC] = 0;
            }

            // Actually Determine Clusters
            Program.Eigenvalue_Methodology = 3;
            if (Program.ContinuousClustering)
                Program.Eigenvalue_Methodology = 2;
            CountPMExtraIterations = 0;

            //  Set up Timing
            PWCUtility.InitializeTiming(7);
            PWCUtility.SetUpMPISubTimers(1, "");
            PWCUtility.SetUpSubTimer(0, "Splitting");

            //  Do Clustering
            Dist runitall = new Dist();
            runitall.getDist();

            //  End Timing
            PWCUtility.EndTiming();

            PWCUtility.SALSAPrint(0, "\n----------------------------- Final Output\n" + PWCUtility.PatternLabel); 
            PWCUtility.SALSAPrint(0, "Labels File: " + _configurationManager.PairwiseSection.ClusterFile);
            PWCUtility.SALSAPrint(0, "Timing Output: " + _configurationManager.PairwiseSection.TimingFile);
            PWCUtility.SALSAPrint(0, "Initial Number of Centers: " + Program.InitialNcent.ToString());
            PWCUtility.SALSAPrint(0,  "Maximum Number of Centers: " + Program.maxNcent.ToString());
            PWCUtility.SALSAPrint(0, "Actual Number of Centers: " + Dist.RunningPWC.Ncent.ToString());
            PWCUtility.SALSAPrint(0, "Minimum Number Points in Final Clusters " + Program.minimumclustercount.ToString());
            PWCUtility.SALSAPrint(0, "Converge Intermediate Clusters " + Program.ConvergeIntermediateClusters.ToString());
            PWCUtility.SALSAPrint(0,  "Initial Cooling Factor in Annealing: " + Program.InitialCoolingFactor.ToString("F4"));
            PWCUtility.SALSAPrint(0, "Refined Cooling Factor in Annealing: " + Program.FineCoolingFactor.ToString("F4"));
            PWCUtility.SALSAPrint(0, "Converging(Final) Cooling Factor in Annealing: " + Program.ConvergingCoolingFactor.ToString("F4"));
            PWCUtility.SALSAPrint(0, "Continuous Clustering " + Program.ContinuousClustering.ToString());
            PWCUtility.SALSAPrint(0,  "Eigenvalue Methodology: " + Program.Eigenvalue_Methodology.ToString());
            PWCUtility.SALSAPrint(0,  "Do not split Clusters smaller than this: " + Program.ToosmalltoSplit.ToString());
            PWCUtility.SALSAPrint(0,  "Pass 1 Eigenvalue Fractional Test: " + Program.MinEigtest.ToString());
            PWCUtility.SALSAPrint(0,  "Wait stages between splits: " + Program.Waititerations.ToString());
            PWCUtility.SALSAPrint(0, "Maximum Number of Simultaneous Cluster Splits " + Program.MaxNumberSplitClusters.ToString());
            PWCUtility.SALSAPrint(0,  "Jiggle and Split Perturbation Vehicle: " + Program.PerturbationVehicle.ToString());
            PWCUtility.SALSAPrint(0,  "Split and Split Perturbation Factor: " + Program.SplitPerturbationFactor.ToString());
            PWCUtility.SALSAPrint(0,  "Jiggle Option:" + Program.JiggleOption.ToString());
            PWCUtility.SALSAPrint(0,  Program.JigglePerturbation.ToString() + " Jiggle Perturbation Method");
            PWCUtility.SALSAPrint(0, Program.JigglePerturbationFactor.ToString() + " Jiggle Perturbation Factor");

            // Calculate Cluster Statistics
            int[] counts = new int[Dist.RunningPWC.Ncent];

            if (ClusterCountOutput > 0)
            {
                // Note Program.ClusterAssignments 
                Program.OutputClusterLabels(counts);
                if (ClusterCountOutput == 2)
                {
                    string file = "CenterFile-M" + Program.maxNcent.ToString() + "-C" + Dist.RunningPWC.Ncent.ToString() + ".txt";
                    string directory1 = Path.GetDirectoryName(_configurationManager.PairwiseSection.ClusterFile);
                    string CenterFileName = Path.Combine(directory1, file);
                    string place = Path.GetDirectoryName(_configurationManager.PairwiseSection.DistanceMatrixFile);
                    string MDSFileName = Path.Combine(place, PWCUtility.addMDSfile);
                    string FullLabelFileName = PWCUtility.Labelfile;
                    FindCenters.FindGroupCenters(Program.ClusterAssignments, Dist.RunningPWC.Ncent, 0, CenterFileName, MDSFileName, FullLabelFileName);
                }
                
            }
            else
            {
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                {
                    counts[ClusterIndex] = 0;
                }
            }

            PWCUtility.SALSAPrint(0, "\n******************\nT " + Dist.RunningPWC.Temperature.ToString("F4") + " Cluster " + Dist.RunningPWC.Ncent.ToString() + 
                " Iter " + Dist.EMIterationCount.ToString() + " Extra Iter " + Dist.Extra_EMIterationCount.ToString() + " Max Iterations per Step " + Program.EMIterationStepCountMax.ToString() );
            PWCUtility.SALSAPrint(0, "Extra P-M Iterations " + Program.CountPMExtraIterations.ToString() + " Total Number of Power Iterations "
                + Program.CountTotalPowerIterations.ToString() + " with errors " + Program.CountTotalPowerErrors.ToString());
            PWCUtility.SALSAPrint(0, "Iterations at End " + Program.Iterationatend.ToString() + " Iteration Limit to converge a fixed EM iteration " + Program.ConvergenceLoopLimit.ToString().ToString() 
                + " Count in Pairwise converging P-M for Continuous Clustering " + Program.CountPMExtraIterations.ToString() 
                + " Maximum Iterations per loop " + Program.EMIterationStepCountMax.ToString() ); 
            PWCUtility.SALSAPrint(0, "Freezing Convergence Limit "  + Program.FreezingLimit.ToString("F5")
                + " Convergence test condition for change in epsilon " + Program.Epsi_max_change.ToString("F5") );

            string nextline = " Center Averages (Counts) [Width] Width here and earlier is sum of d(i,j) weighted with M and divided by squared of occupation count C";
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                double meandist = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
                nextline += Dist.RunningPWC.C_k_[ClusterIndex].ToString("F2") + " (" + counts[ClusterIndex].ToString() + ") [" + meandist.ToString("F5") + "] ";
            }
            PWCUtility.SALSAPrint(0, nextline);

            // Calculate Distance Correlation Matrix
            runitall.CorrelationCalculation();
            PWCUtility.SALSAPrint(0, "\nDistance Correlation Function");
            double successmeasure1 = 0.0;
            double AverageFreeze = 0.0;
            double MaxFreeze = 0.0;
            for (int ClusterIndex1 = 0; ClusterIndex1 < Dist.RunningPWC.Ncent; ClusterIndex1++)
            {
                AverageFreeze += Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex1];
                MaxFreeze = Math.Max(MaxFreeze, Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex1]);
                
                string anotherline = (ClusterIndex1 + Program.FirstClusterNumber).ToString();
                for (int ClusterIndex2 = 0; ClusterIndex2 < Dist.RunningPWC.Ncent; ClusterIndex2++)
                {
                    double tmp = Dist.DistanceCorrelation[ClusterIndex1, ClusterIndex2] / Math.Sqrt(Dist.DistanceCorrelation[ClusterIndex1, ClusterIndex1] * Dist.DistanceCorrelation[ClusterIndex2, ClusterIndex2]);
                    anotherline += "\t" + tmp.ToString("F2");
                    if (ClusterIndex1 != ClusterIndex2)
                    {
                        successmeasure1 += tmp;
                    }

                }
                anotherline += "   Count " + counts[ClusterIndex1].ToString();
                PWCUtility.SALSAPrint(0, anotherline);
            }

            successmeasure1 = successmeasure1 / ((Dist.RunningPWC.Ncent - 1) * Dist.RunningPWC.Ncent);
            PWCUtility.SALSAPrint(0, "\nOff Diagonal Average " + successmeasure1.ToString("F2")); 
            
            AverageFreeze = AverageFreeze / Dist.RunningPWC.Ncent;
            PWCUtility.SALSAPrint(0, "\nAverage Freezing Coefficient " + AverageFreeze.ToString("F6") + " Maximum Freezing Coefficient " + MaxFreeze.ToString("F6"));

            //  Output Statistics on Perturbations at Cluster Splitting
            PWCUtility.SALSAPrint(0, "\nCluster Splitting Perturbation Statistics    Method 0: " + Program.Countmethodzero.ToString() + " Method 2: " + Program.Countmethodtwo.ToString());
            if (Program.CountdeltaCoverC[0] > 0) 
                PWCUtility.SALSAPrint(0, "Actual Delta C over C Method 1 (Ncent =1)    " + (Program.deltaCoverC[0] / Program.CountdeltaCoverC[0]).ToString("F6") + " Count " + Program.CountdeltaCoverC[0].ToString());
            if (Program.CountdeltaCoverC[2] > 0)
                PWCUtility.SALSAPrint(0, "Well Predicted Delta C/C Method 1 (Ncent =1) " + (Program.deltaCoverC[2] / Program.CountdeltaCoverC[2]).ToString("F6") + " Count " + Program.CountdeltaCoverC[2].ToString());
            if (Program.CountdeltaCoverC[1] > 0)
                PWCUtility.SALSAPrint(0, "Actual Delta C over C Method 0 or 2          " + (Program.deltaCoverC[1] / Program.CountdeltaCoverC[1]).ToString("F6") + " Count " + Program.CountdeltaCoverC[1].ToString());
            if (Program.CountdeltaCoverC[3] > 0) 
                PWCUtility.SALSAPrint(0, "Well Predicted Delta C over C Method 0 or 2  " + (Program.deltaCoverC[3] / Program.CountdeltaCoverC[3]).ToString("F6") + " Count " + Program.CountdeltaCoverC[3].ToString());
            if (Program.CountdeltaCoverC[4] > 0) 
                PWCUtility.SALSAPrint(0, "Naive Predicted Delta C over C Method 0 or 2 " + (Program.deltaCoverC[4] / Program.CountdeltaCoverC[4]).ToString("F6") + " Count " + Program.CountdeltaCoverC[4].ToString());

            /* Compute the duration between the initial and the end time ignoring print out. */
            TimeSpan duration = PWCUtility.endTime - PWCUtility.startTime;
            PWCUtility.SALSAPrint(0,  "\nTotal Time excluding I/O  " + duration.ToString() + " " + (PWCUtility.HPDuration * .001).ToString("F0"));
            nextline = "Partial Times ";
            for (int itimer = 0; itimer < PWCUtility.NumberofSubTimings; itimer++)
            {
                if (PWCUtility.SubTimingEnable[itimer])
                {
                    double tmp = PWCUtility.SubDurations[itimer] / PWCUtility.HPDuration;
                    nextline += PWCUtility.SubTimingNames[itimer] + " " + (PWCUtility.SubDurations[itimer] * .001).ToString("F0") + " " + tmp.ToString("F4") + " ";
                }
            }
            PWCUtility.SALSAPrint(0, nextline);

            if ((PWCUtility.MPIIOStrategy > 0) || (PWCUtility.MPI_Rank == 0))
            {
                PWCUtility.WriteResults_Cluster(_configurationManager.PairwiseSection.SummaryFile, PWCUtility.CosmicOutput);

                WriteTimingFile(_configurationManager.PairwiseSection.TimingFile,
                    duration,
                    PWCUtility.HPDuration,
                    PWCUtility.ThreadCount,
                    PWCUtility.MPIperNodeCount,
                    PWCUtility.NodeCount,
                    PWCUtility.PointCount_Process,
                    PWCUtility.PointCount_Global,
                    Program.maxNcent,
                    PWCUtility.SubDurations[0] * 0.001,  // convert to milliseconds
                    PWCUtility.SubDurations[0] / PWCUtility.HPDuration,
                    PWCUtility.SubDurations[1] * 0.001,  // convert to milliseconds
                    PWCUtility.SubDurations[1] / PWCUtility.HPDuration,
                    PWCUtility.SubDurations[2] * 0.001,  // convert to milliseconds
                    PWCUtility.SubDurations[2] / PWCUtility.HPDuration,
                    PWCUtility.SubDurations[3] * 0.001,  // convert to milliseconds
                    PWCUtility.SubDurations[3] / PWCUtility.HPDuration,
                    PWCUtility.SubDurations[4] * 0.001,  // convert to milliseconds
                    PWCUtility.SubDurations[4] / PWCUtility.HPDuration,
                    _configurationManager.PairwiseSection.DistanceMatrixFile,
                    DateTime.Now,
                    MPI.Environment.ProcessorName);
            }


            PWCParallelism.TearDownParallelism(); //  Finalize MPI

            return;

        }
        // End Main

 
        public static void OutputClusterLabels(int[] OccupationCounts)
        {
            // Generate Cluster Labels (cluster number) from final epsi
            int[] labels = new int[PWCUtility.PointCount_Process];
            int[][] partialsum_OccupationCounts = new int[PWCUtility.ThreadCount][];

            for (int ThreadNo = 0; ThreadNo < PWCUtility.ThreadCount; ThreadNo++)
            {
                partialsum_OccupationCounts[ThreadNo] = new int[Dist.RunningPWC.Ncent + cachelinesize];
            }

            //  Parallel Section setting cluster labels
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                Array.Clear(partialsum_OccupationCounts[ThreadNo], 0, Dist.RunningPWC.Ncent);
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;

                for (int index = beginpoint; index < indexlen + beginpoint; index++)
                {
                    double distmin = 9999999999999.0;
                    int knear = 0;
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        if (Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex] < distmin)
                        {
                            distmin = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                            knear = ClusterIndex;
                        }
                    }
                    labels[index] = knear + Program.FirstClusterNumber;
                    partialsum_OccupationCounts[ThreadNo][knear]++;
                }
            });  // End Parallel Section setting cluster labels


            int[] LocalOccupationCounts = new int[Dist.RunningPWC.Ncent];
            Array.Clear(LocalOccupationCounts, 0, Dist.RunningPWC.Ncent);
            for (int ThreadNo = 0; ThreadNo < PWCUtility.ThreadCount; ThreadNo++)
            {
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                {
                    LocalOccupationCounts[ClusterIndex] += partialsum_OccupationCounts[ThreadNo][ClusterIndex];
                }
            }


            if (PWCUtility.MPI_Size > 1)
            {
                PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming);
                LocalOccupationCounts = PWCUtility.MPI_communicator.Allreduce<int>(LocalOccupationCounts, Operation<int>.Add);
                PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming);
            }

            for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
            {
                OccupationCounts[ClusterCount] = LocalOccupationCounts[ClusterCount];
            }
            PWCUtility.PreciseTimer.Stop();    // temporarily stop timer
            PWCUtility.HPDuration += PWCUtility.PreciseTimer.Duration;

            string directory = Path.GetDirectoryName(_configurationManager.PairwiseSection.ClusterFile);
            string file = Path.GetFileNameWithoutExtension(_configurationManager.PairwiseSection.ClusterFile) + "-M" + Program.maxNcent.ToString() + "-C" + Dist.RunningPWC.Ncent.ToString() + Path.GetExtension(_configurationManager.PairwiseSection.ClusterFile);

            string ClusternumberFileName = Path.Combine(directory, file);

            if (PWCUtility.MPIIOStrategy > 0)
            {
                WriteClusterFile(ClusternumberFileName, ref labels, PWCUtility.PointCount_Process, PWCUtility.PointStart_Process, false);
            }


            int MPItag = 100;
            if (PWCUtility.MPI_Rank == 0)
            {
                WriteClusterFile(ClusternumberFileName, ref labels, PWCUtility.PointCount_Process, PWCUtility.PointStart_Process, false);

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;

                    for (int index = beginpoint; index < indexlen + beginpoint; index++)
                    {
                        Program.ClusterAssignments[index + PWCUtility.PointStart_Process] = labels[index] - Program.FirstClusterNumber;
                    }
                });  // End Parallel Section setting cluster assignments in process 0

                for (int MPISource = 1; MPISource < PWCUtility.MPI_Size; MPISource++)
                {
                    MPIPacket<int> fromsource;
                    fromsource = PWCUtility.MPI_communicator.Receive<MPIPacket<int>>(MPISource, MPItag);
                    if (PWCUtility.MPIIOStrategy == 0)
                        WriteClusterFile(ClusternumberFileName, ref fromsource.Marray, fromsource.NumberofPoints, fromsource.FirstPoint, true);

                    for (int index = 0; index < fromsource.NumberofPoints; index++) 
                    {
                        Program.ClusterAssignments[index + fromsource.FirstPoint] = fromsource.Marray[index] - Program.FirstClusterNumber;
                    }
                }
            }   // End Root Process that receives cluster assignments from afar
            else
            {
                MPIPacket<int> tosend = new MPIPacket<int>(PWCUtility.PointCount_Process);
                tosend.FirstPoint = PWCUtility.PointStart_Process;
                tosend.NumberofPoints = PWCUtility.PointCount_Process;
                labels.CopyTo(tosend.Marray, 0);
                PWCUtility.MPI_communicator.Send<MPIPacket<int>>(tosend, 0, MPItag);
            }
            PWCUtility.MPI_communicator.Broadcast<int>(ref Program.ClusterAssignments, 0);
            PWCUtility.MPI_communicator.Barrier();
            PWCUtility.PreciseTimer.Start();   // Restart Timer
            return;
        }

        public static void ReadControlFile(ConfigurationMgr mgr)
        {
            PairwiseSection section = mgr.PairwiseSection;

            PWCUtility.PointCount_Global = section.DataPoints;
            Program.ProcessingOption = section.ProcessingOption;
            PWCUtility.ParallelPattern = section.Pattern;
            Program.maxNcent = section.MaxNcent;
            PWCUtility.ThreadCount = section.ThreadCount;     // Number of Threads
            PWCUtility.NodeCount = section.NodeCount;       // Number of Nodes
            PWCUtility.MPIperNodeCount = section.MPIperNodeCount; // Number of MPI processes per Node
            PWCUtility.MPIIOStrategy = section.MPIIOStrategy;   // Controls strategy of file handling with MPI =0 is ONE FILE
            Program.ToosmalltoSplit = section.ToosmalltoSplit;
            Program.Waititerations = section.Waititerations;
            Program.MinEigtest = section.MinEigtest;
            Program.Epsi_max_change = section.Epsi_max_change;     //converge test condition for change in epsi
            Program.InitialCoolingFactor = section.InitialCoolingFactor;  // InitialCooling Factor in Annealing
            Program.FineCoolingFactor = section.FineCoolingFactor;    // Refined Cooling Factor in Annealing
            Program.eigenvaluechange = section.Eigenvaluechange;    // convergence test in power method for eigenvalues
            Program.eigenvectorchange = section.Eigenvectorchange;   // convergence test in power method for eigenvector
            Program.Iterationatend = section.Iterationatend;       // Finish up EM Loop with this number of iterations
            Program.ConvergenceLoopLimit = section.ConvergenceLoopLimit; // Limit on EM Convergence
            Program.FreezingLimit = section.FreezingLimit;       // In finish stop when all freezing measures are < FreezingLimit
            Program.PowerIterationLimit = section.PowerIterationLimit;   // Maximum iteration count in power method
            Program.ConvergeIntermediateClusters = section.ConvergeIntermediateClusters;
            PWCUtility.DebugPrintOption = section.DebugPrintOption;
            PWCUtility.ConsoleDebugOutput = section.ConsoleDebugOutput;
            Program.ReadPartialMatrix = section.ReadPartialMatrix;
            Program.ContinuousClustering = section.ContinuousClustering;
            PWCUtility.Labelfile = section.LabelFile;
            PWCUtility.addMDSfile = section.AddMdsFile;
            PWCUtility.ClusterNumberfile = section.ClusterNumberFile;
            PWCUtility.addMDS = section.AddMds;
            PWCUtility.BucketFractions = !string.IsNullOrEmpty(section.BucketFractions) ? section.BucketFractions.Trim().Split(new[] { ',' }).Select(x => double.Parse(x)).ToArray() : null;
            PWCUtility.NumberofBuckets = PWCUtility.BucketFractions != null ? PWCUtility.BucketFractions.Length : 0;
            PWCUtility.NumberofCenters = section.NumberOfCenters;
            PWCUtility.CenterPointsPerCenterTypeInOutput = section.CenterPointsPerCenterTypeInOuput;
            PWCUtility.CenterPlotFile = section.CenterPlotFile;
        }

        public static void WriteControlFile(ConfigurationMgr mgr)
        {
            PairwiseSection section = mgr.PairwiseSection;

            section.DataPoints = PWCUtility.PointCount_Global;
            section.ProcessingOption = Program.ProcessingOption;
            section.Pattern = PWCUtility.ParallelPattern;
            section.MaxNcent = Program.maxNcent;
            section.ThreadCount = PWCUtility.ThreadCount;
            section.NodeCount = PWCUtility.NodeCount;
            section.MPIperNodeCount = PWCUtility.MPIperNodeCount;
            section.MPIIOStrategy = PWCUtility.MPIIOStrategy;
            section.ToosmalltoSplit = Program.ToosmalltoSplit;
            section.Waititerations = Program.Waititerations;
            section.MinEigtest = Program.MinEigtest;
            section.Epsi_max_change = Program.Epsi_max_change;
            section.InitialCoolingFactor = Program.InitialCoolingFactor;
            section.FineCoolingFactor = Program.FineCoolingFactor;
            section.Eigenvaluechange = Program.eigenvaluechange;
            section.Eigenvectorchange = Program.eigenvectorchange;
            section.Iterationatend = Program.Iterationatend;
            section.ConvergenceLoopLimit = Program.ConvergenceLoopLimit;
            section.FreezingLimit = Program.FreezingLimit;
            section.PowerIterationLimit = Program.PowerIterationLimit;
            section.ConvergeIntermediateClusters = Program.ConvergeIntermediateClusters;
            section.DebugPrintOption = PWCUtility.DebugPrintOption;
            section.ConsoleDebugOutput = PWCUtility.ConsoleDebugOutput;
            section.ReadPartialMatrix = Program.ReadPartialMatrix;
            section.ContinuousClustering = Program.ContinuousClustering;
            section.LabelFile = PWCUtility.Labelfile;
            section.AddMdsFile = PWCUtility.addMDSfile;
            section.ClusterNumberFile = PWCUtility.ClusterNumberfile;
            section.AddMds = PWCUtility.addMDS;
            section.BucketFractions = PWCUtility.BucketFractions != null
                                          ? string.Join(",", PWCUtility.BucketFractions)
                                          : string.Empty;
            section.NumberOfCenters = PWCUtility.NumberofCenters;
        }

        public static void WriteClusterFile(string fname, ref int[] labels, int dataPoints, int startposition, bool append)
        {
            FileMode mode = append ? FileMode.Append : FileMode.Create;

            using (FileStream stream = new FileStream(fname, mode, FileAccess.Write, FileShare.ReadWrite))
            {
                using (StreamWriter writer = new StreamWriter(stream))
                {
                    for (int i = 0; i < dataPoints; i++)
                    {
                        writer.WriteLine(String.Format("{0} {1}", i + startposition, labels[i]));
                    }
                }
            }
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

        // Read Cluster Numbers 
        public static void ReadClusterNumbers(string fname, out int NumberofClusters, int[] ClusterNumbers, int BeginPoint, int PointstoRead)
        {
            char[] _sep = new[] { ' ', '\t' };

            bool success = false;
            string line = "NotSet";
            int count = 0;
            NumberofClusters = 0;

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                    sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            int NumberPosition = 1;
                            if (splits.Length > 2)
                                NumberPosition = 4;
                            int index = int.Parse(splits[0]);
                            if (index < BeginPoint)
                                continue;
                            if (index >= BeginPoint + PointstoRead)
                                break;
                            int readnumber = int.Parse(splits[NumberPosition]);
                            ClusterNumbers[index - BeginPoint] = readnumber - Program.FirstClusterNumber;
                            NumberofClusters = Math.Max(NumberofClusters, readnumber);
                            if (readnumber < Program.FirstClusterNumber)
                            {
                                Exception e = PWCUtility.SALSAError("Illegal Cluster Number on Cluster Number / MDS file Index " + index.ToString()
                                    + " Cluster Number Read " + readnumber + " First Cluster # " + Program.FirstClusterNumber.ToString());
                                throw (e);
                            }
                            count++;
                        }
                    }
                    if (count != PointstoRead)
                    {
                        Exception e = PWCUtility.SALSAError("Illegal count on Cluster Number file " + count.ToString() + " " + PointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Cluster Number data: count " + count.ToString() + " Line " + line + "\n" + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = PWCUtility.SALSAError("Cluster Number File read error: count " + count.ToString() + " Line " + line + "\n" + fname);
                throw (e);
            }

            if (Program.FirstClusterNumber == 0)
                ++NumberofClusters;

        }   // End Read Cluster Numbers

        public static ParallelOptions ParallelOptions
        {
            get
            {
                return _parallelOptions;
            }
        }
    }
}
