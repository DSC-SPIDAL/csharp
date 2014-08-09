using System;
using System.Collections.Generic;
using System.IO;
using System.Threading.Tasks;
using MDS;
using SALSALibrary;
using Salsa.Core;
using Salsa.Core.Configuration;
using Salsa.Core.Configuration.Sections;
using Environment = MPI.Environment;
using ManxcatMDSRoutines = MDS.ManxcatMDS;
using ManxcatRotateMDSRoutines = RotateMDS.RotateManxcatMDS;

namespace Manxcat
{
    public class ManxcatCentral
    {
        //  Input Variables               
        public static ConfigurationMgr ConfMgr;
        public static ManxcatSection Configuration;
        public static string ConfigurationFileName = string.Empty; // Configuration file name
        public static int ProcessingOption = 0; // Processing Option

        //  Utility variables
        public static MPI2DDoubleVectorPacket TogoDistributed2DDoubleVector;
                                              // Parameter Data to be sent on MPI transfers

        public static MPI2DDoubleVectorPacket FromAfar2DDoubleVector; // Data received by MPI transfers

        public static MPI1DStringVectorPacket TogoDistributed1DStringVector;
                                              // Parameter Data to be sent on MPI transfers

        public static MPI1DStringVectorPacket FromAfar1DStringVector; // Data received by MPI transfers

        public static MPI2DDoubleVectorPacket TogoDiagVector; // diag Data to be sent on MPI transfers
        public static MPI2DDoubleVectorPacket TogoSqDgInvVector; // sqdiag Data to be sent on MPI transfers

        public static bool violat;
                           // error return for Calcfg (= false i,mplies parameter values changed as violated boundary conditions)

        public static double ChisqFunctionCalcMultiplier = 10.0;
                             // Multiplicative Constant for Function (square this for Chisq) Calculation)

        public static double ChisqPrintConstant = 1.0;
                             // Extrs Multiplicative Constant for Chisq Printing (divided by number of degrees of freedom)

        public static string ResultDirectoryName = ""; // Full String for Result Directory
        public static string ActualDataFileName = ""; // Full Input File name  for Data File Name
        public static string ActualOutputFileName = ""; // Full Output File namee

        //  Note Best solution alwasys stored with Scaling or Addition of Marquardt Factor if this requested
        public static bool BestSolutionSet = false; // If true there is a best solutiom
        public static bool BestxshiftSet = false; // If true there is a best xshift set
        public static bool GlobalParameterSet = false; // If True Global parameters set for Current Solution

        public static int IterationtorecalculateQLimits = 0;
                          // Iteration after and including which calculating Q limits warranted

        public static double LineFactorGuess = 1.0; // Initial Guess at Line Factor

        private static void Main(string[] args)
        {
            // Load the command line args into our helper class which allows us to name arguments
            var pargs = new Arguments(args);
            pargs.Usage = "Usage: Salsa.PairwiseClustering.exe /configFile=<string> /nodeCount=<int> /threadCount=<int>";

            if (pargs.CheckRequired(new[] {"configFile", "nodeCount", "threadCount"}) == false)
            {
                Console.WriteLine(pargs.Usage);
                return;
            }


            ConfigurationFileName = pargs.GetValue<string>("configFile"); // Name of file containing configuration data
            ReadConfiguration();

            if (args.Length >= 2)
            {
                SALSAUtility.NodeCount = pargs.GetValue<int>("nodeCount");
                Configuration.NodeCount = SALSAUtility.NodeCount;
            }

            SALSAUtility.ThreadCount = Configuration.ThreadCount; // Number of Threads

            if (args.Length == 3)
            {
                SALSAUtility.ThreadCount = pargs.GetValue<int>("threadCount");
                Configuration.ThreadCount = SALSAUtility.ThreadCount;
            }

            SALSAUtility.ConsoleDebugOutput = Configuration.ConsoleDebugOutput;
            SALSAUtility.DebugPrintOption = Configuration.DebugPrintOption;

            //  Set up Threading
            SALSAUtility.SetupParallelOptions();

            //  Set up Parallelism
            SALSAParallelism.SetupParallelism(ref args);

            // Define Points to be used
            SALSA_ProcessVariedandFixed.setupFixedandVaried();

            // Set up Decomposition of USED points
            SALSAParallelism.SetParallelDecomposition();

            //  Redistribute points so an equal number (upto 1) varied in each thread
            SALSA_ProcessVariedandFixed.RedistributePoints();

            try
            {
                // Set up Normalizations
                ChisqPrintConstant = Configuration.ChisqPrintConstant; // This can be customized in set up routines
                ChisqFunctionCalcMultiplier = Configuration.FunctionErrorCalcMultiplier;

                ProcessingOption = Configuration.ProcessingOption; // Determines type of fit

                //  Set up Hotsum and DoubleStar initial parameters
                if (ProcessingOption <= 100)
                    ManxcatMDSRoutines.SetupHotsunforMDS();

                else if (ProcessingOption <= 200)
                    ManxcatRotateMDSRoutines.SetupHotsunforRotateMDS();
                Hotsun.SetupManxcat();

                if (Hotsun.DecomposeParameters)
                {
                    // Set up MPI if parallel parameter
                    SetupMPIPacket(out FromAfar2DDoubleVector);
                    SetupMPIPacket(out TogoDistributed2DDoubleVector);
                    SetupMPIPacket(out FromAfar1DStringVector);
                    SetupMPIPacket(out TogoDistributed1DStringVector);

                    SetupMPIPacket(out TogoDiagVector);
                    SetupMPIPacket(out TogoSqDgInvVector);
                }

                //  Set up Timing
                SALSAUtility.InitializeTiming(8);
                SALSAUtility.SetUpMPISubTimers(3, "");
                SALSAUtility.SetUpSubTimer(0, "MatrixSolve");
                SALSAUtility.SetUpSubTimer(1, "EigenSolve");
                SALSAUtility.SetUpSubTimer(2, "Calcfg");
                Hotsun.TotalTimeUsed = 0.0;

                // Set up general Manxcat Application structure
                Hotsun.FullSecondDerivative = false;

                //  Increase Run Number
                if ((ProcessingOption == 0) || (ProcessingOption == 100))
                    ++Configuration.RunNumber;

                //  Construct One Line Label
                Configuration.Pattern = SALSAUtility.ParallelPattern;
                SALSAUtility.PatternLabel =
                    String.Format(
                        "==== MDS {0} ==== Option:{1} PointCount_Global:{2} ==== Run {3} {4} == File {5} == {6} ==== ",
                        SALSAUtility.ParallelPattern, Configuration.ProcessingOption, SALSAUtility.PointCount_Global,
                        Configuration.RunSetLabel, Configuration.RunNumber, Configuration.DistanceMatrixFile,
                        SALSAUtility.startTime.ToLocalTime());
                SALSAUtility.SALSAPrint(0, SALSAUtility.PatternLabel);


                string startstring = "";

                if (Configuration.Comment != "")
                    startstring = "\n";
                Configuration.Comment += startstring + SALSAUtility.PatternLabel;

                // Initial Processing Complete
                SALSAUtility.MPI_communicator.Barrier();
                    // Make certain all processes have processed original data before writing updated

                // Set results directory

                // Begin Changes smbeason 8-21-2009
                // This change is to support a "Standard" naming schema for all our results.  It support the performance runs.
                string timestamp = DateTime.Now.ToString("'D['yyyy'-'MM'-'dd'] T['HH'-'mm'-'ss']'");
                string pattern = string.Format("{0}x{1}x{2}", SALSAUtility.ThreadCount, SALSAUtility.MPIperNodeCount,
                                               SALSAUtility.NodeCount);

                if (ProcessingOption != 12)
                {
                    Configuration.ResultDirectoryExtension = string.Format("ManxcatMDS {0} {1} P[{2}]",
                                                                           Configuration.RunSetLabel, timestamp, pattern);
                    ResultDirectoryName = Configuration.BaseResultDirectoryName + "\\" +
                                          Configuration.ResultDirectoryExtension;
                }

                // ManxcatCentral.MetadataforRun.ResultDirectoryExtension = "Results-" + ManxcatCentral.MetadataforRun.RunSetLabel + "-R" + ManxcatCentral.MetadataforRun.RunNumber.ToString();
                // ManxcatCentral.ResultDirectoryName = ManxcatCentral.MetadataforRun.BaseResultDirectoryName + "\\" + ManxcatCentral.MetadataforRun.ResultDirectoryExtension;
                // End Changes smbeason 8-21-2009

                if (SALSAUtility.MPIIOStrategy > 0)
                    ResultDirectoryName += "-Unit" + SALSAUtility.MPI_Size.ToString();

                if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0))
                {
                    if (Directory.Exists(ResultDirectoryName))
                    {
                        SALSAUtility.SALSAPrint(0, "The directory " + ResultDirectoryName + " exists");
                    }
                    else
                    {
                        DirectoryInfo di = Directory.CreateDirectory(ResultDirectoryName);
                        SALSAUtility.SALSAPrint(0,
                                                "The directory " + ResultDirectoryName + " was created successfully at "
                                                + Directory.GetCreationTime(ResultDirectoryName));
                    }
                }

                if (ProcessingOption == 12)
                {
                    ManxcatMDSDataProcessing.UpdateManxcatMDS_Option12();

                    if ((SALSAUtility.MPIIOStrategy > 1) || (SALSAUtility.MPI_Rank == 0))
                    {
                        string[] split = ConfigurationFileName.Split(new[] {'\\', '/'});
                        string ControlFileLastName = split[split.Length - 1];
                        SALSAUtility.SALSAUpdateMetaData(Configuration);
                        // ConfMgr.SaveAs(Configuration.ControlDirectoryName + "\\" + ControlFileLastName, true);
                        // MetaDataIO.WriteControl_Cluster(ManxcatCentral.Configuration.ControlDirectoryName + "\\" + ControlFileLastName, ref ManxcatCentral.Configuration); // Updated Control data
                    }
                    Console.WriteLine("Completed");
                    return;
                }

                // Begin Changes smbeason 8-21-2009
                // for our performance runs, the control file lives on the headnode, whereas the data file lives on the local compute node.  Thus, we need to specify
                // the full path to the data file in the control file rather than building it here.
                ActualDataFileName = Configuration.DistanceMatrixFile;
                // ActualDataFileName = ManxcatCentral.MetadataforRun.ControlDirectoryName + "\\" + ManxcatCentral.MetadataforRun.DataFileName;
                // End Changes smbeason 8-21-2009

                SALSAUtility.MPI_communicator.Barrier();
                    // Make certain all processes have processed original data before writing updated

                //  Write out Current Summary and Updated Control File
                if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0))
                {
                    ConfMgr.SaveAs(Configuration.SummaryOutputFileName);
                }

                if ((SALSAUtility.MPIIOStrategy > 1) || (SALSAUtility.MPI_Rank == 0))
                {
                    string[] split = ConfigurationFileName.Split(new[] {'\\', '/'});
                    string ControlFileLastName = split[split.Length - 1];
                    // ConfMgr.SaveAs(Configuration.ControlDirectoryName + "\\" + ControlFileLastName, true);                    
                }


                // Actually run application
                SALSAUtility.MPI_communicator.Barrier();
                    // Make certain all processes have processed original data before writing updated
                if (ProcessingOption <= 90)
                {
                    ManxcatControl(ManxcatMDSRoutines.Calcfg, ManxcatMDSRoutines.SetupMDSasChisq,
                                   ManxcatMDSRoutines.InitializeParameters,
                                   ManxcatMDSRoutines.SolveMatrix, ManxcatMDSWriteSolution.WriteMDS,
                                   ManxcatMDSRoutines.FindQlimits,
                                   ManxcatMDSRoutines.GlobalMatrixVectorProduct,
                                   ManxcatMDSRoutines.Sequel);
                }
                else if (ProcessingOption <= 100)
                {
                    /* Special case of <= 100 where actualy ManxcatAsChisq is ommitted. 
                     * Original distances are read in from file along with previously
                     * calculated coordinates specified by initialization file and density
                     * graphs along with html page are created */
                    ManxcatMDSRoutines.SetupMDSasChisq();
                    // fill hotsun.globalparam with initialization file
                    ManxcatMDSRoutines.FillupHotsun();
                    // Todo (html+density) - something is different when using built in manxcat population
//                    ManxcatMDSRoutines.InitializeParameters(Hotsun.CurrentSolution, Hotsun.InitializationLoops);
                    ManxcatMDSRoutines.Sequel();
                }
                else if (ProcessingOption <= 200)
                {
                    ManxcatControl(ManxcatRotateMDSRoutines.Calcfg,
                                   ManxcatRotateMDSRoutines.SetupRotateMDS,
                                   ManxcatRotateMDSRoutines.InitializeParameters,
                                   GenericManxcat.SolveMatrix, ManxcatRotateMDSRoutines.WriteRotateMDS,
                                   GenericManxcat.FindQlimits, GenericManxcat.GlobalMatrixVectorProduct,
                                   ManxcatRotateMDSRoutines.Sequel);
                }
                SALSAUtility.MPI_communicator.Barrier();

                string TimingMessage = "\nTiming ";
                for (int subtimer = 0; subtimer < SALSAUtility.NumberofSubTimings; subtimer++)
                {
                    if (!SALSAUtility.SubTimingEnable[subtimer])
                        continue;
                    double SubTime = SALSAUtility.SubDurations[subtimer];
                    TimingMessage += " " + SALSAUtility.SubTimingNames[subtimer] + " " + SubTime.ToString("F0") + " (" +
                                     (SubTime/SALSAUtility.HPDuration).ToString("F3") + ")";
                }
                SALSAUtility.SALSAPrint(0, TimingMessage);

                if ((SALSAUtility.MPIIOStrategy > 0) || (SALSAUtility.MPI_Rank == 0))
                {
                    MetaDataIO.WriteResults_Cluster(Configuration.SummaryOutputFileName, SALSAUtility.CosmicOutput);
                    SALSAUtility.WriteTiming_Cluster(Configuration.TimingOutputFileName, Configuration.RunSetLabel,
                                                     Configuration.RunNumber, Configuration.DistanceMatrixFile,
                                                     Environment.ProcessorName);
                    ManxcatHtmlUtility.WriteHTML();
                }
            }
            finally
            {
                // Begin Changes smbeason 8-21-2009
                SALSAParallelism.TearDownParallelism();
                // End Changes smbeason 8-21-2009
            }
        }

        private static void ReadConfiguration()
        {
            // load configuration from file
            ConfMgr = ConfigurationMgr.LoadConfiguration(ConfigurationFileName, true);
            Configuration = ConfMgr.ManxcatSection;

            SALSAUtility.NumberOriginalPoints = Configuration.DataPoints;
            SALSAUtility.CalcFixedCrossFixed = Configuration.CalcFixedCrossFixed;
            SALSAUtility.StoredDistanceOption = Configuration.StoredDistanceOption;
            SALSAUtility.DiskDistanceOption = Configuration.DiskDistanceOption;
            SALSAUtility.UndefinedDistanceValue = Configuration.UndefinedDistanceValue;

            SALSAUtility.DistanceCut = Configuration.DistanceCut;
            SALSAUtility.LinkCut = Configuration.LinkCut;
            SALSAUtility.TransformMethod = Configuration.TransformMethod;
            SALSAUtility.TransformParameter = Configuration.TransformParameter;

            SALSAUtility.NodeCount = Configuration.NodeCount; // Number of Nodes

            SALSAUtility.Xmaxbound = Configuration.XmaxBound > 0 ? Configuration.XmaxBound : 1.5;
            SALSAUtility.Ymaxbound = Configuration.YmaxBound > 0 ? Configuration.YmaxBound : 1.5;
            SALSAUtility.Xres = Configuration.Xres > 0 ? Configuration.Xres : 50;
            SALSAUtility.Yres = Configuration.Yres > 0 ? Configuration.Yres : 50;
            SALSAUtility.Alpha = Configuration.Alpha > 0 ? Configuration.Alpha : 2;
            SALSAUtility.Pcutf = Configuration.Pcutf > 0 ? Configuration.Pcutf : 0.85;
            SALSAUtility.Normalize = Configuration.Normalize;
            SALSAUtility.ClusterFile = Configuration.ClusterFile;
            SALSAUtility.SelectedClusters = new HashSet<int>();
            if (!string.IsNullOrEmpty(Configuration.SelectedClusters) && !"none".Equals(Configuration.SelectedClusters))
            {
                var sep = new[] {','};
                string[] splits = Configuration.SelectedClusters.Trim().Split(sep);
                foreach (string split in splits)
                {
                    int cnum = int.Parse(split);
                    if (!SALSAUtility.SelectedClusters.Contains(cnum))
                    {
                        SALSAUtility.SelectedClusters.Add(cnum);
                    }
                }
            }

            if (SALSAUtility.SelectedClusters.Count > 0 && !File.Exists(SALSAUtility.ClusterFile))
            {
                throw new Exception("Cluster file not specified to decide selected clusters");
            }

            SALSAUtility.IsClustersSelected = SALSAUtility.SelectedClusters.Count > 0;

            SALSAUtility.ManxcatRunName = !string.IsNullOrEmpty(Configuration.ManxcatRunName)
                                              ? Configuration.ManxcatRunName
                                              : "Unspecified Run";
            SALSAUtility.ManxcatRunDescription = !string.IsNullOrEmpty(Configuration.ManxcatRunDescription)
                                                     ? Configuration.ManxcatRunDescription
                                                     : "Description not specified";
        }

        public static void ManxcatControl(Hotsun.CalcfgSignature Calcfg, Hotsun.IntializeSignature InitializeApplication,
                                          Hotsun.InitializeParametersSignature InitializeParameters,
                                          Hotsun.SolveMatrixSignature SolveMatrix,
                                          Hotsun.WriteSolutionSignature WriteSolution,
                                          Hotsun.FindQlimitsSignature FindQlimits,
                                          Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct,
                                          Hotsun.SequelSignature Sequel)
        {
            //  Set up Specific Manxcat Applicatin
            InitializeApplication(); // Reads in Data

            while (Hotsun.InitializationLoopCount < Hotsun.InitializationLoops)
            {
                // Iterate over Initializations
                // First call to user routine
                Hotsun.numit = 0;
                InitializeParameters(Hotsun.CurrentSolution, Hotsun.InitializationLoopCount);
                ZeroSolution(Hotsun.CurrentSolution);
                MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
                GlobalParameterSet = true; // Set Indicator  that Global Parameters are set

                SALSAUtility.StartSubTimer(2);
                violat = Calcfg(Hotsun.CurrentSolution);
                SALSAUtility.StopSubTimer(2);

                Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
                Hotsun.tsolve = 0.0;
                Hotsun.teigen = 0.0;
                Hotsun.NumberofCGIterations = -2;
                Hotsun.EigenvalueIndicator1 = 0;
                Hotsun.EigenvalueIndicator2 = 0;

                launchQlimits(FindQlimits); // Find limits on Q

                Hotsun.ichsav = -1; // changed from = 1 in old Fortran which seems wrong
                Hotsun.chsave[0] = Hotsun.zerocr; // Unnecessary
                Hotsun.expchg = 0.0;
                Hotsun.delchi = 0.0;
                Hotsun.isaved = 0;
                Hotsun.bnderr = 0;
                Hotsun.igood = 0;
                Hotsun.CountToSD = 0;
                Hotsun.CountQSmall = 0;
                Hotsun.CountQLarge = 0;
                Hotsun.CountToRandom = 0;
                Hotsun.CountExercusion = 0;
                Hotsun.QgoodFactor = 1.0;

                Hotsun.Q = Hotsun.QHighInitialFactor*Hotsun.Qhigh;
                bool QlessthanQlow = Hotsun.Q < Hotsun.Qlow;
                SALSAUtility.SynchronizeMPIvariable(ref QlessthanQlow);

                if (QlessthanQlow)
                    Hotsun.Q = Hotsun.Qlow;
                Hotsun.Qbest = Hotsun.Q;

                Hotsun.Qgood = -1.0;
                Hotsun.isdtry = 0;

                if (Hotsun.InitialSteepestDescents > 0)
                    Hotsun.isdtry = 1;

                // Print Initial Trial
                SALSAUtility.SALSAStatus(ResultDirectoryName, " Initial Chisq "
                                                              +
                                                              (ChisqPrintConstant*Hotsun.zerocr).ToString
                                                                  ("F3"));

                if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
                    SALSAUtility.SALSAPrint(1,
                                            "\n-----------------------------------------------\n" +
                                            SALSAUtility.startTime.ToLocalTime() + " Iterations "
                                            + Hotsun.numit.ToString() + " Chisq "
                                            + (ChisqPrintConstant*Hotsun.zerocr).ToString("F3") + " Parameters " +
                                            Hotsun.npar.ToString() + " Data Points " + Hotsun.ndata.ToString());

                Hotsun.DecisionMethod = 0; //   Initial Conditions on Decisions
                Hotsun.DecisionMethod_1 = 0;
                Hotsun.DecisionMethod_2 = 0;
                Hotsun.DecisionMethod_3 = 0;
                Hotsun.DecisionMethod_4 = 0;
                Hotsun.DecisionMethod_5 = 0;
                Hotsun.DecisionMethod_LineSearch = 0;

                //  Save initial guess as best solution
                SaveBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
                BestSolutionSet = true;
                BestxshiftSet = false;

                int CoordinateWriteCount = 0;
                while (true)
                {
                    // Iterate Chisq Solution
                    //  Save diagonal elements removed as we no longer overwrite by adding Matquardt's good idea Q

                    // Hotsun.PreviousSolution.param holds parameters of Starting Solution
                    // while Hotsun.PreviousSolution.xshift was shift that was added to make this param
                    CopySolution(Hotsun.PreviousSolution, Hotsun.CurrentSolution);

                    //  Start code used each iteration whether or not last iteration succeeded
                    Hotsun.materr = 0;

                    bool wefailed = false;
                    Hotsun.DecisionMethod_LineSearch = 0;

                    while (true)
                    {
                        // Solve Matrix and Find Good Q value if fails

                        //  Process Options for Q getting Stuck
                        if (Hotsun.CountQSmall < 0)
                            Hotsun.CountQSmall = 0;

                        if (Hotsun.CountQLarge < 0)
                            Hotsun.CountQLarge = 0;
                        double Qtest = 0.05*Hotsun.Qhigh;
                        Qtest = Math.Max(2.0*Hotsun.Qlow, Qtest);
                        bool Qisaboveaverage = Hotsun.Q > Qtest;
                        SALSAUtility.SynchronizeMPIvariable(ref Qisaboveaverage);

                        if (Qisaboveaverage)
                        {
                            Hotsun.CountQSmall = 0;
                            ++Hotsun.CountQLarge;
                        }
                        else
                        {
                            Hotsun.CountQLarge = 0;
                            ++Hotsun.CountQSmall;
                        }

                        if (Hotsun.CountQSmall > Hotsun.ResetQLimit)
                        {
                            Hotsun.CountQSmall = -1;
                            Hotsun.CountQLarge = 0;
                            Hotsun.Q = Qtest*4.0;
                            Hotsun.DecisionMethod_5 = 1;
                        }

                        if (Hotsun.CountQLarge > Hotsun.ResetQLimit)
                        {
                            Hotsun.CountQLarge = -1;
                            Hotsun.CountQSmall = 0;
                            Hotsun.Q = Math.Max(Hotsun.Qlow, 0.1*Qtest);
                            Hotsun.DecisionMethod_5 = 2;
                        }

                        // Process Steepest Descent Options
                        if ((Hotsun.isdtry == 2) || (Hotsun.isdtry == 4))
                            Hotsun.isdtry = 0;

                        if ((Hotsun.isdtry == 1) || (Hotsun.isdtry == 3))
                        {
                            // Steepest Descent Solution  -- Q is NOT added
                            Hotsun.AddMarquardtQDynamically = false;
                            SteepestDescentSolution(Hotsun.CurrentSolution, GlobalMatrixVectorProduct);
                            ++Hotsun.isdtry;
                            Hotsun.Qsd = Hotsun.Q;
                            Hotsun.materr = 0;
                            break;
                        }
                        else
                        {
                            // Typical Full solution

                            Hotsun.UseDiagonalScaling = Hotsun.UseDiagonalScalinginSolvers;

                            Hotsun.FullSecondDerivative = false;

                            if (Configuration.FullSecondDerivativeOption == 1)
                                Hotsun.FullSecondDerivative = true;

                            if (Hotsun.UseDiagonalScaling)
                                SetupDiagonalScaling(Hotsun.CurrentSolution); // This does NOT change matrices

                            // Add in Marquardt Parameter (removed afterwards)
                            Hotsun.AddMarquardtQDynamically = !Hotsun.AddMarquardtQExplicitly;
                            AddSubtractMarquardtFactor(Hotsun.CurrentSolution, 1.0);

                            //  Invert Modified Chisq Matrix
                            SALSAUtility.StartSubTimer(0);
                            bool SolveMatrixSuccess = SolveMatrix(Hotsun.CurrentSolution.xshift, Hotsun.CurrentSolution);
                            Hotsun.tsolve = SALSAUtility.StopSubTimer(0);
                            AddSubtractMarquardtFactor(Hotsun.CurrentSolution, -1.0);

                            if (SolveMatrixSuccess)
                            {
                                Hotsun.materr = 0;
                                break;
                            }
                            else
                            {
                                // error in Matrix Solver -- either give up or throw away trial and increase Q
                                if (Hotsun.materr > 1)
                                {
                                    // Too many failures
                                    EndupManxcat(6, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                                    wefailed = true;
                                    break;
                                }

                                ++Hotsun.materr;

                                // First try in matrix failure is increasing Q
                                Hotsun.Q = 4.0*Hotsun.Q;

                                bool Qessentiallyzero = Hotsun.Q <= (2.0*Hotsun.QHighInitialFactor);
                                SALSAUtility.SynchronizeMPIvariable(ref Qessentiallyzero);
                                if (Qessentiallyzero)
                                    Hotsun.Q = 1.0;

                                // Second try in matrix failure is steepest descent
                                if (Hotsun.materr == 2)
                                    Hotsun.isdtry = 3;
                                continue;
                            }
                        }
                    } // End while Solving Matrix and finding good Q value if failure

                    if (wefailed)
                        break;

                    // We only reach here if Matrix Solved Correctly by Steepest Descent or Marquardt
                    // Old Manxcat used two step solution (LU and Solve)
                    //  Here we just use direct Conjugate Gradient

                    //  Calculate Predicted Chisq Change from this solution
                    FindPredictedChisqChange(Hotsun.CurrentSolution, GlobalMatrixVectorProduct);
                    Hotsun.expchg = Hotsun.pred1 + Hotsun.pred2;

                    if (Hotsun.pred2 > 0)
                    {
                        if (SALSAUtility.MPI_Rank == 0)
                            SALSAUtility.SALSAPrint(2,
                                                    "\nIllegal pred2 " +
                                                    (Hotsun.pred2*ChisqPrintConstant).ToString("F3") + " CurrentChisq "
                                                    +
                                                    (Hotsun.CurrentSolution.Chisquared*ChisqPrintConstant).ToString("F3"));
                    }

                    // For derivative test avoid first iteration that could be special
                    if (Hotsun.derivtest && (Hotsun.numit == 3))
                    {
                        DoDerivTest(Calcfg, GlobalMatrixVectorProduct);
                        Hotsun.derivtest = false;
                    }
                    Hotsun.doingderivtest = false;

                    //  Set up parameters for next call to Calcfg
                    //  This adds in estimated change to param and zeros first derivative, Second Derivative and Chisq (zerocr) in a  CurrentSolution
                    // xshift was actually calculated from PREVIOUS Current Solution (now stored in Previous Solution)
                    SALSABLAS.LinearCombineVector(Hotsun.CurrentSolution.param, 1.0, Hotsun.CurrentSolution.param, -1.0,
                                                  Hotsun.CurrentSolution.xshift);
                    MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
                    GlobalParameterSet = true;
                    ZeroSolution(Hotsun.CurrentSolution); // Does not vary param or xshift

                    // Save Current best Chisq in circular list
                    ++Hotsun.ichsav;

                    if (Hotsun.ichsav >= Hotsun.nbadgo)
                        Hotsun.ichsav = 0;
                    Hotsun.chsave[Hotsun.ichsav] = Hotsun.zeromn;

                    Hotsun.numit++; // Increment Iteration count
                    Hotsun.idata = 0;


                    //  Call Calcfg to calculate Taylor expansion
                    SALSAUtility.StartSubTimer(2);
                    violat = Calcfg(Hotsun.CurrentSolution);
                    SALSAUtility.StopSubTimer(2);

                    //  Cope with boundary errors
                    // In this case parameter values will have been changed by Calcfg
                    if (violat)
                    {
                        Hotsun.bnderr++;
                        Findxnorm(2); // Reset xshift and xnorm to reflect Calcfg resetting
                    }
                    else
                    {
                        Hotsun.bnderr = 0;
                    }

                    //   set change in chisq
                    Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;

                    //  Look at case where actual change significantly larger than expected
                    // Only do search if last solution was also a success
                    bool LineSearchUsed = false;
                    bool randomized = false;
                    bool searchsuccess = false;
                    double LineFactor = 1.0;
                    double ExtraDecrease = 0.0;
                    int LineIterations = 0;
                    bool Delchipositive = Hotsun.delchi >= 0.0;
                    SALSAUtility.SynchronizeMPIvariable(ref Delchipositive);

                    if (Delchipositive)
                    {
                        if ((Hotsun.OmegaOption > 0) && (Hotsun.igood > 0))
                        {
                            bool LookatOmegaOption = Hotsun.delchi > (Hotsun.Omega*Hotsun.expchg);
                            SALSAUtility.SynchronizeMPIvariable(ref LookatOmegaOption);

                            if (LookatOmegaOption)
                            {
                                // Do a line search exploring larger shifts
                                LineSearchUsed = true;
                                LineFactorGuess = -1.0; // Explore larger solutions
                                searchsuccess = LineSearch(Hotsun.CurrentSolution, Hotsun.PreviousSolution,
                                                           out Hotsun.DecisionMethod_LineSearch, out LineFactor,
                                                           out ExtraDecrease, out LineIterations, Calcfg, -1.0);
                                Hotsun.DecisionMethod_LineSearch += 100;
                                Hotsun.TotalSearchIterations += LineIterations;
                            }
                        }
                    }
                    else
                    {
                        if ((Hotsun.pred1 > 0.0) && (Hotsun.pred2 < 0.0))
                        {
                            double dchi1 = Hotsun.pred1*ChisqPrintConstant;

                            if (dchi1 < 0.5*Hotsun.zerocr)
                            {
                                if (Hotsun.CountExercusion == 0)
                                    ++Hotsun.CountToRandom;

                                if ((Hotsun.CountToRandom >= Hotsun.RandomLimit) && ((Hotsun.maxit - Hotsun.numit) > 20))
                                {
                                    Hotsun.CountToRandom = 0;
                                    randomized = true;
                                }
                                else
                                {
                                    double dchi4 = Hotsun.delchi*ChisqPrintConstant;
                                    double alpha = dchi1/(2.0*(dchi1 - dchi4));
                                    LineFactorGuess = alpha;
                                    searchsuccess = LineSearch(Hotsun.CurrentSolution, Hotsun.PreviousSolution,
                                                               out Hotsun.DecisionMethod_LineSearch, out LineFactor,
                                                               out ExtraDecrease, out LineIterations, Calcfg, alpha);
                                    Hotsun.TotalSearchIterations += LineIterations;
                                    Hotsun.zerocr = Hotsun.CurrentSolution.Chisquared;
                                    Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;
                                    Delchipositive = Hotsun.delchi >= 0.0;
                                    LineSearchUsed = true;
                                }
                            }
                        }
                    }

                    // Change Information
                    double size1 = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.xshift,
                                                                 Hotsun.CurrentSolution.xshift);
                    double size2 = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first,
                                                                 Hotsun.CurrentSolution.first);
                    double Changexshiftfirst = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first,
                                                                             Hotsun.CurrentSolution.xshift);
                    Changexshiftfirst = Changexshiftfirst/Math.Sqrt(size1*size2);
                    double Changeinfirst = 0.0;
                    double Changeinxshift = 0.0;
                    double xshiftRatio = 0.0;
                    double firstRatio = 0.0;

                    // Changeinfirst is scalar product of first derivative at I and First Derivative at II / Size ( I * II)
                    // firstratio is Ratio( Size first for II / Size first for I)
                    // In most algorithms I is Best solution and II current solution
                    // If exploring away from best I is previous solution
                    double size3 = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.first,
                                                                 Hotsun.PreviousSolution.first);
                    Changeinfirst = SALSABLAS.VectorScalarProduct(Hotsun.CurrentSolution.first,
                                                                  Hotsun.PreviousSolution.first);
                    Changeinfirst = Changeinfirst/Math.Sqrt(size2*size3);
                    firstRatio = Math.Sqrt(size2/size3);

                    // Changeinxshift is scalar product of xshift at I and xshift at II / Size ( I * II)
                    // xshiftratio is Ratio( Size xshift for II / Size xshift for I)
                    // In most algorithms I is Best solution and II current solution
                    // If exploring away from best I is previous solution
                    bool shiftchangeset = false;
                    if (Hotsun.numit > 1)
                    {
                        shiftchangeset = true;
                        double size4 = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.xshift,
                                                                     Hotsun.PreviousSolution.xshift);
                        Changeinxshift = SALSABLAS.VectorScalarProduct(Hotsun.PreviousSolution.xshift,
                                                                       Hotsun.CurrentSolution.xshift);
                        Changeinxshift = Changeinxshift/Math.Sqrt(size1*size4);
                        xshiftRatio = Math.Sqrt(size1/size4);
                    }

                    //   Reset change in chisq
                    Hotsun.delchi = Hotsun.PreviousSolution.Chisquared - Hotsun.zerocr;

                    //  Timing Information
                    // Hotsun.tsolve and Hotsun.teigen set in StopTimer
                    Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
                    SALSAUtility.InterimTiming();
                    Hotsun.TotalTimeUsed = SALSAUtility.HPDuration;

                    //  Print Summary of this solution
                    SALSAUtility.SALSAStatus(ResultDirectoryName,
                                             "In Loop " + Hotsun.InitializationLoopCount.ToString() + " Iteration "
                                             + Hotsun.numit.ToString() + " Chisq " +
                                             (ChisqPrintConstant*Hotsun.zerocr).ToString("F3") + " Best Chisq " +
                                             (ChisqPrintConstant*Hotsun.BestSolution.Chisquared).ToString("F3"));

                    if ((SALSAUtility.DebugPrintOption > 1) && (SALSAUtility.MPI_Rank == 0))
                    {
                        string shift1 = "Unset";
                        string shift2 = "unset";
                        if (shiftchangeset)
                        {
                            shift1 = Changeinxshift.ToString("F4");
                            shift2 = xshiftRatio.ToString("F4");
                        }
                        SALSAUtility.SALSAPrint(2, "\nLoop "
                                                   + Hotsun.InitializationLoopCount.ToString() + " Iteration "
                                                   + Hotsun.numit.ToString() + " Q "
                                                   + Hotsun.Q.ToString("E4") + " Qlow "
                                                   + Hotsun.Qlow.ToString("F3") + " ("
                                                   + Hotsun.EigenvalueIndicator2.ToString() + ") Qhigh "
                                                   + Hotsun.Qhigh.ToString("F3") + " ("
                                                   + Hotsun.EigenvalueIndicator1.ToString() + ") Qgood "
                                                   + Hotsun.Qgood.ToString("E4") + " Trace Q "
                                                   + Hotsun.ChisqMatrixTrace.ToString("F3") + " Norm Q " +
                                                   Hotsun.ChisqMatrixNorm.ToString("F3") +
                                                   "\nDecision " + Hotsun.DecisionMethod.ToString() + " 1:" +
                                                   Hotsun.DecisionMethod_1.ToString() + " 2:" +
                                                   Hotsun.DecisionMethod_2.ToString()
                                                   + " 3:" + Hotsun.DecisionMethod_3.ToString() + " 4:" +
                                                   Hotsun.DecisionMethod_4.ToString() + " 5:" +
                                                   Hotsun.DecisionMethod_5.ToString() + " Counts Q Low "
                                                   + Hotsun.CountQSmall.ToString() + " Q High "
                                                   + Hotsun.CountQLarge.ToString() + " SD Count "
                                                   + Hotsun.CountToSD.ToString() + " Randoms "
                                                   + Hotsun.CountToRandom.ToString() + " Exercusions " +
                                                   Hotsun.CountExercusion.ToString()
                                                   + " SD Opt " + Hotsun.isdtry.ToString() + " Solver Errs "
                                                   + Hotsun.materr.ToString() + " CG Iters "
                                                   + Hotsun.NumberofCGIterations.ToString()
                                                   + "\nScalProd Shift.First " + Changexshiftfirst.ToString("F4") +
                                                   " Change Shift from Prev DotProd "
                                                   + shift1 + " and Change Shift in norm " + shift2 +
                                                   "\nChange in first deriv from Prev DotProd "
                                                   + Changeinfirst.ToString("F4") + " and in first deriv norm " +
                                                   firstRatio.ToString("F4"));

                        string chisqmethod = "Marquardt's Algorithm";

                        if (Hotsun.FullSecondDerivative)
                            chisqmethod = "Corrected Marquardt";

                        if ((Hotsun.isdtry == 2) || (Hotsun.isdtry == 4))
                        {
                            chisqmethod = "Steepest Descent";
                            if (Hotsun.isdtry == 4)
                                chisqmethod += " from Singularity";
                            if (Hotsun.FullSecondDerivative)
                                chisqmethod = "Corrected Steep Desc";
                        }

                        string chisqstatus = "OK Parameters";

                        if (violat)
                            chisqstatus = "Out of Range Parameter";
                        double dchi1 = Hotsun.pred1*ChisqPrintConstant;
                        double dchi2 = Hotsun.pred2*ChisqPrintConstant;
                        double dchi3 = Hotsun.pred3*ChisqPrintConstant;
                        double dchi4 = Hotsun.delchi*ChisqPrintConstant;

                        if (LineSearchUsed)
                            SALSAUtility.SALSAPrint(2, "Linesearch Chisq gains "
                                                       + (ChisqPrintConstant*ExtraDecrease).ToString("F3") +
                                                       " Shift Factor "
                                                       + LineFactor.ToString("F3") + " with guess " +
                                                       LineFactorGuess.ToString("F3") + " Iterations " +
                                                       LineIterations.ToString() + " Method " +
                                                       Hotsun.DecisionMethod_LineSearch.ToString());

                        // Output parameter values if not too many
                        string parametervalues = "";

                        if ((!Hotsun.DecomposeParameters) && (Hotsun.npar < 20))
                        {
                            parametervalues = "\nParameters:";

                            for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                            {
                                for (int LocalVectorIndex = 0;
                                     LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                     LocalVectorIndex++)
                                {
                                    double tmp = Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
                                    parametervalues += " " + tmp.ToString("E4");
                                }
                            }
                        }

                        string labelCG = "Not Scaled ";
                        if (Hotsun.NumberofCGIterations > 0)
                        {
                            labelCG = "Scaled ";
                            Hotsun.tsolve = Hotsun.tsolve/Hotsun.NumberofCGIterations;
                        }
                        string labeleigen = "Not Scaled ";
                        if ((Hotsun.EigenvalueIndicator1 > 0) && (Hotsun.EigenvalueIndicator2 > 0))
                        {
                            labeleigen = "Scaled ";
                            Hotsun.teigen = Hotsun.teigen/(Hotsun.EigenvalueIndicator1 + Hotsun.EigenvalueIndicator2);
                        }

                        SALSAUtility.SALSAPrint(2, "Chisq "
                                                   + (ChisqPrintConstant*Hotsun.zerocr).ToString("F3") + " Best Chisq " +
                                                   (ChisqPrintConstant*Hotsun.BestSolution.Chisquared).ToString("F3") +
                                                   " " + chisqmethod + " " + chisqstatus + " Calcfg Time " +
                                                   Hotsun.tcalcfg.ToString("F0") +
                                                   " Matrix Solve " + labelCG + Hotsun.tsolve.ToString("F1") +
                                                   " Eigen " + labeleigen + Hotsun.teigen.ToString("F1") +
                                                   " Total " + Hotsun.TotalTimeUsed.ToString("F0") +
                                                   "\n1st Deriv Reduction " + dchi1.ToString("F4") +
                                                   " Exp 2nd " + dchi2.ToString("F4") + " Other 2nd " +
                                                   dchi3.ToString("F4") +
                                                   " Actual Reduction " + dchi4.ToString("F4") +
                                                   parametervalues);

                        Hotsun.tsolve = 0.0;
                        Hotsun.teigen = 0.0;
                        Hotsun.tcalcfg = 0.0;
                        Hotsun.NumberofCGIterations = -2;
                        Hotsun.EigenvalueIndicator1 = 0;
                        Hotsun.EigenvalueIndicator2 = 0;
                    } // End Print Out

                    // Save Results periodically
                    CoordinateWriteCount++;
                    if (CoordinateWriteCount == Hotsun.CoordinateWriteFrequency)
                    {
                        CoordinateWriteCount = 0;
                        CopySolution(Hotsun.SearchSolution1, Hotsun.CurrentSolution);
                        CopySolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
                        double chisave = Hotsun.zerocr;
                        Hotsun.zerocr = Hotsun.zeromn;

                        SALSABLAS.zrword(Hotsun.perr);
                        MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
                        GlobalParameterSet = true; // Set Indicator  that Global Parameters are set

                        WriteCoordinates(0, WriteSolution);

                        CopySolution(Hotsun.CurrentSolution, Hotsun.SearchSolution1);
                        Hotsun.zerocr = chisave;
                        MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
                    }

                    //  Examine Nature of Solution
                    Hotsun.DecisionMethod = 1; // Default non Initialization
                    Hotsun.DecisionMethod_1 = 0;
                    Hotsun.DecisionMethod_2 = 0;
                    Hotsun.DecisionMethod_3 = 0;
                    Hotsun.DecisionMethod_4 = 0;
                    Hotsun.DecisionMethod_5 = 0;

                    if (Hotsun.CountExercusion > 0)
                    {
                        ++Hotsun.CountExercusion;
                        Hotsun.DecisionMethod_1 = 1;

                        if (Hotsun.CountExercusion >= Hotsun.ExercusionLimit)
                        {
                            Hotsun.DecisionMethod_1 = 2;
                            bool CurrentSolutionBest = Hotsun.CurrentSolution.Chisquared <
                                                       Hotsun.BestSolution.Chisquared;
                            SALSAUtility.SynchronizeMPIvariable(ref CurrentSolutionBest);

                            if (!CurrentSolutionBest)
                            {
                                RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
                                Hotsun.isaved = 2;
                            }
                            Hotsun.CountExercusion = 0;
                            Hotsun.CountToRandom = 0;
                        }
                    }

                    if (Delchipositive)
                    {
                        // Chisq Decreased - Success

                        // If necessary, save current best parameters
                        // We save true (diagonal of) second derivative matrix -- not one overwritten by Marquardt
                        bool CurrentSolutionBest = Hotsun.CurrentSolution.Chisquared < Hotsun.BestSolution.Chisquared;
                        SALSAUtility.SynchronizeMPIvariable(ref CurrentSolutionBest);

                        if (CurrentSolutionBest)
                        {
                            SaveBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
                            BestSolutionSet = true;
                            BestxshiftSet = true;
                            Hotsun.CountExercusion = 0;
                        }

                        Hotsun.DecisionMethod = 2;
                        Hotsun.succ = true;
                        Hotsun.isaved = 0;

                        if (Hotsun.isdtry == 2)
                        {
                            Hotsun.DecisionMethod = 3;
                        }

                        if (Hotsun.isdtry == 4)
                        {
                            Hotsun.DecisionMethod = 4;
                            bool QsdTest = Hotsun.Q < (Hotsun.Qsd - 0.0000000001);
                            SALSAUtility.SynchronizeMPIvariable(ref QsdTest);

                            if (QsdTest)
                            {
                                Hotsun.DecisionMethod = 5;
                            }
                        }

                        if (LineSearchUsed)
                        {
                            Hotsun.DecisionMethod += 20;

                            if (!searchsuccess)
                                Hotsun.DecisionMethod += 20;
                        }
                    }
                    else
                    {
                        // Chisq Increases -- Failure
                        Hotsun.DecisionMethod = 11;
                        Hotsun.isaved = 1;

                        if (randomized)
                        {
                            Hotsun.DecisionMethod = 12;
                            Hotsun.CountExercusion = 1;
                        }
                        Hotsun.succ = false;

                        if (Hotsun.CountExercusion == 0)
                        {
                            RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
                            Hotsun.isaved = 2;
                        }
                    }

                    //  Decide if to give up ghost and call this the last Iteration
                    int readinstruction = -1;
                    SALSAUtility.SALSAGracefulend(ResultDirectoryName, out readinstruction);
                    SALSAUtility.SynchronizeMPIvariable(ref readinstruction);
                    if (readinstruction == 0)
                    {
                        Hotsun.InitializationLoops = Hotsun.InitializationLoopCount + 1;
                        EndupManxcat(7, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                        break;
                    }
                    else if (readinstruction > 0)
                    {
                        Hotsun.maxit = readinstruction;
                    }

                    bool TimeExceeded = (Hotsun.timmax > 0.0) && (Hotsun.TotalTimeUsed > Hotsun.timmax);
                    SALSAUtility.SynchronizeMPIvariable(ref TimeExceeded);

                    if (TimeExceeded)
                    {
                        // Time limit reached
                        EndupManxcat(3, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                        break;
                    }

                    if (Hotsun.numit >= Hotsun.maxit)
                    {
                        // Iteration limit reached
                        EndupManxcat(1, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                        break;
                    }

                    bool SmallExpectedChisqChange = (Hotsun.expchg*ChisqPrintConstant) <= Hotsun.dellim;
                    SALSAUtility.SynchronizeMPIvariable(ref SmallExpectedChisqChange);

                    if (SmallExpectedChisqChange)
                    {
                        // Expected change in chisq <= preset limit
                        EndupManxcat(2, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                        break;
                    }

                    if (Hotsun.bnderr > Hotsun.bnderrLimit)
                    {
                        // Boundary Value limit reached
                        EndupManxcat(5, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                        break;
                    }

                    // Stop if too little progress made in nbadggo steps
                    // Do this by comparing current best chisq with that found nbadgo iterations before
                    if (Hotsun.numit > Hotsun.nbadgo)
                    {
                        int itest = Hotsun.ichsav + 1;

                        if (itest >= Hotsun.nbadgo)
                            itest = 0;

                        bool TooLittleProgress = (Hotsun.chsave[itest] - Hotsun.zeromn) < Hotsun.dellim;
                        SALSAUtility.SynchronizeMPIvariable(ref TooLittleProgress);

                        if (TooLittleProgress)
                        {
                            EndupManxcat(4, WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                            break;
                        }
                    }

                    //  Another Iteration called for -- Decide on Strategy

                    //  Reasonable Success -- Decrease Marquardt Parameter if going very well
                    //  Namely actual change is larger that Fletcher's rho * Expected change estimated from Taylor expansion
                    bool GoodenoughChisqChange = Hotsun.delchi >= (Hotsun.rho*Hotsun.expchg);
                    SALSAUtility.SynchronizeMPIvariable(ref GoodenoughChisqChange);

                    if (GoodenoughChisqChange)
                    {
                        launchQlimits(FindQlimits); // Reset Q limits every now and then
                        Hotsun.DecisionMethod_4 = 2;

                        Hotsun.Qgood = Hotsun.Q;
                        ++Hotsun.igood;
                        bool Qessentiallyzero = Hotsun.Q < 0.0000000001;
                        SALSAUtility.SynchronizeMPIvariable(ref Qessentiallyzero);

                        if (Qessentiallyzero)
                        {
                            Hotsun.DecisionMethod_2 = 1;
                            continue; // Next Iteration
                        }

                        //  Modest change are those "good" iterations that change Chisq by less than Feltcher's Sigma *   Expected change estimated from Taylor expansion
                        //  Good but modest chisq changes leave Q unchanged
                        bool ChisqChangemodest = Hotsun.delchi < (Hotsun.sigma*Hotsun.expchg);
                        SALSAUtility.SynchronizeMPIvariable(ref ChisqChangemodest);

                        if (ChisqChangemodest)
                        {
                            Hotsun.DecisionMethod_2 = 2;
                            ++Hotsun.CountToSD;

                            if (Hotsun.CountToSD >= Hotsun.SDLimit)
                            {
                                Hotsun.CountToSD = 0;

                                if (Hotsun.isdtry == 0)
                                    Hotsun.isdtry = 1;

                                Hotsun.DecisionMethod_2 = 3;
                            }
                            continue; // Next Iteration
                        }

                        //  This is case of really good Chisq change -- decrease Q
                        Hotsun.CountToSD = 0;
                        bool Qalreadysmall = Hotsun.Q <= (Hotsun.Qscale*Hotsun.Qlow);
                        SALSAUtility.SynchronizeMPIvariable(ref Qalreadysmall);

                        if (Qalreadysmall)
                        {
                            if (Hotsun.igood >= 3)
                            {
                                // Reduce Q below Qlow if really good convergence
                                Hotsun.QgoodFactor *= Hotsun.QgoodReductionFactor;
                                Hotsun.Q = Hotsun.Qlow*Hotsun.QgoodFactor;
                                Hotsun.DecisionMethod_2 = 4;
                            }
                            else
                            {
                                Hotsun.DecisionMethod_2 = 5;
                            }
                            continue; // Next Iteration
                        }
                        else
                        {
                            Hotsun.DecisionMethod_2 = 6;
                            Hotsun.Q = Math.Max(Math.Sqrt(Hotsun.Q*Hotsun.Qlow), Hotsun.Q*0.33);
                            continue; // Next Iteration
                        }
                    } // End case when reasonable success -- possibly decrease Q

                        //  Unexpectedly small or negative change -- Increase Marquardt Parameter
                        // Reset limits Qlow and Qhigh on allowed range of Q values if needed
                    else
                    {
                        Hotsun.DecisionMethod_3 = 1;

                        Hotsun.igood = 0; // Restart count of good solutions
                        Hotsun.QgoodFactor = 1.0; // Restart factor for lowering Q below Qlow
                        ++Hotsun.CountToSD;

                        if (Hotsun.CountToSD >= Hotsun.SDLimit)
                        {
                            Hotsun.CountToSD = 0;

                            if (Hotsun.isdtry == 0)
                                Hotsun.isdtry = 1;

                            Hotsun.DecisionMethod_3 = 2;
                            continue; // Next Iteration
                        }
                        bool Qlookswrong = false;
                        bool qlimitsvalid = true;
                        bool QnearQhigh = false;

                        if (Hotsun.CurrentSolution.IterationCalculated != Hotsun.IterationforQlimit)
                            qlimitsvalid = false;
                        bool totaltest = qlimitsvalid || (Hotsun.numit < IterationtorecalculateQLimits);

                        if (!totaltest)
                        {
                            // Qlow and Qhigh need to be reset
                            Hotsun.DecisionMethod_4 = 1;
                            launchQlimits(FindQlimits);
                        } // End case Qlow and Qhigh need to be reset

                        QnearQhigh = Hotsun.Qhigh <= (Hotsun.Qscale*Hotsun.Q);
                        SALSAUtility.SynchronizeMPIvariable(ref QnearQhigh);

                        if (QnearQhigh)
                        {
                            // Q is near enough Qhigh and it still didn't work very well! That's very disappointing

                            Hotsun.DecisionMethod_3 = 3;

                            if (Hotsun.isdtry != 0)
                            {
                                // Previous Round was not a full Chisq calculation
                                Hotsun.Q *= 0.5;
                                Qlookswrong = false;
                                Hotsun.DecisionMethod_3 = 4;
                                Hotsun.isdtry = 2;
                            }
                            else
                            {
                                // This will select Steepest Descent next time
                                Qlookswrong = true;
                            }
                        } // End case Q is near enough Qhigh
                        else
                        {
                            Qlookswrong = false;
                            Boolean Qgoodsmall = Hotsun.Qgood <= (Hotsun.Qscale*Hotsun.Q);
                            SALSAUtility.SynchronizeMPIvariable(ref Qgoodsmall);

                            if (Qgoodsmall)
                            {
                                Hotsun.Q = Math.Sqrt(Hotsun.Qlow*Hotsun.Qhigh);
                                Hotsun.DecisionMethod_3 = 5;
                            }
                            else
                            {
                                Hotsun.Q = Math.Sqrt(Hotsun.Q*Hotsun.Qgood);
                                Hotsun.DecisionMethod_3 = 6;
                            }
                        }

                        SALSAUtility.SynchronizeMPIvariable(ref Qlookswrong);

                        if (Qlookswrong)
                        {
                            Hotsun.DecisionMethod_3 += 10;

                            if (Hotsun.isdtry == 0)
                                Hotsun.isdtry = 1;

                            Hotsun.Q = 0.5*Hotsun.Qbest;
                        }
                    } // End Case when unexpectedly poor performance -- increase Q

                    //  Another Iteration called for -- End decision on Strategy
                } // End Iterate Chisq Solution with while(true)
            } // End while over Initialization Iteration Loops

            // Perform follow-up
            Sequel();
        }

        // End Control

        // Explain Decision codes
        // DecisionMethod = 0  	Initial Conditions
        // DecisionMethod = 1	Default case if not initialization
        // DecisionMethod = 2	Default success -- chisq decreases
        // DecisionMethod = 3	chisq decreases: isdtry=2 (Steepest Descent used with no matrix failure) reset isdtry=0
        // DecisionMethod = 4	chisq decreases: isdtry=3 (Steepest Descent used with matrix failure) Q >=than previous steepest value
        // DecisionMethod = 5	chisq decreases: isdtry=3 (Steepest Descent used with matrix failure) Q < than previous steepest value; reset isdtry=04
        // DecisionMethod = 11	Default failure -- chisq increases
        // DecisionMethod = 12	Forced failure -- chisq increases and Exercusion started

        // DeltaDecisionMethod = 20	Success Chisq decreases and use Line Search Successfully
        // DeltaDecisionMethod = 40	Success Chisq decreases and use Line Search Unsuccessfully

        // DecisionMethod_1 = 0 No Exercusion
        // DecisionMethod_1 = 1 End Exercusion
        // DecisionMethod_1 = 2 Continue Exercusion

        // DecisionMethod_2 = 1	Good success as measured by rho and current Q essentially zero -- Final value
        // DecisionMethod_2 = 2	Good success as measured by rho and Modest success as measured by sigma -- Final Value
        // DecisionMethod_2 = 3	Good Success as measured by sigma and rho -- Q already small and Hotsun.igood >= 3
        // DecisionMethod_2 = 4	Good Success as measured by sigma and rho -- Q already small and Hotsun.igood < 3
        // DecisionMethod_2 = 5	Good Success as measured by sigma and rho -- Q not already small
        // DecisionMethod_2 = 6	Good success as measured by rho and Modest success as measured by sigma -- Try Steepest descent as happened too often

        // DecisionMethod_3 = 1 Case of Chisq increasing or very small change
        // DecisionMethod_3 = 2 Q looks too small -- Case of Chisq increasing or very small change
        // DecisionMethod_3 = 3 Q is near Qhigh Use Steepest Descent as not used last time -- Case of Chisq increasing or very small change
        // DecisionMethod_3 = 4 Q is near Qhigh Do not use Steepest Descent as used last time -- Case of Chisq increasing or very small change
        // DecisionMethod_3 = 5 Q NOT near enough Qhigh and current Qgood quite small -- Case of Chisq increasing or very small change
        // DecisionMethod_3 = 6 Q NOT near enough Qhigh and current Qgood quite large -- Case of Chisq increasing or very small change

        // DeltaDecisionMethod_3 = 10 Chose Steepest Descent - Case of Chisq increasing or very small change

        // DecisionMethod_4 = 1 Need to reset limits on Q -- Case of Chisq increasing or very small change
        // DecisionMethod-4 = 2 Need to reset limits on Q -- Every now and then

        // DecisionMethod_5 = 1 Q too small for too long so reset SET whether success or failure
        // DecisionMethod_5 = 2 Q too big for too long so reset SET whether success or failure

        // Iteration has ended -- Clean up in EndupManxcat
        //  ReasonToStop = 1 Iteration Limit
        //  ReasonToStop = 2 Expected Change in Chisq too small
        //  ReasonToStop = 3 Time exhausted
        //  ReasonToStop = 4 Progress too
        //  ReasonToStop = 5 Boundary value limit reached
        //  ReasonToStop = 6 Matrix Singular even though Q added

        public static void EndupManxcat(int ReasonToStop, Hotsun.WriteSolutionSignature WriteSolution,
                                        Hotsun.CalcfgSignature Calcfg,
                                        Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
        {
            if (Hotsun.isaved != 0)
            {
                // Copy residuals func from BestSolution to CurrentSolution -- This Logic NOT IMPLEMENTED in Current Manxcat
                Hotsun.zerocr = Hotsun.zeromn;
                if (Hotsun.isaved == 1)
                    RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestSolution);
            }
            // Hotsun.tsolve and Hotsun.teigen set in StopTimer
            Hotsun.tcalcfg = SALSAUtility.SubDurations[2];
            SALSAUtility.EndTiming();
            Hotsun.TotalTimeUsed = SALSAUtility.HPDuration;

            SALSAUtility.SALSAStatus(ResultDirectoryName,
                                     "End Loop " + Hotsun.InitializationLoopCount.ToString() + " Chisq "
                                     + (ChisqPrintConstant*Hotsun.zerocr).ToString("F3"));

            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
            {
                SALSAUtility.SALSAPrint(1, "\n-------------------"
                                           + SALSAUtility.endTime.ToLocalTime() + "\nIterations "
                                           + Hotsun.numit.ToString() + " Chisq "
                                           + (ChisqPrintConstant*Hotsun.zerocr).ToString("F3") + " Q "
                                           + Hotsun.Q.ToString("F3") + " Qlow "
                                           + Hotsun.Qlow.ToString("F3") + " Qhigh "
                                           + Hotsun.Qhigh.ToString("F3") + " Qgood "
                                           + Hotsun.Qgood.ToString("F3") + " Trace Q "
                                           + Hotsun.ChisqMatrixTrace.ToString("F3") + " Norm Q " +
                                           Hotsun.ChisqMatrixNorm.ToString("F3"));
                string bnderrlabel = "";

                if (Hotsun.bnderr != 0)
                    bnderrlabel = " Boundary Violations " + Hotsun.bnderr.ToString() + " ";
                string EndReason = "Correct Convergence";

                if (ReasonToStop == 1)
                    EndReason = "Iteration Limit";

                if (ReasonToStop == 3)
                    EndReason = "Time Cut";

                if (ReasonToStop == 4)
                    EndReason = "No Progress in " + Hotsun.nbadgo.ToString() + " Iterations";

                if (ReasonToStop == 5)
                    EndReason = "Boundary Violations in Calcfg";

                if (ReasonToStop == 6)
                    EndReason = "Matrix Singular";

                if (ReasonToStop == 7)
                    EndReason = "User End";

                // Output parameter values if not too many
                string parametervalues = "";

                if ((!Hotsun.DecomposeParameters) && (Hotsun.npar < 20))
                {
                    parametervalues = "\nParameters:";

                    for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                    {
                        for (int LocalVectorIndex = 0;
                             LocalVectorIndex < Hotsun.ParameterVectorDimension;
                             LocalVectorIndex++)
                        {
                            double tmp = Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
                            parametervalues += " " + tmp.ToString("E4");
                        }
                    }
                }

                SALSAUtility.SALSAPrint(1, "Loop "
                                           + Hotsun.InitializationLoopCount.ToString() + " "
                                           + bnderrlabel + EndReason + " CG Failures "
                                           + Hotsun.TotalCGFailures + " Iterations "
                                           + Hotsun.TotalCGIterations + " Eigenvalue Failures "
                                           + Hotsun.TotalPowerFailures + " Iterations "
                                           + Hotsun.TotalPowerIterations + " Search Iterations "
                                           + Hotsun.TotalSearchIterations + " Total Time "
                                           + Hotsun.TotalTimeUsed.ToString("F1") + " Calcfg Time " +
                                           Hotsun.tcalcfg.ToString("F1") +
                                           " Per Iter Solve " +
                                           (SALSAUtility.SubDurations[0]/Hotsun.TotalCGIterations).ToString("F1")
                                           + " Per Iter Eigen " +
                                           (SALSAUtility.SubDurations[1]/Hotsun.TotalPowerIterations).ToString("F1") +
                                           parametervalues);
                Hotsun.tsolve = 0.0;
                Hotsun.teigen = 0.0;
                Hotsun.tcalcfg = 0.0;
                Hotsun.NumberofCGIterations = -2;
                Hotsun.EigenvalueIndicator1 = 0;
                Hotsun.EigenvalueIndicator2 = 0;
            }

            // Output Chisqcomponent
            if (SALSAUtility.NumberFixedPoints > 0)
                CalculateChisqComponents(Calcfg, GlobalMatrixVectorProduct, Hotsun.CurrentSolution, ChisqPrintConstant);

            //      Process Iteration over Initialization parameters
            if (Hotsun.InitializationLoopCount == 0)
            {
                CopySolution(Hotsun.BestLoopedSolution, Hotsun.CurrentSolution);
                Hotsun.BestChisqLoop = 0;
                if (Hotsun.InitializationLoopCount != (Hotsun.InitializationLoops - 1))
                    WriteOutParameters(WriteSolution, Calcfg, GlobalMatrixVectorProduct);
            }
            else
            {
                if (Hotsun.CurrentSolution.Chisquared < Hotsun.BestLoopedSolution.Chisquared)
                {
                    CopySolution(Hotsun.BestLoopedSolution, Hotsun.CurrentSolution);
                    Hotsun.BestChisqLoop = Hotsun.InitializationLoopCount;

                    if (Hotsun.InitializationLoopCount != (Hotsun.InitializationLoops - 1))
                        WriteOutParameters(WriteSolution, Calcfg, GlobalMatrixVectorProduct);
                }
            }
            Hotsun.InitLoopChisq[Hotsun.InitializationLoopCount] = Hotsun.CurrentSolution.Chisquared;
            Hotsun.InitializationLoopCount++;

            if (Hotsun.InitializationLoopCount < Hotsun.InitializationLoops)
                return;
            WriteOutParameters(WriteSolution, Calcfg, GlobalMatrixVectorProduct);
        }

        // End EndupMancat(int ReasonToStop)

        public static void WriteOutParameters(Hotsun.WriteSolutionSignature WriteSolution, Hotsun.CalcfgSignature Calcfg,
                                              Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
        {
            //  Set up really best Solution
            RestoreBestSolution(Hotsun.CurrentSolution, Hotsun.BestLoopedSolution);

            Hotsun.HotsunComment = ""; // Build final comment here

            //  Output Information on Initialization Loop
            if (Hotsun.InitializationLoopCount == Hotsun.InitializationLoops)
            {
                if (Hotsun.InitializationLoops > 1)
                {
                    string firstcomment = "Best Initial Condition " + Hotsun.BestChisqLoop.ToString() + " out of " +
                                          Hotsun.InitializationLoops.ToString();
                    if (SALSAUtility.MPI_Rank == 0)
                        SALSAUtility.SALSAPrint(1, firstcomment);

                    string ChisqList = "Chisq List";

                    for (int localloop = 0; localloop < Hotsun.InitializationLoops; localloop++)
                    {
                        ChisqList += " " + (ChisqPrintConstant*Hotsun.InitLoopChisq[localloop]).ToString("F3");
                    }
                    firstcomment += "\n" + ChisqList;
                    if (SALSAUtility.MPI_Rank == 0)
                        SALSAUtility.SALSAPrint(1, ChisqList);
                    Hotsun.HotsunComment = firstcomment;
                }
                else
                {
                    Hotsun.HotsunComment = (Hotsun.CurrentSolution.Chisquared*ChisqPrintConstant).ToString("F3");
                }
                // Output Chisqcomponent
                if (SALSAUtility.NumberFixedPoints > 0)
                    CalculateChisqComponents(Calcfg, GlobalMatrixVectorProduct, Hotsun.CurrentSolution,
                                             ChisqPrintConstant);
            }

            //  Final Initialization Loop -- save results
            if (Hotsun.errcal > 0)
            {
                // In this version we assume that we cannot invert Chisq matrix so we estimate errors from Chisq matrix -- not its inverse
                CalculateParameterErrors(Hotsun.CurrentSolution);
                Hotsun.errcal = 1;
            }
            else
            {
                SALSABLAS.zrword(Hotsun.perr);
                Hotsun.errcal = 0;
            }
            WriteCoordinates(1, WriteSolution);
            return;
        }

        // End WriteOutParameters

        public static void WriteCoordinates(int Outputtype, Hotsun.WriteSolutionSignature WriteSolution)
        {
            string outputfiletype = "";
            if (Outputtype == 0)
                outputfiletype = "RESTART";

            //  Output Parameter Values and Errors
            if (!GlobalParameterSet)
            {
                GlobalParameterSet = true;
                MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
            }
            MakeVectorGlobal(Hotsun.perr, Hotsun.UtilityGlobalVector1);

            if (SALSAUtility.MPI_Rank == 0)
            {
                // Set Output File Name and write file with errors
                string labelfilename = Configuration.ReducedVectorOutputFileName;

                // Todo. removing this as it's unnecessary and buggy in Linux (when run with Mono)
//                if (!labelfilename.Contains(":") && !labelfilename.Contains("$"))
//                    ActualOutputFileName = ManxcatCentral.ResultDirectoryName + "\\" + labelfilename;
//                else
                ActualOutputFileName = labelfilename;

                WriteSolution(ActualOutputFileName + outputfiletype, 0, Hotsun.GlobalParameter,
                              Hotsun.UtilityGlobalVector1);

                // Todo. removing this as it's unnecessary and buggy in Linux (when run with Mono)
//                if (!labelfilename.Contains(":") && !labelfilename.Contains("$"))
//                    ActualOutputFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + labelfilename;
//                else
                {
                    int slashpos = labelfilename.LastIndexOf(Path.DirectorySeparatorChar);
                    string endpart = labelfilename.Substring(slashpos + 1);
                    labelfilename = labelfilename.Remove(slashpos + 1);
                    ActualOutputFileName = labelfilename + "SIMPLE" + endpart;
                }
                WriteSolution(ActualOutputFileName + outputfiletype, 1, Hotsun.GlobalParameter,
                              Hotsun.UtilityGlobalVector1);
            }
            return;
        }

        public static void launchQlimits(Hotsun.FindQlimitsSignature FindQlimits)
        {
            if ((Hotsun.numit == 0) || (Hotsun.numit >= IterationtorecalculateQLimits))
            {
                IterationtorecalculateQLimits = Hotsun.numit + Hotsun.QLimitscalculationInterval;
                Hotsun.FullSecondDerivative = false;

                if (Configuration.FullSecondDerivativeOption == 1)
                    Hotsun.FullSecondDerivative = true;

                Hotsun.UseDiagonalScaling = Hotsun.UseDiagonalScalinginSolvers;

                if (Hotsun.UseDiagonalScaling)
                    SetupDiagonalScaling(Hotsun.CurrentSolution);
                Hotsun.AddMarquardtQDynamically = false;

                SALSAUtility.StartSubTimer(1);
                FindQlimits(Hotsun.CurrentSolution, ref Hotsun.Qhigh, ref Hotsun.Qlow, ref Hotsun.EigenvalueIndicator1,
                            ref Hotsun.EigenvalueIndicator2);
                Hotsun.teigen = SALSAUtility.StopSubTimer(1);

                Hotsun.IterationforQlimit = Hotsun.CurrentSolution.IterationCalculated;
            }
            return;
        }

        // End launchQlimits

        public static void CalculateChisqComponents(Hotsun.CalcfgSignature Calcfg,
                                                    Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct,
                                                    Desertwind FromSolution, double MultiplicationFactor)
        {
            int savecomponentflag = SALSAUtility.chisqcomponent;
            SALSAUtility.chisqcomponent = 1;
            double chisq1 = MultiplicationFactor*StandaAloneChisq(Calcfg, GlobalMatrixVectorProduct, FromSolution);
            SALSAUtility.chisqcomponent = 2;
            double chisq2 = MultiplicationFactor*StandaAloneChisq(Calcfg, GlobalMatrixVectorProduct, FromSolution);
            SALSAUtility.chisqcomponent = 3;
            double chisq3 = MultiplicationFactor*StandaAloneChisq(Calcfg, GlobalMatrixVectorProduct, FromSolution);
            SALSAUtility.chisqcomponent = 4;
            double chisq4 = MultiplicationFactor*StandaAloneChisq(Calcfg, GlobalMatrixVectorProduct, FromSolution);
            string componentcomment = " V-V " + chisq1.ToString("F3") + " V-F " + chisq2.ToString("F3") + " F-V " +
                                      chisq3.ToString("F3") + " F-F " + chisq4.ToString("F3");
            if (SALSAUtility.MPI_Rank == 0)
                SALSAUtility.SALSAPrint(1, componentcomment);
            SALSAUtility.chisqcomponent = savecomponentflag;
            Hotsun.HotsunComment += "\n" + componentcomment;
            return;
        }

        // End CalculateChisqComponents()

        public static double StandaAloneChisq(Hotsun.CalcfgSignature Calcfg,
                                              Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct,
                                              Desertwind FromSolution)
        {
            double savechisq = Hotsun.zerocr;
            Hotsun.zerocr = 0.0;
            bool SaveFullSecondDerivative = Hotsun.FullSecondDerivative;
            Hotsun.FullSecondDerivative = true;
            bool SaveUseDiagonalScaling = Hotsun.UseDiagonalScaling;
            Hotsun.UseDiagonalScaling = false;
            bool SaveAddMarquardtQDynamically = Hotsun.AddMarquardtQDynamically;
            Hotsun.AddMarquardtQDynamically = false;
            bool saveCalcFixedCrossFixed = SALSAUtility.CalcFixedCrossFixed;
            SALSAUtility.CalcFixedCrossFixed = true;
            if (SALSAUtility.StoredDistanceOption == 3)
                SALSAUtility.CalcFixedCrossFixed = false;

            ZeroSolution(Hotsun.SearchSolution1);
            int Numberparms = Hotsun.Number_VectorParameters;
            if (Hotsun.DecomposeParameters)
                Numberparms = SALSAUtility.PointCount_Process;
            SALSABLAS.CopyVector(Hotsun.SearchSolution1.param, FromSolution.param, 0, Numberparms);
            MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);

            bool localviolat = Calcfg(Hotsun.SearchSolution1);
            double toreturn = Hotsun.zerocr;

            SALSAUtility.CalcFixedCrossFixed = saveCalcFixedCrossFixed;
            Hotsun.FullSecondDerivative = SaveFullSecondDerivative;
            Hotsun.UseDiagonalScaling = SaveUseDiagonalScaling;
            Hotsun.AddMarquardtQDynamically = SaveAddMarquardtQDynamically;
            Hotsun.zerocr = savechisq;
            MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
            return toreturn;
        }

        // End standalone chisq calculation

        public static void DoDerivTest(Hotsun.CalcfgSignature Calcfg,
                                       Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
        {
            if (SALSAUtility.MPI_Size > 1)
                return;
            double savechisq = Hotsun.zerocr;

            int PointtoVary1 = Hotsun.Number_VectorParameters - 25;
            if (PointtoVary1 >= SALSAUtility.NumberVariedPoints)
                PointtoVary1 = SALSAUtility.NumberVariedPoints - 5;
            if (PointtoVary1 < 0)
                PointtoVary1 = 5;

            int IndextoVary1 = 2;
            if (IndextoVary1 > Hotsun.ParameterVectorDimension)
                IndextoVary1 = 0;

            int PointtoVary2 = 25;
            if (PointtoVary2 >= SALSAUtility.NumberVariedPoints)
                PointtoVary2 = 4;

            int IndextoVary2 = 2;
            if (IndextoVary2 > Hotsun.ParameterVectorDimension)
                IndextoVary2 = 0;
            if (Hotsun.npar < 10)
            {
                PointtoVary1 = 6;
                PointtoVary2 = 3;
            }
            int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[PointtoVary1];
            PointtoVary1 = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex1];
            int OriginalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[PointtoVary2];
            PointtoVary2 = SALSAUtility.OriginalPointtoUsedPointMap[OriginalPointIndex2];

            SALSAUtility.SALSAPrint(1, "Point1 " + PointtoVary1.ToString() + "," + IndextoVary1.ToString() +
                                       " Point2 " + PointtoVary2.ToString() + "," + IndextoVary2.ToString());
            Hotsun.doingderivtest = true;
            bool SaveFullSecondDerivative = Hotsun.FullSecondDerivative;
            Hotsun.FullSecondDerivative = true;
            bool SaveUseDiagonalScaling = Hotsun.UseDiagonalScaling;
            Hotsun.UseDiagonalScaling = false;
            bool SaveAddMarquardtQDynamically = Hotsun.AddMarquardtQDynamically;
            Hotsun.AddMarquardtQDynamically = false;

            ZeroSolution(Hotsun.SearchSolution1);
            int Numberparms = Hotsun.Number_VectorParameters;

            if (Hotsun.DecomposeParameters)
                Numberparms = SALSAUtility.PointCount_Global;
            SALSABLAS.CopyVector(Hotsun.SearchSolution1.param, Hotsun.CurrentSolution.param, 0, Numberparms);
            Hotsun.zerocr = 0.0;
            bool localviolat = Calcfg(Hotsun.SearchSolution1);
            double startingpoint = Hotsun.zerocr;

            double variableincrement = 0.0001*Math.Abs(Hotsun.CurrentSolution.xshift[PointtoVary1][IndextoVary1]);

            double FirstDerivative1 = 2.0*Hotsun.SearchSolution1.first[PointtoVary1][IndextoVary1];
            double FirstDerivative2 = 2.0*Hotsun.SearchSolution1.first[PointtoVary2][IndextoVary2];
            double SecondDerivative11;

            if (Hotsun.fullmatrixset)
                SecondDerivative11 = 2.0*
                                     Hotsun.SearchSolution1.ExactFullMatrix[PointtoVary1, PointtoVary1][
                                         IndextoVary1, IndextoVary1];
            else
                SecondDerivative11 = 2.0*
                                     Hotsun.SearchSolution1.ExactDiagonalofMatrix[PointtoVary1][
                                         IndextoVary1, IndextoVary1];
            MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);

            // Calculate Second Derivative
            SALSABLAS.zrword(Hotsun.UtilityLocalVector1);
            Hotsun.UtilityLocalVector1[PointtoVary1][IndextoVary1] = 1.0;
            MakeVectorGlobal(Hotsun.UtilityLocalVector1, Hotsun.UtilityGlobalVector1);
            GlobalMatrixVectorProduct(Hotsun.UtilityLocalVector2, Hotsun.SearchSolution1, true, Hotsun.GlobalParameter,
                                      Hotsun.UtilityGlobalVector1);
            double SecondDerivative21 = 2.0*Hotsun.UtilityLocalVector2[PointtoVary2][IndextoVary2];
            double SecondDerivative11A = 2.0*Hotsun.UtilityLocalVector2[PointtoVary1][IndextoVary1];
            SALSABLAS.zrword(Hotsun.UtilityLocalVector1);
            Hotsun.UtilityLocalVector1[PointtoVary2][IndextoVary2] = 1.0;
            MakeVectorGlobal(Hotsun.UtilityLocalVector1, Hotsun.UtilityGlobalVector1);
            GlobalMatrixVectorProduct(Hotsun.UtilityLocalVector2, Hotsun.SearchSolution1, true, Hotsun.GlobalParameter,
                                      Hotsun.UtilityGlobalVector1);
            double SecondDerivative12 = 2.0*Hotsun.UtilityLocalVector2[PointtoVary1][IndextoVary1];

            if (Hotsun.npar < 10)
            {
                for (int i1 = 0; i1 < Hotsun.npar; i1++)
                {
                    string linetoprint = "Row " + i1.ToString() + " ";
                    for (int i2 = 0; i2 < Hotsun.npar; i2++)
                    {
                        double a;
                        if (Hotsun.fullmatrixset)
                            a = 2.0*Hotsun.SearchSolution1.ExactFullMatrix[i1, i2][IndextoVary1, IndextoVary1];
                        else
                            a = 2.0*Hotsun.SearchSolution1.ExactDiagonalofMatrix[i1][IndextoVary1, IndextoVary1];
                        linetoprint += a.ToString("E4") + " * ";
                    }
                    SALSAUtility.SALSAPrint(1, linetoprint);
                }
            }

            double variablestep = 0.0;
            SALSAUtility.SALSAPrint(1, " Variable " + PointtoVary1 + " " + IndextoVary1 + " Value "
                                       + Hotsun.SearchSolution1.param[PointtoVary1][IndextoVary1].ToString("E4") +
                                       " Increment " + variableincrement.ToString("E4"));
            for (int step = 0; step < 2; step++)
            {
                Hotsun.SearchSolution1.param[PointtoVary1][IndextoVary1] += variableincrement;
                MakeVectorGlobal(Hotsun.SearchSolution1.param, Hotsun.GlobalParameter);
                Hotsun.zerocr = 0.0;
                ZeroSolution(Hotsun.SearchSolution1);
                localviolat = Calcfg(Hotsun.SearchSolution1);
                variablestep += variableincrement;
                double change = (Hotsun.zerocr - startingpoint);
                double newzerocr = Hotsun.zerocr;
                double NumericalFirstDerivative1 = change/variablestep;
                double NewFirstDerivative1 = 2.0*Hotsun.SearchSolution1.first[PointtoVary1][IndextoVary1];
                double NumericalSecondDerivative11 = (NewFirstDerivative1 - FirstDerivative1)/variablestep;
                    // Diagonal Derivative
                double NumericalSecondDerivative21 = (2.0*Hotsun.SearchSolution1.first[PointtoVary2][IndextoVary2] -
                                                      FirstDerivative2)/variablestep; // Off-Diagonal Derivative
                double FuncChangeestimate = FirstDerivative1*variablestep +
                                            0.5*SecondDerivative11*variablestep*variablestep;
                SALSAUtility.SALSAPrint(1, step + " step "
                                           + variablestep.ToString("E4") + " chi "
                                           + newzerocr.ToString("E4") + " Fact "
                                           + change.ToString("E4") + " Fest "
                                           + FuncChangeestimate.ToString("E4") + " D1est "
                                           + NumericalFirstDerivative1.ToString("E4") + " D1act "
                                           + FirstDerivative1.ToString("E4") + " S11est "
                                           + NumericalSecondDerivative11.ToString("E4") + " S11act "
                                           + SecondDerivative11.ToString("E4") + " S11actA "
                                           + SecondDerivative11A.ToString("E4") + " S21est "
                                           + NumericalSecondDerivative21.ToString("E4") + " S21act "
                                           + SecondDerivative21.ToString("E4") + " S12act " +
                                           SecondDerivative12.ToString("E4"));
            }

            Hotsun.FullSecondDerivative = SaveFullSecondDerivative;
            Hotsun.UseDiagonalScaling = SaveUseDiagonalScaling;
            Hotsun.AddMarquardtQDynamically = SaveAddMarquardtQDynamically;
            Hotsun.doingderivtest = false;
            Hotsun.zerocr = savechisq;
            MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
            return;
        }

        // End DoDerivTest()


        //  Set up 2D Double MPI Packet
        public static void SetupMPIPacket(out MPI2DDoubleVectorPacket TogoVector)
        {
            TogoVector = new MPI2DDoubleVectorPacket(SALSAUtility.PointCount_Largest, Hotsun.ParameterVectorDimension);
            TogoVector.FirstPoint = SALSAUtility.PointStart_Process;
            TogoVector.NumberofPoints = SALSAUtility.PointCount_Process;

            if (SALSAUtility.PointCount_Process == SALSAUtility.PointCount_Largest)
                return;

            for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
            {
                TogoVector.Marray[SALSAUtility.PointCount_Largest - 1][LocalVectorIndex] = 0.0;
            }
        }

        // End Setup2DDoubleMPIPacket

        //  Set up 1D String MPI Packet
        public static void SetupMPIPacket(out MPI1DStringVectorPacket TogoVector)
        {
            TogoVector = new MPI1DStringVectorPacket(SALSAUtility.PointCount_Largest, Hotsun.ParameterVectorDimension);
            TogoVector.FirstPoint = SALSAUtility.PointStart_Process;
            TogoVector.NumberofPoints = SALSAUtility.PointCount_Process;

            if (SALSAUtility.PointCount_Process == SALSAUtility.PointCount_Largest)
                return;
            TogoVector.Marray[SALSAUtility.PointCount_Largest - 1] = "";
        }

        // End Setup1DStringMPIPacket


        // Globalize distributed local copies with 2D Double
        public static void MakeVectorGlobal(double[][] DistributedVector, double[][] GlobalVector)
        {
            if ((SALSAUtility.MPI_Size <= 1) || (!Hotsun.DecomposeParameters))
            {
                SALSABLAS.CopyVector(GlobalVector, DistributedVector, 0, Hotsun.Number_VectorParameters);
            }
            else
            {
                SALSABLAS.CopyVector(TogoDistributed2DDoubleVector.Marray, DistributedVector, 0,
                                     SALSAUtility.PointCount_Process);

                SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

                for (int MPICommunicationSteps = 0;
                     MPICommunicationSteps < SALSAUtility.MPI_Size;
                     MPICommunicationSteps++)
                {
                    if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
                        MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector,
                                                                            TogoDistributed2DDoubleVector);
                    SALSAUtility.MPI_communicator.Broadcast(ref FromAfar2DDoubleVector, MPICommunicationSteps);
                    SALSABLAS.CopyVector(GlobalVector, FromAfar2DDoubleVector.Marray, FromAfar2DDoubleVector.FirstPoint,
                                         FromAfar2DDoubleVector.NumberofPoints);
                } // End loop over MPIrankcount
                SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
            }
        }

        // End MakeVectorGlobal() 2D Double

        // Globalize distributed local copies with 1D String
        public static void MakeVectorGlobal(string[] DistributedVector, string[] GlobalVector)
        {
            if (SALSAUtility.MPI_Size <= 1)
            {
                SALSABLAS.CopyVector(GlobalVector, DistributedVector, 0, Hotsun.Number_VectorParameters);
            }
            else
            {
                SALSABLAS.CopyVector(TogoDistributed1DStringVector.Marray, DistributedVector, 0,
                                     SALSAUtility.PointCount_Process);

                SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

                for (int MPICommunicationSteps = 0;
                     MPICommunicationSteps < SALSAUtility.MPI_Size;
                     MPICommunicationSteps++)
                {
                    if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
                        MPI1DStringVectorPacket.CopyMPI1DStringVectorPacket(FromAfar1DStringVector,
                                                                            TogoDistributed1DStringVector);
                    SALSAUtility.MPI_communicator.Broadcast(ref FromAfar1DStringVector, MPICommunicationSteps);
                    SALSABLAS.CopyVector(GlobalVector, FromAfar1DStringVector.Marray, FromAfar1DStringVector.FirstPoint,
                                         FromAfar1DStringVector.NumberofPoints);
                } // End loop over MPIrankcount
                SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
            }
        }

        // End MakeVectorGlobal() with 1D String

        //  Set up GLOBAL arrays Hotsun.diag and Hotsun.sqdginv
        public static void SetupDiagonalScaling(Desertwind Solution)
        {
            if (!Hotsun.DecomposeParameters)
            {
                for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                {
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        double tmp;

                        if (Hotsun.fullmatrixset)
                        {
                            if (Hotsun.FullSecondDerivative)
                                tmp =
                                    Math.Abs(
                                        Solution.ExactFullMatrix[LongIndex, LongIndex][
                                            LocalVectorIndex, LocalVectorIndex]);
                            else
                                tmp =
                                    Math.Abs(
                                        Solution.FullMatrix[LongIndex, LongIndex][LocalVectorIndex, LocalVectorIndex]);
                        }
                        else
                        {
                            if (Hotsun.FullSecondDerivative)
                                tmp =
                                    Math.Abs(
                                        Solution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex]);
                            else
                                tmp = Math.Abs(Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex]);
                        }

                        Hotsun.diag[LongIndex][LocalVectorIndex] = tmp;
                        if (Hotsun.FixedParameter[LongIndex][LocalVectorIndex])
                            Hotsun.sqdginv[LongIndex][LocalVectorIndex] = 0.0;
                        else
                            Hotsun.sqdginv[LongIndex][LocalVectorIndex] = 1.0/Math.Sqrt(tmp);
                    }
                }
                return;
            }

            // Parallel Local Calculation of Diagonal Scaling Elements
            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LocalToProcessIndex = beginpoint;
                                      LocalToProcessIndex < indexlen + beginpoint;
                                      LocalToProcessIndex++)
                                 {
                                     int GlobalIndex = LocalToProcessIndex + SALSAUtility.PointStart_Process;
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         double tmp;

                                         if (Hotsun.fullmatrixset)
                                         {
                                             if (Hotsun.FullSecondDerivative)
                                                 tmp =
                                                     Math.Abs(
                                                         Solution.ExactFullMatrix[
                                                             LocalToProcessIndex, LocalToProcessIndex][
                                                                 LocalVectorIndex, LocalVectorIndex]);
                                             else
                                                 tmp =
                                                     Math.Abs(
                                                         Solution.FullMatrix[LocalToProcessIndex, LocalToProcessIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                         }
                                         else
                                         {
                                             if (Hotsun.FullSecondDerivative)
                                                 tmp =
                                                     Math.Abs(
                                                         Solution.ExactDiagonalofMatrix[LocalToProcessIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                             else
                                                 tmp =
                                                     Math.Abs(
                                                         Solution.DiagonalofMatrix[LocalToProcessIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                         }

                                         TogoDiagVector.Marray[LocalToProcessIndex][LocalVectorIndex] = tmp;
                                         if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
                                             TogoSqDgInvVector.Marray[LocalToProcessIndex][LocalVectorIndex] = 0.0;
                                         else
                                             TogoSqDgInvVector.Marray[LocalToProcessIndex][LocalVectorIndex] = 1.0/
                                                                                                               Math.Sqrt
                                                                                                                   (tmp);
                                     }
                                 }
                             }); // End loop over Point dependent quantities

            // Convert into Global arrays
            if (SALSAUtility.MPI_Size <= 1)
            {
                // No MPI
                SALSABLAS.CopyVector(Hotsun.diag, TogoDiagVector.Marray, TogoDiagVector.FirstPoint,
                                     TogoDiagVector.NumberofPoints);
                SALSABLAS.CopyVector(Hotsun.sqdginv, TogoSqDgInvVector.Marray, TogoSqDgInvVector.FirstPoint,
                                     TogoSqDgInvVector.NumberofPoints);
            }
            else
            {
                for (int MPICommunicationSteps = 0;
                     MPICommunicationSteps < SALSAUtility.MPI_Size;
                     MPICommunicationSteps++)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPISENDRECEIVETiming);

                    if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
                        MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector, TogoDiagVector);

                    SALSAUtility.MPI_communicator.Broadcast(ref FromAfar2DDoubleVector, MPICommunicationSteps);
                    SALSABLAS.CopyVector(Hotsun.diag, FromAfar2DDoubleVector.Marray, FromAfar2DDoubleVector.FirstPoint,
                                         FromAfar2DDoubleVector.NumberofPoints);

                    if (MPICommunicationSteps == SALSAUtility.MPI_Rank)
                        MPI2DDoubleVectorPacket.CopyMPI2DDoubleVectorPacket(FromAfar2DDoubleVector, TogoSqDgInvVector);
                    SALSAUtility.MPI_communicator.Broadcast(ref FromAfar2DDoubleVector, MPICommunicationSteps);
                    SALSABLAS.CopyVector(Hotsun.sqdginv, FromAfar2DDoubleVector.Marray,
                                         FromAfar2DDoubleVector.FirstPoint, FromAfar2DDoubleVector.NumberofPoints);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPISENDRECEIVETiming);
                } // End loop over MPIrankcount
            }
        }

        // End SetupDiagonalScaling()

        //  Approximate error estimate in paramters
        public static void CalculateParameterErrors(Desertwind Solution)
        {
            // Parallel Calculation of Parameter Errors
            double NormalizeSQRTChisq = Math.Sqrt(Hotsun.zeromn/(Hotsun.ndata - Hotsun.npar));

            if (!Hotsun.DecomposeParameters)
            {
                for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                {
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        double tmp;

                        if (Hotsun.fullmatrixset)
                        {
                            tmp = Solution.FullMatrix[LongIndex, LongIndex][LocalVectorIndex, LocalVectorIndex];
                        }
                        else
                        {
                            tmp = Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex];
                        }
                        Hotsun.perr[LongIndex][LocalVectorIndex] = NormalizeSQRTChisq/Math.Sqrt(tmp);
                    }
                }
                return;
            }

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         double tmp;

                                         if (Hotsun.fullmatrixset)
                                         {
                                             tmp =
                                                 Solution.FullMatrix[LongIndex, LongIndex][
                                                     LocalVectorIndex, LocalVectorIndex];
                                         }
                                         else
                                         {
                                             tmp =
                                                 Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex
                                                     ];
                                         }
                                         Hotsun.perr[LongIndex][LocalVectorIndex] = NormalizeSQRTChisq/
                                                                                    Math.Sqrt(
                                                                                        Solution.DiagonalofMatrix[
                                                                                            LongIndex][
                                                                                                LocalVectorIndex,
                                                                                                LocalVectorIndex]);
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // End CalculateParameterErrors()

        public static void AddSubtractMarquardtFactor(Desertwind Solution, double factor)
        {
            if (Hotsun.AddMarquardtQDynamically)
                return;

            if (!Hotsun.DecomposeParameters)
            {
                for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                {
                    int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        if (Hotsun.fullmatrixset)
                        {
                            Solution.FullMatrix[LongIndex, LongIndex][LocalVectorIndex, LocalVectorIndex] += factor*
                                                                                                             Hotsun.Q*
                                                                                                             Hotsun.diag
                                                                                                                 [
                                                                                                                     GlobalIndex
                                                                                                                 ][
                                                                                                                     LocalVectorIndex
                                                                                                                 ];
                            Solution.ExactFullMatrix[LongIndex, LongIndex][LocalVectorIndex, LocalVectorIndex] +=
                                factor*Hotsun.Q*Hotsun.diag[GlobalIndex][LocalVectorIndex];
                        }
                        else
                        {
                            Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex] += factor*Hotsun.Q*
                                                                                                        Hotsun.diag[
                                                                                                            GlobalIndex]
                                                                                                            [
                                                                                                                LocalVectorIndex
                                                                                                            ];
                            Solution.ExactDiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex] += factor*
                                                                                                             Hotsun.Q*
                                                                                                             Hotsun.diag
                                                                                                                 [
                                                                                                                     GlobalIndex
                                                                                                                 ][
                                                                                                                     LocalVectorIndex
                                                                                                                 ];
                        }
                    }
                }
                return;
            }

            // Parallel Addition (factor = 1.0) or Subtraction (factor = -1.0) of Marquardt A
            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         if (Hotsun.fullmatrixset)
                                         {
                                             Solution.FullMatrix[LongIndex, LongIndex][
                                                 LocalVectorIndex, LocalVectorIndex] += factor*Hotsun.Q*
                                                                                        Hotsun.diag[GlobalIndex][
                                                                                            LocalVectorIndex];
                                             Solution.ExactFullMatrix[LongIndex, LongIndex][
                                                 LocalVectorIndex, LocalVectorIndex] += factor*Hotsun.Q*
                                                                                        Hotsun.diag[GlobalIndex][
                                                                                            LocalVectorIndex];
                                         }
                                         else
                                         {
                                             Solution.DiagonalofMatrix[LongIndex][LocalVectorIndex, LocalVectorIndex] +=
                                                 factor*Hotsun.Q*Hotsun.diag[GlobalIndex][LocalVectorIndex];
                                             Solution.ExactDiagonalofMatrix[LongIndex][
                                                 LocalVectorIndex, LocalVectorIndex] += factor*Hotsun.Q*
                                                                                        Hotsun.diag[GlobalIndex][
                                                                                            LocalVectorIndex];
                                         }
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // End AddinMarquardtFactor

        public static void SteepestDescentSolution(Desertwind Solution,
                                                   Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
        {
            if (!Hotsun.DecomposeParameters)
            {
                // Sequential Calculation of Steepest Descent
                for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                {
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        double tmp = Hotsun.sqdginv[LongIndex][LocalVectorIndex];
                        Solution.xshift[LongIndex][LocalVectorIndex] = Solution.first[LongIndex][LocalVectorIndex]*tmp;
                    }
                }
            }
            else
            {
                // Parallel Calculation of Steepest Descent
                Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                             (ThreadNo) =>
                                 {
                                     int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                     int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                      SALSAUtility.PointStart_Process;

                                     for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                     {
                                         int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

                                         for (int LocalVectorIndex = 0;
                                              LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                              LocalVectorIndex++)
                                         {
                                             double tmp = Hotsun.sqdginv[GlobalIndex][LocalVectorIndex];
                                             Solution.xshift[LongIndex][LocalVectorIndex] =
                                                 Solution.first[LongIndex][LocalVectorIndex]*tmp;
                                         }
                                     }
                                 }); // End loop over Point dependent quantities
            }

            // Convert xshift into a Global Vector
            MakeVectorGlobal(Solution.xshift, Hotsun.UtilityGlobalVector1);
            Hotsun.UseDiagonalScaling = false;
            Hotsun.AddMarquardtQDynamically = false;

            //  Calculate Normalization
            // First Calculate ChisqMatrix dotted into xshift
            Hotsun.FullSecondDerivative = false;
            GlobalMatrixVectorProduct(Hotsun.UtilityLocalVector1, Solution, false, Hotsun.GlobalParameter,
                                      Hotsun.UtilityGlobalVector1);
            double Fdotxshift = SALSABLAS.VectorScalarProduct(Solution.xshift, Solution.first);
            double xshiftdotGdotxshift = SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);
            double factor = Fdotxshift/xshiftdotGdotxshift;
            SALSABLAS.LinearCombineVector(Solution.xshift, factor, Solution.xshift, 0.0, Solution.xshift);
        }

        // End SteepestDescentSolution()

        //  Find xnorm
        //  If method = 1 take existing xshift
        //  If method = 2 first set xshift from Best - Current Solution
        public static void Findxnorm(int method)
        {
            if (!Hotsun.DecomposeParameters)
            {
                double localnorm = 0.0;

                for (int LongIndex = 0; LongIndex < Hotsun.Number_VectorParameters; LongIndex++)
                {
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        if (Hotsun.FixedParameter[LongIndex][LocalVectorIndex])
                        {
                            Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = 0;
                            continue;
                        }

                        if (method == 2)
                            Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] =
                                Hotsun.BestSolution.param[LongIndex][LocalVectorIndex] -
                                Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
                        double tmp;

                        if (Hotsun.fullmatrixset)
                        {
                            if (Hotsun.FullSecondDerivative)
                                tmp =
                                    Math.Abs(
                                        Hotsun.CurrentSolution.ExactFullMatrix[LongIndex, LongIndex][
                                            LocalVectorIndex, LocalVectorIndex]);
                            else
                                tmp =
                                    Math.Abs(
                                        Hotsun.CurrentSolution.FullMatrix[LongIndex, LongIndex][
                                            LocalVectorIndex, LocalVectorIndex]);
                        }
                        else
                        {
                            if (Hotsun.FullSecondDerivative)
                                tmp =
                                    Math.Abs(
                                        Hotsun.CurrentSolution.ExactDiagonalofMatrix[LongIndex][
                                            LocalVectorIndex, LocalVectorIndex]);
                            else
                                tmp =
                                    Math.Abs(
                                        Hotsun.CurrentSolution.DiagonalofMatrix[LongIndex][
                                            LocalVectorIndex, LocalVectorIndex]);
                        }
                        localnorm += Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex]*
                                     Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex]*tmp;
                    }
                }
                Hotsun.xnorm = localnorm;
                return;
            }

            var FindxNorm = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 double localnorm = 0.0;
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
                                         {
                                             Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] = 0;
                                             continue;
                                         }

                                         if (method == 2)
                                             Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex] =
                                                 Hotsun.BestSolution.param[LongIndex][LocalVectorIndex] -
                                                 Hotsun.CurrentSolution.param[LongIndex][LocalVectorIndex];
                                         double tmp;

                                         if (Hotsun.fullmatrixset)
                                         {
                                             if (Hotsun.FullSecondDerivative)
                                                 tmp =
                                                     Math.Abs(
                                                         Hotsun.CurrentSolution.ExactFullMatrix[LongIndex, LongIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                             else
                                                 tmp =
                                                     Math.Abs(
                                                         Hotsun.CurrentSolution.FullMatrix[LongIndex, LongIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                         }
                                         else
                                         {
                                             if (Hotsun.FullSecondDerivative)
                                                 tmp =
                                                     Math.Abs(
                                                         Hotsun.CurrentSolution.ExactDiagonalofMatrix[LongIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                             else
                                                 tmp =
                                                     Math.Abs(
                                                         Hotsun.CurrentSolution.DiagonalofMatrix[LongIndex][
                                                             LocalVectorIndex, LocalVectorIndex]);
                                         }
                                         localnorm += Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex]*
                                                      Hotsun.CurrentSolution.xshift[LongIndex][LocalVectorIndex]*tmp;
                                     }
                                 }
                                 FindxNorm.addapoint(ThreadNo, localnorm);
                             }); // End loop over Point dependent quantities

            FindxNorm.sumoverthreadsandmpi();
            Hotsun.xnorm = FindxNorm.Total;
        }

        // End Findxnorm

        // Even if diagonal scaling used in solver, this is not stored in DiagonalofMatrix or first
        public static void FindPredictedChisqChange(Desertwind Solution,
                                                    Hotsun.GlobalMatrixVectorProductSignature GlobalMatrixVectorProduct)
        {
            // Calculate xnorm pred1 and pred2
            Findxnorm(1);

            Hotsun.pred1 = 2.0*SALSABLAS.VectorScalarProduct(Solution.xshift, Solution.first); // Linear Term pred1

            // Convert xshift into a Global Vector
            MakeVectorGlobal(Solution.xshift, Hotsun.UtilityGlobalVector1);
            Hotsun.UseDiagonalScaling = false;
            Hotsun.AddMarquardtQDynamically = false;
            bool saveflag = Hotsun.FullSecondDerivative;
            Hotsun.FullSecondDerivative = false;
            GlobalMatrixVectorProduct(Hotsun.UtilityLocalVector1, Solution, false, Hotsun.GlobalParameter,
                                      Hotsun.UtilityGlobalVector1);
            Hotsun.pred2 = -SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);
            Hotsun.FullSecondDerivative = true;
            GlobalMatrixVectorProduct(Hotsun.UtilityLocalVector1, Solution, true, Hotsun.GlobalParameter,
                                      Hotsun.UtilityGlobalVector1);
            double FullPred = -SALSABLAS.VectorScalarProduct(Solution.xshift, Hotsun.UtilityLocalVector1);

            if (saveflag)
            {
                Hotsun.pred3 = Hotsun.pred2;
                Hotsun.pred2 = FullPred;
            }
            else
            {
                Hotsun.pred3 = FullPred;
            }
            Hotsun.FullSecondDerivative = saveflag;
        }

        // End FindPredictedChisqChange()

        //  CurrentSolution is current configuration wth value derivative
        //  OriginalSolution was previous configuration with value and  derivative
        //  LineFactor * xshift is returned solution that decreases Chisq by an additional Extradecrease in LineIteration iterations

        //  LineSearchMethod    1 No chisq increase. Method should not have been called. Shift is 1
        //  LineSearchMethod    40 + DerivativeSearchMethod Call DerivativeSearch for case where shift < 1; DerivativeSearch is method indicator returned by DerivativeSearch
        //  LineSearchMethod    3 Request to find a shift > 1 but Chisq decreases for shift = 1 so rejected. Return Shift = 1
        //  LineSearchMethod    80 + DerivativeSearchMethod Call DerivativeSearch for case where requested shift > 1 BUT derivative positive at shift of 1; So Shift < 1 returned; DerivativeSearch is method indicator returned by DerivativeSearch
        //  LineSearchMethod    120 + DerivativeSearchMethod Call DerivativeSearch for case where requested shift > 1 AND derivative negative at shift of 1; Shift > 1 possible; DerivativeSearch is method indicator returned by DerivativeSearch

        //  DerivativeSearchMethod 30 + GoldenSectionMethod if Right Derivative < 0 This is a value search between Middle and Right 
        //  DerivativeSearchMethod 20 + GoldenSectionMethod if Left Solution bad This is a value search between Left and Middle
        //  or DerivativeSearchMethod = SimpleDerivativeSearchMethod      search between Middle and Right        
        //  or DerivativeSearchMethod = SimpleDerivativeSearchMethod + 5     search between Right and Middle
        //  or DerivativeSearchMethod = SimpleDerivativeSearchMethod + 10     search between Right and Left

        //  SimpleDerivativeSearchMethod 1 Failure New point higher Take Middle  
        //  SimpleDerivativeSearchMethod 2 Solution is found point with lower chisq and Positive derivative 
        //  SimpleDerivativeSearchMethod 3 Solution is found point possibly with lower chisq -- Cut on position change
        //  SimpleDerivativeSearchMethod 4 Solution is found point possibly with lower chisq -- Cut on value change 
        //  SimpleDerivativeSearchMethod 5 Iteration Cut


        //  GoldenSectionMethod 1   Null Solution handed to routine
        //  GoldenSectionMethod 2   Inconsistent Data   Right Lowest and also Left < Middle
        //  GoldenSectionMethod 3   Inconsistent Data   Left Lowest
        //  GoldenSectionMethod 4   Inconsistent Data   Right Lowest
        //  GoldenSectionMethod 5   Inconsistent Data   Positions inconsistent
        //  GoldenSectionMethod 6   Consistent Data Cut on position values
        //  GoldenSectionMethod 7   Consistent Data Cut on Chisq values
        //  GoldenSectionMethod 8   Consistent Data Iteration Cut 

        public static bool LineSearch(Desertwind CurrentSolution, Desertwind OriginalSolution, out int LineSearchMethod,
                                      out double LineFactor, out double ExtraDecrease,
                                      out int LineIterations, Hotsun.CalcfgSignature Calcfg, double EstimatedLineFactor)
        {
            var Left = new SearchStuff();
            var Middle = new SearchStuff();
            var Right = new SearchStuff();
            var BestFromLineSearch = new SearchStuff();
            bool success;
            int LocalIterations;

            LineFactor = 1.0;
            ExtraDecrease = 0.0;
            LineIterations = 0;
            LineSearchMethod = 0;
            int DerivSearchMethod = 0;

            // BeginningLinePositionSolution; Set to starting point of search (holds param of starting point and value of first derivative there)
            // EndingLinePositionSolution: Set to starting point of search (holds xshift from starting point and value and first derivative of end of line )
            Hotsun.BeginningLinePositionSolution = OriginalSolution;
            Hotsun.EndingLinePositionSolution = CurrentSolution;

            if (EstimatedLineFactor > 0.0)
            {
                // Case of LineFactor < 1 from an initial step that failed
                bool NoChisqIncrease = CurrentSolution.Chisquared <= OriginalSolution.Chisquared;
                SALSAUtility.SynchronizeMPIvariable(ref NoChisqIncrease);

                if (NoChisqIncrease)
                {
                    LineSearchMethod = 1;
                    return false; // Method should not have been called
                }
                double alpha = Math.Min(0.5, 2.0*EstimatedLineFactor);
                CopySolution(Hotsun.SearchSolution1, OriginalSolution);
                CopySolution(Hotsun.SearchSolution3, CurrentSolution);

                ExistingCalcDeriv_LineSearch(0.0, Hotsun.SearchSolution1, ref Left);
                NewCalcDeriv_LineSearch(alpha, Hotsun.SearchSolution2, ref Middle, Calcfg);
                ExistingCalcDeriv_LineSearch(1.0, Hotsun.SearchSolution3, ref Right);

                bool DerivSearch3 = DerivativeSearch(Left, Middle, Right, out DerivSearchMethod, out BestFromLineSearch,
                                                     out LocalIterations, Calcfg);
                ++LocalIterations;
                LineSearchMethod = 40 + DerivSearchMethod;
                SALSAUtility.SynchronizeMPIvariable(ref DerivSearch3);
                success = DerivSearch3;
            }
            else
            {
                // Case of LineFactor > 1 from a good step that can be improved
                bool NoChisqDecrease = CurrentSolution.Chisquared >= OriginalSolution.Chisquared;
                SALSAUtility.SynchronizeMPIvariable(ref NoChisqDecrease);

                if (NoChisqDecrease)
                {
                    LineSearchMethod = 3;
                    return false; // Method should not have been called
                }
                Left.Solution = null;
                CopySolution(Hotsun.SearchSolution1, OriginalSolution);
                CopySolution(Hotsun.SearchSolution2, CurrentSolution);

                ExistingCalcDeriv_LineSearch(0.0, Hotsun.SearchSolution1, ref Middle);
                ExistingCalcDeriv_LineSearch(1.0, Hotsun.SearchSolution2, ref Right);

                bool RightDerivativePositive = Right.deriv > 0;
                SALSAUtility.SynchronizeMPIvariable(ref RightDerivativePositive);

                if (RightDerivativePositive)
                {
                    // Inconsistent derivatives

                    bool DerivSearch1 = DerivativeSearch(Left, Middle, Right, out DerivSearchMethod,
                                                         out BestFromLineSearch, out LocalIterations, Calcfg);
                    LineSearchMethod = 80 + DerivSearchMethod;
                    SALSAUtility.SynchronizeMPIvariable(ref DerivSearch1);
                    success = DerivSearch1;
                }
                else
                {
                    // Correct Sign of derivative -- explore larger values of alpha
                    double alpha = 1.0;

                    while (LineIterations <= 20)
                    {
                        // Iterate over increasing alpha
                        ++LineIterations;
                        alpha += 1.0;
                        Desertwind SaveSolutionReference;

                        if (Left.Solution == null)
                        {
                            SaveSolutionReference = Hotsun.SearchSolution3;
                        }
                        else
                        {
                            SaveSolutionReference = Left.Solution;
                        }
                        Left = Middle;
                        Middle = Right;
                        NewCalcDeriv_LineSearch(alpha, SaveSolutionReference, ref Right, Calcfg);

                        //  Right.deriv > 0 is best case but also stop if value increases at Right even with unexpected negative derivative
                        bool TimetoBreak = (Right.deriv > 0) || (Right.value > Middle.value);
                        SALSAUtility.SynchronizeMPIvariable(ref TimetoBreak);

                        if (TimetoBreak)
                            break;
                    }
                    bool DerivSearch2 = DerivativeSearch(Left, Middle, Right, out DerivSearchMethod,
                                                         out BestFromLineSearch, out LocalIterations, Calcfg);
                    LineSearchMethod = 120 + DerivSearchMethod;
                    SALSAUtility.SynchronizeMPIvariable(ref DerivSearch2);
                    success = DerivSearch2;
                }
            } // End case LineFactor > 1.0

            //  Copy Best Solution if needed
            LineIterations += LocalIterations;

            if (!success)
            {
                MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
                GlobalParameterSet = true;
                return false;
            }

            // Success -- Change Current Solution to Best from Line Search
            LineFactor = BestFromLineSearch.alpha;
            ExtraDecrease = CurrentSolution.Chisquared - BestFromLineSearch.value;
            SALSABLAS.LinearCombineVector(BestFromLineSearch.Solution.xshift, LineFactor,
                                          Hotsun.EndingLinePositionSolution.xshift, 0.0,
                                          Hotsun.EndingLinePositionSolution.xshift);
            CopySolution(CurrentSolution, BestFromLineSearch.Solution);
            MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
            GlobalParameterSet = true;
            Hotsun.pred1 *= LineFactor;
            Hotsun.pred2 *= LineFactor*LineFactor;
            Hotsun.pred3 *= LineFactor*LineFactor;

            return true;
        }

        // End Line Search

        public static bool DerivativeSearch(SearchStuff Left, SearchStuff Middle, SearchStuff Right,
                                            out int SearchMethod, out SearchStuff Best, out int Iterations,
                                            Hotsun.CalcfgSignature Calcfg)
        {
            // Do a derivative Secant search but never stray outside interval by requiring End point derivatives opposite sign
            //  Revert to Value search if problems


            Iterations = 0;
            Best = Middle;
            bool WrongSignofRightDeriv = Right.deriv <= 0;
            SALSAUtility.SynchronizeMPIvariable(ref WrongSignofRightDeriv);

            if (WrongSignofRightDeriv)
            {
                // Inconsistent Derivatives -- Use value Search
                int GoldenSectionMethod = 0;
                bool GoldenResult = GoldenSectionValueSearch(Left, Middle, Right, out GoldenSectionMethod, out Best,
                                                             out Iterations, Calcfg);
                SearchMethod = 20 + GoldenSectionMethod;
                SALSAUtility.SynchronizeMPIvariable(ref GoldenResult);
                return GoldenResult;
            }

            int SimpleSearchMethod = 0;
            int LocalIterations = 0;

            bool UnexpectedSignofMiddleDeriv = Middle.deriv > 0;
            SALSAUtility.SynchronizeMPIvariable(ref UnexpectedSignofMiddleDeriv);

            if (UnexpectedSignofMiddleDeriv)
            {
                // Best solution betwenn Left and Middle// First see if Deriv search possible. That needs Left to be defined with a negative derivative
                bool noleftsolution = (Left.Solution == null);
                bool derivsearchforbidden = noleftsolution;
                if (!derivsearchforbidden)
                {
                    if (Left.deriv > 0)
                        derivsearchforbidden = true;
                }
                SALSAUtility.SynchronizeMPIvariable(ref derivsearchforbidden);

                if (derivsearchforbidden)
                {
                    int GoldenSectionMethod = 0;
                    bool GoldenResult = GoldenSectionValueSearch(Left, Middle, Right, out GoldenSectionMethod, out Best,
                                                                 out Iterations, Calcfg);
                    SearchMethod = 30 + GoldenSectionMethod;
                    SALSAUtility.SynchronizeMPIvariable(ref GoldenResult);
                    return GoldenResult;
                }
                if (Middle.value > Right.value)
                {
                    bool result2 = SimpleDerivativeSearch(Left, Right, out SimpleSearchMethod, out Best,
                                                          out LocalIterations, Calcfg, Middle);
                    Iterations += LocalIterations;
                    SearchMethod = SimpleSearchMethod + 5;
                    return result2;
                }
                bool result1 = SimpleDerivativeSearch(Left, Middle, out SimpleSearchMethod, out Best,
                                                      out LocalIterations, Calcfg, Right);
                Iterations += LocalIterations;
                SearchMethod = SimpleSearchMethod + 10;
                return result1;
            }
            bool result = SimpleDerivativeSearch(Middle, Right, out SimpleSearchMethod, out Best, out LocalIterations,
                                                 Calcfg, Left);
            Iterations += LocalIterations;
            SearchMethod = SimpleSearchMethod;
            return result;
        }

        // End DerivativeSearch

        public static bool SimpleDerivativeSearch(SearchStuff OneSide, SearchStuff OtherSide, out int SearchMethod,
                                                  out SearchStuff Best, out int Iterations,
                                                  Hotsun.CalcfgSignature Calcfg, SearchStuff NotUsed)
        {
            double x2 = OneSide.alpha;
            double x3 = OtherSide.alpha;
            double InitialRange = x3 - x2;
            double TargetDelta = 0.01*(OtherSide.value - OneSide.value);
            double TargetValueChange = 0.01*Math.Abs(OtherSide.value - OneSide.value);
            SearchMethod = 5;
            Iterations = 0;

            while (Iterations < 10)
            {
                ++Iterations;
                double x4 = x3 - (x3 - x2)*OtherSide.deriv/(OtherSide.deriv - OneSide.deriv);
                var NewPoint = new SearchStuff();
                Desertwind NewSolution = matchingSolution(NotUsed.Solution, OneSide.Solution, OtherSide.Solution);
                NewCalcDeriv_LineSearch(x4, NewSolution, ref NewPoint, Calcfg);
                double oldtestvalue = Math.Min(OneSide.value, OtherSide.value);

                bool NewPointValueLarger = (NewPoint.value > OneSide.value);
                SALSAUtility.SynchronizeMPIvariable(ref NewPointValueLarger);

                if (NewPointValueLarger)
                {
                    if (NewPoint.deriv < 0)
                    {
                        Best = OneSide;
                        SearchMethod = 1;
                        return true;
                    }
                    OtherSide = NewPoint;
                }
                else
                {
                    if (NewPoint.deriv > 0)
                    {
                        Best = NewPoint;
                        SearchMethod = 2;
                        return true;
                    }
                    OneSide = NewPoint;
                }
                x3 = OtherSide.alpha;
                x2 = OneSide.alpha;

                bool casetobreak1 = ((OtherSide.value - OneSide.value) < TargetDelta) || ((x3 - x2) < 0.01*InitialRange);
                SALSAUtility.SynchronizeMPIvariable(ref casetobreak1);

                if (casetobreak1)
                {
                    SearchMethod = 3;
                    break;
                }

                double newtestvalue = Math.Min(OneSide.value, OtherSide.value);
                bool casetobreak2 = (newtestvalue > oldtestvalue - TargetValueChange);
                SALSAUtility.SynchronizeMPIvariable(ref casetobreak2);

                if (casetobreak2)
                {
                    SearchMethod = 4;
                    break;
                }
            }
            Best = OneSide;
            return true;
        }


        public static bool GoldenSectionValueSearch(SearchStuff Left, SearchStuff Middle, SearchStuff Right,
                                                    out int GoldenSectionMethod, out SearchStuff Best,
                                                    out int Iterations, Hotsun.CalcfgSignature Calcfg)
        {
            // Do a value search to find minimum assuming three values are set
            //  Return End point if this is minimum (Derivative search will not do this if End Point has positive derivative

            double goldennumber = 1.618033989;
            GoldenSectionMethod = 0;
            Iterations = 0;
            Best = Middle;

            if (Left.Solution == null)
            {
                GoldenSectionMethod = 1;
                return false;
            }

            //  First Check for Consistency of data
            bool leftlow = (Left.value <= Middle.value);
            SALSAUtility.SynchronizeMPIvariable(ref leftlow);

            if (leftlow)
            {
                bool leftbiggerthanright = (Left.value > Right.value);
                SALSAUtility.SynchronizeMPIvariable(ref leftbiggerthanright);

                if (leftbiggerthanright)
                {
                    Best = Right;
                    GoldenSectionMethod = 2;
                    return true;
                }
                else
                {
                    Best = Left;
                    GoldenSectionMethod = 3;
                    return false;
                }
            }

            bool rightlow = (Right.value <= Middle.value);
            SALSAUtility.SynchronizeMPIvariable(ref rightlow);

            if (rightlow)
            {
                GoldenSectionMethod = 4;
                Best = Right;
                return true;
            }

            double x1 = Left.alpha;
            double x2 = Middle.alpha;
            double x3 = Right.alpha;
            bool inconsistentpositions = (((x2 - x1) <= 0.0) || ((x3 - x2) <= 0.0));
            SALSAUtility.SynchronizeMPIvariable(ref inconsistentpositions);

            if (inconsistentpositions)
            {
                // Data Error in positions
                GoldenSectionMethod = 5;
                return false;
            }
            double InitialRange = x3 - x1;
            double TargetDelta = 0.01*Math.Max(Left.value - Middle.value, Right.value - Middle.value);
            double LastChange = Math.Min(Left.value - Middle.value, Right.value - Middle.value);

            //  Data Consistent
            GoldenSectionMethod = 8;
            while (Iterations <= 10)
            {
                var NewPoint = new SearchStuff();
                Desertwind NewSolution = matchingSolution(Left.Solution, Middle.Solution, Right.Solution);
                double x4;
                SearchStuff save;

                bool leftintervallarger = ((x2 - x1) > (x3 - x2));
                SALSAUtility.SynchronizeMPIvariable(ref leftintervallarger);

                if (leftintervallarger)
                {
                    // Left Interval Bigger -- place New Point to left
                    x4 = x2 - (x3 - x2)/goldennumber;
                    NewCalcDeriv_LineSearch(x4, NewSolution, ref NewPoint, Calcfg);

                    bool NewPointsmallerthanMiddle = (NewPoint.value < Middle.value);
                    SALSAUtility.SynchronizeMPIvariable(ref NewPointsmallerthanMiddle);

                    if (NewPointsmallerthanMiddle)
                    {
                        // Its x1 x4 x2
                        save = Middle;
                        Middle = NewPoint;
                        Right = save;
                    }
                    else
                    {
                        // Its x4 x2 x3
                        Left = NewPoint;
                    }
                }
                else
                {
                    // // Right Interval Bigger -- place New Point to right
                    x4 = x2 + (x2 - x1)/goldennumber;
                    NewCalcDeriv_LineSearch(x4, NewSolution, ref NewPoint, Calcfg);

                    bool NewPointsmallerthanMiddle = (NewPoint.value < Middle.value);
                    SALSAUtility.SynchronizeMPIvariable(ref NewPointsmallerthanMiddle);

                    if (NewPointsmallerthanMiddle)
                    {
                        // Its x2 x4 x3
                        save = Middle;
                        Middle = NewPoint;
                        Left = save;
                    }
                    else
                    {
                        // Its x1 x2 x4
                        Right = NewPoint;
                    }
                }
                ++Iterations;
                x1 = Left.alpha;
                x2 = Middle.alpha;
                x3 = Right.alpha;

                bool casetobreakonsectionsize = ((x3 - x1) < 0.01*InitialRange);
                SALSAUtility.SynchronizeMPIvariable(ref casetobreakonsectionsize);

                if (casetobreakonsectionsize)
                {
                    GoldenSectionMethod = 6;
                    break;
                }


                double CurrentChange = Math.Min(Left.value - Middle.value, Right.value - Middle.value);
                bool casetobreakonvaluechange = ((CurrentChange <= TargetDelta) && (LastChange <= TargetDelta));
                SALSAUtility.SynchronizeMPIvariable(ref casetobreakonvaluechange);

                if (casetobreakonvaluechange)
                {
                    GoldenSectionMethod = 7;
                    break;
                }
                LastChange = CurrentChange;
            }
            Best = Middle;


            return true;
        }

        // End GoldenSectionValueSearch

        // Find available Solution from set of pre allocated ones
        public static Desertwind matchingSolution(Desertwind a, Desertwind b, Desertwind c)
        {
            var Possibles = new[]
                                {
                                    Hotsun.SearchSolution1, Hotsun.SearchSolution2, Hotsun.SearchSolution3,
                                    Hotsun.SearchSolution4
                                };
            var used = new int[4];

            for (int i = 0; i < 4; i++)
            {
                used[i] = 0;
            }
            used[findsolutionindex(a, Possibles)] = 1;
            used[findsolutionindex(b, Possibles)] = 1;
            used[findsolutionindex(c, Possibles)] = 1;

            for (int i = 0; i < 4; i++)
            {
                if (used[i] == 0)
                    return Possibles[i];
            }
            Exception e = SALSAUtility.SALSAError("Invalid Solution Set");

            throw (e);
        }

        // End matchingSolution

        //  Find which of pre allocated Solutions correspondds to a particular one stored in trial
        public static int findsolutionindex(Desertwind trial, Desertwind[] Possibles)
        {
            int arraysize = Possibles.GetLength(0);

            for (int i = 0; i < arraysize; i++)
            {
                if (trial.Equals(Possibles[i]))
                    return i;
            }

            Exception e = SALSAUtility.SALSAError("Invalid Solution Index");

            throw (e);
        }

        // End findsolutionindex

        public static void ExistingCalcDeriv_LineSearch(double alpha, Desertwind SearchSolution,
                                                        ref SearchStuff SearchPoint)
        {
            // Calculate value and derivative at a point where Calcfg has been called

            SearchPoint.alpha = alpha;
            SearchPoint.value = SearchSolution.Chisquared;
            SearchPoint.Solution = SearchSolution;
            SearchPoint.deriv = -2.0*
                                SALSABLAS.VectorScalarProduct(SearchPoint.Solution.first,
                                                              Hotsun.EndingLinePositionSolution.xshift);
            return;
        }

        // End CalcDeriv_LineSearch at a point where Calcfg has been called

        public static void NewCalcDeriv_LineSearch(double alpha, Desertwind NewSolution, ref SearchStuff SearchPoint,
                                                   Hotsun.CalcfgSignature Calcfg)
        {
            // Calculate value and derivative at a point where Calcfg has NOT been called

            SearchPoint.Solution = NewSolution;

            //  Set up parameters for next call to Calcfg
            //  This adds in estimated change to param and zeros first derivative, Second Derivative and Chisq (zerocr)
            SALSABLAS.LinearCombineVector(NewSolution.param, 1.0, Hotsun.BeginningLinePositionSolution.param, -alpha,
                                          Hotsun.EndingLinePositionSolution.xshift);
            MakeVectorGlobal(NewSolution.param, Hotsun.GlobalParameter);
            GlobalParameterSet = true;

            Hotsun.idata = 0;
            ZeroSolution(NewSolution);
            //  Call Calcfg to calculate Taylor expansion
            SALSAUtility.StartSubTimer(2);
            violat = Calcfg(NewSolution);
            SALSAUtility.StopSubTimer(2);

            SearchPoint.alpha = alpha;
            SearchPoint.value = NewSolution.Chisquared;
            SearchPoint.deriv = -2.0*
                                SALSABLAS.VectorScalarProduct(NewSolution.first,
                                                              Hotsun.EndingLinePositionSolution.xshift);
            return;
        }

        // End CalcDeriv_LineSearch at a point where Calcfg has NOT been called

        // End Little struct SearchStuff

        public static void CopySolution(Desertwind Solution2, Desertwind Solution1)
        {
            int Numberparms = Hotsun.Number_VectorParameters;

            if (Hotsun.DecomposeParameters)
            {
                Numberparms = SALSAUtility.PointCount_Process;
            }
            SALSABLAS.CopyVector(Solution2.param, Solution1.param, 0, Numberparms);
            SALSABLAS.CopyVector(Solution2.first, Solution1.first, 0, Numberparms);
            SALSABLAS.CopyVector(Solution2.xshift, Solution1.xshift, 0, Numberparms);

            if (Hotsun.fullmatrixset)
            {
                SALSABLAS.CopyMatrix(Solution2.FullMatrix, Solution1.FullMatrix);
                SALSABLAS.CopyMatrix(Solution2.ExactFullMatrix, Solution1.ExactFullMatrix);
            }
            else
            {
                SALSABLAS.CopyVector(Solution2.DiagonalofMatrix, Solution1.DiagonalofMatrix, 0, Numberparms);
                SALSABLAS.CopyVector(Solution2.ExactDiagonalofMatrix, Solution1.ExactDiagonalofMatrix, 0, Numberparms);
            }

            Solution2.Chisquared = Solution1.Chisquared;
            Solution2.IterationCalculated = Solution1.IterationCalculated;
        }

        public static void SaveBestSolution(Desertwind CurrentSolution, Desertwind BestSolution)
        {
            CopySolution(BestSolution, CurrentSolution);
            Hotsun.zeromn = Hotsun.zerocr;
            Hotsun.Qbest = Hotsun.Q;
            return;
        }

        // End SaveBestSolution()

        public static void RestoreBestSolution(Desertwind CurrentSolution, Desertwind BestSolution)
        {
            CopySolution(CurrentSolution, BestSolution);
            Hotsun.zerocr = Hotsun.zeromn;

            MakeVectorGlobal(Hotsun.CurrentSolution.param, Hotsun.GlobalParameter);
            GlobalParameterSet = true; // Set Indicator  that Global Parameters are set
        }

        // End RestoreBestSolution()

        public static void ZeroSolution(Desertwind CurrentSolution)
        {
            //  xshift and param are NOT zeroed
            SALSABLAS.zrword(CurrentSolution.first);

            if (Hotsun.fullmatrixset)
            {
                SALSABLAS.zrword(CurrentSolution.FullMatrix);
                SALSABLAS.zrword(CurrentSolution.ExactFullMatrix);
            }
            else
            {
                SALSABLAS.zrword(CurrentSolution.DiagonalofMatrix);
                SALSABLAS.zrword(CurrentSolution.ExactDiagonalofMatrix);
            }
            CurrentSolution.Chisquared = 0.0;
            CurrentSolution.IterationCalculated = Hotsun.numit + 1;
            Hotsun.zerocr = 0.0;
        }

        #region Nested type: SearchStuff

        public struct SearchStuff
        {
            public Desertwind Solution; // Full Solution Detail -- xshift NOT set typically
            public double alpha; // Free Parameter in Line Search
            public double deriv; // First Derivative wrt alpha
            public double value; // Chisq Value
        }

        #endregion

        // End ZeroSolution()
    }

    // End Class ManxcatCentral
}

// End namespace Manxcat