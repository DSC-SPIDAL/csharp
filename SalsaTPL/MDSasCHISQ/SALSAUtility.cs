#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using HPC.Utilities;
using MPI;
using Manxcat;
using Salsa.Core.Configuration.Sections;
using Environment = MPI.Environment;
using TDistance = System.Int16;

#else
using TDistance = System.Double;
#endif

namespace SALSALibrary
{
    public class SALSAUtility
    {
        public static ParallelOptions _parallelOptions;

        public static int checkerboard = 0;
                          // If 0 store full matrix, If 1 lower triangular, If 2 use load balanced checkerboard pattern

        public static Environment MPI_Environment; // MPI Environment
        public static Intracommunicator MPI_communicator = null; //MPI communicator
        public static int MPI_Rank = 0; // Rank of process
        public static int MPI_Size = 1; // Number of MPI Processes
        public static string ParallelPattern = ""; // Basic Parallel Pattern
        public static string PatternLabel = ""; // One line Label for Output

        public static int DistanceProcessingOption = 0;
                          // Control Processing of Distance on Input if 0 or 1, leave distance unchanged; if =2 square input distance

        public static bool CalcFixedCrossFixed = true;
                           // If true calculate fixed-fixed terms in Chisq. Must be false if StoredDistanceOption =3;

        public static int chisqcomponent = -1;
                          // -1 ignore, =1 calculate Varied-Varied only, =2 Calculate Varied Fixed, =3 Calculate Fixed Varied, =4 Calculate Fixed Fixed

        public static int StoredDistanceOption = 2;
                          // Specify how distance data stored in memory =1 Illegal, =2 Used by Used (Used = Varied+Fixed), =3 Varied by Used

        public static int DiskDistanceOption = 2;
                          // Specify whats in distance data on disk =1 Original by Original, =2 Used by Used (Used = Varied+Fixed), =3 Varied by Used

        public static double DistanceCut = 0.96;
                             // Ignore if negative. If positive undefine all distances larger than this

        public static int LinkCut = 5;
                          // Delete all points with this # Links <= this (Must have at least 3 for possibly well determined problem)

        public static double AllowedDeletedFraction = 0.25;
                             // Veto run if Deleted Points (from missing distances) greater than this fraction

        public static int TransformMethod = 0;
                          // Define transformation Method if nonzero; 11 Special Blast centered at one

        public static double TransformParameter = 0.125;
                             // Define Parameter to control transformation. It is power of (1-d) for TransformMethod 11

        public static double UndefinedDistanceValue = -1.0;
                             // If positive replace undefined distances by this value; use -1.0 if want to unset these distances

        public static double[] DistanceWeightCuts;
                               // List of Distance Cuts for Weights. These define upper limits of distance bins with "infinity" as upper and 0 of course as lower

        public static int NumberDistanceWeightCuts = 0;
                          // Number of Distance Weight Cuts. This is one less than number of bins

        public static double[] ActualWeightCuts;
                               // List of Weights averaging to 1. Current alkgorithm makes each distance bin have equal total weight so bins with lowest number of points have largest weight

        public static int PointCount_Global = 0;
                          // Total number of used points summed over all threads and processes; same as Hotsun.Number_Vectors

        // This includes varied and fixed points but not ignored points-- it is = NumberFixedPoints + NumberVariedPoints
        // It can include deleted points if some removed 
        public static int PointCount_Process = 0; // Total number of points summed over all threads in this process
        public static int PointCount_Largest = 0; // Largest number of points in all processes
        public static int PointStart_Process = 0; //	First data point in this process
        public static int VariedPointStart_Process = 0; //  First Varied Point (on global count) in Process      
        public static int VariedPointCount_Process = 0; //  Number of Varied points in process

        public static int NumberFixedPoints = 0; // Number of Fixed Points
        public static int NumberVariedPoints = 0; // Number of Varied Points
        public static int NumberOriginalPoints = 0; // Total Number of Points specified in Input
        public static int SALSASHIFT = 100; // Allow room for flags in OriginalPointDisposition -- Should be 1 or higher

        public static int[] OriginalPointDisposition;
                            // = 0 Original Point Not Used; = SALSASHIFT + n This is Varied Point n; = -SALSASHIFT -m This is Fixed Point m

        public static int[] FixedPointOriginal; // Value of Original Point Index Corresponding to this fixed point
        public static int[] VariedPointOriginal; // Value of Original Point Index Corresponding to this varied point

        public static int[] UsedPointtoOriginalPointMap;
                            // Value of Original Point Index Corresponding to this Used point

        public static int[] OriginalPointtoUsedPointMap;
                            // Value of Used Point Index Corresponding to this Original point; -1 if OriginalPoint not used

        public static int[] ActualtoNaiveUsedOrder;
                            // Value of original used position of reordered (for load balancing) used points

        public static int[] NaivetoActualUsedOrder; // Value of actual used position for original used position

        public static bool Usedreordered = false; // If true the used are reordered

        public static SALSADataPointProperties[] GlobalPointProperties; // Properties of USED Point System
        public static SALSAFileProperties GlobalFileProperties; // Properties of Current USED data set

        public static string[] OriginalDataLabels; // Labels of Current ORIGINAL data set
        public static SALSADataPointProperties[] OriginalPointProperties; // Properties of ORIGINAL Point System
        public static SALSAFileProperties OriginalFileProperties; // Properties of Original data set

        //  Within a job data points will be divided into MPI_Size parts -- each part is assigned to a separate MPI Process
        public static int[] PointsperProcess = null; //how many data points each process will take care
        public static int[][] PointsperThreadperProcess = null; // Number of data points in each process-thread 

        //	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
        public static int[] PointsperThread = null; //how many data points each thread will take care
        public static int[] StartPointperThread = null; //the starting point that a thread will take care

        public static int ThreadCount = 1; // maximum number of parallel threads in a process
        public static int NodeCount = 1; // maximum number of separate nodes in run
        public static int MPIperNodeCount = 1; // Number of MPI processes per node
        public static int MPIIOStrategy = 0; // MPI I/O Strategy

        public static bool sequentialBLAS = false;
                           // If true calculate BLAS sequentially using imput lengths; if false use SALSAParallelism and do using Threads/MPI

        public static TDistance[][] PointDistances; // Point Distances
        public static int MatrixBreakFactor = 1;
        public static int[] ArrayDivision;
        public static int[] Indexsubtraction;
        public static int DivisionSize;

        public static DateTime startTime;
        public static DateTime endTime;

        // Todo. revert to HiPerfTimer if necessary
//        public static HiPerfTimer PreciseTimer;   //    Hold Precise Timing
//        public static HiPerfTimer[] SubTimers;   // Timing Objects
        public static HiPerfTimer PreciseTimer;
        public static HiPerfTimer[] SubTimers;

        public static int NumberofSubTimings = 0; // Number of subtimings
        public static double HPDuration = 0.0; // Time with Precision
        public static double[] SubDurations; // Hold partial timing
        public static string[] SubTimingNames; //  Labels of partial timing
        public static bool[] SubTimingEnable;
        public static int MPIREDUCETiming = -1;
        public static int MPIREDUCETiming1 = -1;
        public static int MPISENDRECEIVETiming = -1;
        public static int MPIBROADCASTTiming = -1;
        public static int ThreadTiming = -1;

        /* Parameters for density graph creation */
        public static double Xmaxbound = 1.5; // bounding max x value
        public static double Ymaxbound = 1.5; // bounding max y value
        public static int Xres = 50; // resolution of x axis 
        public static int Yres = 50; // resolution of y axis
        public static double Alpha = 2;
        public static double Pcutf = 0.85;
        public static bool Normalize = true;
        public static string ClusterFile = string.Empty;
        public static HashSet<int> SelectedClusters; // Selected clusters to plot density graph
        public static bool IsClustersSelected = false;

        /* Parameters for html creation */
        public static string ManxcatRunName = string.Empty;
        public static string ManxcatRunDescription = string.Empty;


        //  These are general parameters for C# codes
        public static ArrayList CosmicOutput = new ArrayList(1000); // Monitoring Output
        public static bool ConsoleDebugOutput = true; // If true send Monitoring output to console
        public static int DebugPrintOption = 2; // Control Printing (= 0 None, ==1 Summary, = 2 Full)

        //  Set up Parallel Threading

        public static ParallelOptions ParallelOptions
        {
            get { return _parallelOptions; }
        }

        public static void SetupParallelOptions()
        {
            _parallelOptions = new ParallelOptions();
            _parallelOptions.MaxDegreeOfParallelism = ThreadCount;
        }


        public static void SetupDistanceWeights()
        {
            // Set up Distance Weights
            if (!string.IsNullOrEmpty(ManxcatCentral.Configuration.DistanceWeightsCuts))
            {
                var sep = new[] {','};
                string[] weightCuts = ManxcatCentral.Configuration.DistanceWeightsCuts.Trim().Split(sep);
                NumberDistanceWeightCuts = weightCuts.Length;
                DistanceWeightCuts = new double[NumberDistanceWeightCuts];
                ActualWeightCuts = new double[NumberDistanceWeightCuts + 1];
                for (int i = 0; i < NumberDistanceWeightCuts; i++)
                {
                    DistanceWeightCuts[i] = double.Parse(weightCuts[i]);
                }
            }
        }

        public static Exception SALSAError(string message)
        {
            Console.WriteLine("SALSA Error " + message);
            var e = new Exception(message);
            return e;
        }

        // end SALSAError

        // PrintOption = 0 Essential Printout
        // PrintOption = 1 Summary Printout
        // PrintOption = 2 Only if full print out requested
        public static void SALSAPrint(int PrintOption, string StufftoPrint)
        {
            if (MPI_Rank != 0)
                return;
            if (DebugPrintOption < PrintOption)
                return;
            CosmicOutput.Add(StufftoPrint);

            if (ConsoleDebugOutput)
                Console.WriteLine(StufftoPrint);
            return;
        }

        // End SALSAPrint

        public static void SALSAGracefulend(string directory, out int readinstruction)
        {
            readinstruction = -1;
            if (MPI_Rank != 0)
                return;

            string filename = directory + "\\GracefulEnd.txt";
            if (!File.Exists(filename))
                return;

            try
            {
                using (StreamReader sr = File.OpenText(filename))
                {
                    // Read first line of file
                    String inputLineStr;
                    if ((inputLineStr = sr.ReadLine()) == null)
                    {
                        sr.Close();
                        return;
                    }
                    if (inputLineStr.Length > 0)
                    {
                        inputLineStr = inputLineStr.Trim(new[] {' ', '\t'});
                        string[] split = inputLineStr.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                        if ((split.Length > 0) && (split[0].Length > 0))
                        {
                            readinstruction = Convert.ToInt32(split[0]);
                            SALSAPrint(0, " Termination Instruction " + readinstruction.ToString());
                        }
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAError(" Failure reading " + filename + " " + e);
                throw (e);
            }
        }

        public static void SALSAStatus(string directory, string Message)
        {
            if (MPI_Rank != 0)
                return;
            string statusfilename = directory + "\\Status.txt";

            try
            {
                StreamWriter sw = null;

                if (!string.IsNullOrEmpty(statusfilename))
                {
                    sw = new StreamWriter(statusfilename, false, Encoding.UTF8);
                }

                if (sw != null)
                {
                    sw.WriteLine(Message);
                }

                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing status data" + e);
                throw (e);
            }
        }

        // End SALSAStatus

        // todo: saliya - Possible Improvement: handle comments in this function to use the List instead of the string
        //                you will have to search for the places where comment is manipulated in other
        //                places as well
        public static void SALSAUpdateMetaData(ManxcatSection Configuration)
        {
            string Comment = Configuration.Comment;
            for (int i = 0; i < CosmicOutput.Count; i++)
            {
                if (i == 0)
                    continue;
                Comment += "\n" + CosmicOutput[i];
            }
            Comment = Comment.Replace(':', '-');
            Configuration.Comment = Comment;
        }

        // end SALSAUpdateMetaData(string Comment)

        public static void InitializeTiming(int InputNumberofTimers)
        {
            NumberofSubTimings = InputNumberofTimers;
            SubTimers = new HiPerfTimer[NumberofSubTimings]; // Timing Objects

            SubDurations = new double[NumberofSubTimings]; // Hold partial timing
            SubTimingEnable = new bool[NumberofSubTimings];
            SubTimingNames = new string[NumberofSubTimings];

            for (int itimer = 0; itimer < NumberofSubTimings; itimer++)
            {
                SubTimers[itimer] = new HiPerfTimer();
                SubDurations[itimer] = 0.0;
                SubTimingEnable[itimer] = true;
            }

            PreciseTimer = new HiPerfTimer();
            PreciseTimer.Start();

            startTime = DateTime.Now; // Set the initial time.
        }

        // End InitializeTiming

        public static void SetUpSubTimer(int TimingIndex, string TimingLabel)
        {
            if (TimingIndex >= NumberofSubTimings)
            {
                SALSAError("Error in Timing Index " + TimingIndex.ToString() + " Max " + NumberofSubTimings.ToString());
                return;
            }
            SubTimingNames[TimingIndex] = TimingLabel;
            SubTimingEnable[TimingIndex] = true;
            SubDurations[TimingIndex] = 0.0;
        }

        // End SetUpSubTimer

        public static void SetUpMPISubTimers(int StartTimingIndex, string MPISetLabel)
        {
            if ((StartTimingIndex + 4) >= NumberofSubTimings)
            {
                SALSAError("Error in  MPI Timing Index " + (StartTimingIndex + 2).ToString() + " Max " +
                           NumberofSubTimings.ToString());
                return;
            }
            SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce");
            SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI SR");
            SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI Bcast");
            SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Global Reductions");
            SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "Thread Global Reductions");

            if ((MPISetLabel != "") && (MPISetLabel != "Lib "))
                return;
            MPIREDUCETiming = StartTimingIndex;
            MPISENDRECEIVETiming = StartTimingIndex + 1;
            MPIBROADCASTTiming = StartTimingIndex + 2;
            MPIREDUCETiming1 = StartTimingIndex + 3;
            ThreadTiming = StartTimingIndex + 4;
        }

        // End SetUpMPISubTimers

        public static void InterimTiming()
        {
            PreciseTimer.Stop();
            // Todo. revert to HiPerfTimer if necessary
//            HPDuration += PreciseTimer.Duration;
            HPDuration += PreciseTimer.Duration*0.001;
            PreciseTimer.Start();
        }

        // end EndTiming

        public static void EndTiming()
        {
            endTime = DateTime.Now; // Set the final time.
            PreciseTimer.Stop();
            HPDuration += PreciseTimer.Duration*0.001;
        }

        // end EndTiming

        public static void StartSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Start();
            }
        }

        // End StartSubTimer

        public static double StopSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return 0.0;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Stop();
                double tmp = SubTimers[TimingIndex].Duration*0.001;
                SubDurations[TimingIndex] += tmp;
                return tmp;
            }
            return 0.0;
        }

        // End StopSubTimer


        // Synchronize a boolean across MPI processes
        public static void SynchronizeMPIvariable(ref bool cosmicboolean)
        {
            if (MPI_Size > 1)
            {
                StartSubTimer(MPIBROADCASTTiming);
                MPI_communicator.Broadcast(ref cosmicboolean, 0);
                StopSubTimer(MPIBROADCASTTiming);
            }
            return;
        }

        // End synchronizeboolean(bool cosmicboolean)

        // Synchronize a Double across MPI processes
        public static void SynchronizeMPIvariable(ref double cosmicdouble)
        {
            if (MPI_Size > 1)
            {
                StartSubTimer(MPIBROADCASTTiming);
                MPI_communicator.Broadcast(ref cosmicdouble, 0);
                StopSubTimer(MPIBROADCASTTiming);
            }
            return;
        }

        // End synchronizeboolean(double cosmicdouble)
        public static void SynchronizeMPIvariable(ref int cosmicint)
        {
            if (MPI_Size > 1)
            {
                StartSubTimer(MPIBROADCASTTiming);
                MPI_communicator.Broadcast(ref cosmicint, 0);
                StopSubTimer(MPIBROADCASTTiming);
            }
            return;
        }

        // End synchronizeboolean(int cosmicint)

        public static MPIReducePlusIndex MinwithIndex(MPIReducePlusIndex one, MPIReducePlusIndex two)
        {
            if (one.index < 0)
                return two;
            if (two.index < 0)
                return one;
            if (one.value < two.value)
                return one;
            else
                return two;
        }

        public static MPIReducePlusIndex MaxwithIndex(MPIReducePlusIndex one, MPIReducePlusIndex two)
        {
            if (one.index < 0)
                return two;
            if (two.index < 0)
                return one;
            if (one.value > two.value)
                return one;
            else
                return two;
        }

        // Return Value and Index corresponding to minimum over all MPI Processes
        //  Replace Input values by minimum values
        public static void AllReduceMinWithIndex(ref double ProcessValue, ref int ProcessIndex)
        {
            if (MPI_Size > 1)
            {
                var LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = MPI_communicator.Allreduce(LocalStructure, MinwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        }

        // End AllReduceMinWithIndex

        // Return Value and Index corresponding to maximum over all MPI Processes
        //  Replace Input values by maximum values
        public static void AllReduceMaxWithIndex(ref double ProcessValue, ref int ProcessIndex)
        {
            if (MPI_Size > 1)
            {
                var LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = MPI_communicator.Allreduce(LocalStructure, MaxwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        }

        // End AllReduceMaxWithIndex

        // Write performance results of process 0 into a file
        public static void WriteTiming_Cluster(string fileName, string RunLabel, int RunNumber, string DataFileName,
                                               string ProcessorName)
        {
            try
            {
                // Print results
                StreamWriter sw = null;

                // New file flag
                Boolean newfile = false;

                if (!File.Exists(fileName))
                    newfile = true;

                if (!string.IsNullOrEmpty(fileName))
                {
                    sw = new StreamWriter(fileName, true, Encoding.UTF8); //Append

                    if (newfile)
                    {
                        // Formats must match those in following data write
                        sw.Write(
                            String.Format(
                                "{0,-15}{1,-15}{2,-8}{3,-12}{4,-8}{5,-12}{6,-12}{7,-16}{8,-12}{9,-25}{10,-25}{11,-25}",
                                "Duration(ms)", "HPDuration(ms)", "Thread#", "MPIperNode", "Node#", "Pt#/Process",
                                "Pt#/Global", "Run Label", "Run Number", "DataFileName", "CurrentTime", "ProcessorName"));

                        for (int subtimer = 0; subtimer < NumberofSubTimings; subtimer++)
                        {
                            if (!SubTimingEnable[subtimer])
                                continue;

                            sw.Write(String.Format("{0,-15}{1,-8}", SubTimingNames[subtimer], // Subset Time    Name
                                                   "Ratio"));
                        }
                        sw.WriteLine(" "); // Blank line
                        sw.WriteLine(" "); // Blank line
                    }
                }

                if (sw != null)
                {
                    TimeSpan duration = endTime - startTime;

                    // Format syntax (Variable Number , Spaces with - left justified and no dash right aligned
                    sw.Write(
                        String.Format(
                            "{0,-15}{1,-15}{2,-8}{3,-12}{4,-8}{5,-12}{6,-12}{7,-16}{8,-12}{9,-25}{10,-25}{11,-25}",
                            Math.Round(duration.TotalMilliseconds), // Total time
                            // Todo. revert to HiPerfTimer if necessary
//                    Math.Round(SALSAUtility.HPDuration * 0.001, 0),                                                                                                        //High performance timer
                            Math.Round(HPDuration, 0), //High performance timer
                            ThreadCount, // Thread# per process
                            MPIperNodeCount, // Process# per node
                            NodeCount, // Node# aquired
                            PointCount_Process, // Local points
                            PointCount_Global, // Global pointsn
                            RunLabel, // Label for Run
                            RunNumber, // Run Number
                            DataFileName, // Name of input data file
                            startTime.ToLocalTime(), // Current Time
                            ProcessorName // Processor Used
                            ));

                    for (int subtimer = 0; subtimer < NumberofSubTimings; subtimer++)
                    {
                        if (!SubTimingEnable[subtimer])
                            continue;

                        // Todo. revert to HiPerfTimer if necessary
//                        double SubTime = SALSAUtility.SubDurations[subtimer] * 0.001;
                        double SubTime = SubDurations[subtimer];

//                        double SubRatio = SubTime / (HPDuration * 0.001);
                        double SubRatio = SubTime/(HPDuration);
                        sw.Write(String.Format("{0,-15}{1,-8}", Math.Round(SubTime), // Subset Time
                                               SubRatio.ToString("F4")));
                    }
                    sw.WriteLine(" "); // Blank line
                }

                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing data" + e);
                throw (e);
            }
        }
    }

    // End class SalsaUtility

    [Serializable]
    public class MPIReducePlusIndex
    {
        public int index;
        public double value;

        public MPIReducePlusIndex(int _index, double _value)
        {
            index = _index;
            value = _value;
        }
    }

    // End MPIMinPlusIndex
}

// end Namespace SALSALibrary