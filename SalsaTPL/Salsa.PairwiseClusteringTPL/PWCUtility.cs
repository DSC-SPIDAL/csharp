using System;
using System.Collections;
using System.IO;
using System.Text;
using MPI;
using Salsa.Core.Blas;
using System.Diagnostics;
using Salsa.Core;

using Salsa.PairwiseClusteringTPL;

#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using TDistance = System.Int16;
#else
using TDistance = System.Double;
#endif

namespace SALSALibrary
{
	public class PWCUtility
	{
        public static int PointCount_Global = 0; // Total number of points summed over all threads and processes
        public static int PointCount_Process = 0;	// Total number of points summed over all threads in this process
        public static int PointCount_Largest = 0;	// Largest number of points in all processes
        public static int PointStart_Process = 0;	//	First data point in this process  
        
        // Parallel Parameters
        public static int MPI_Rank = 0;                          // Rank of process
        public static int MPI_Size = 1;                          // Number of MPI Processes
        public static MPI.Environment MPI_Environment;           // MPI Environment
        public static Intracommunicator MPI_communicator = null; //MPI communicator

        public static string ParallelPattern = " "; // Pattern of parallel execution (e.g. 8x1x8 indicates [threads/process][processes/node][nodes])
        public static string PatternLabel = " "; // Title line for print

        //  Within a job data points will be divided into MPI_Size parts -- each part is assigned to a separate MPI Process
        public static int[] PointsperProcess = null;     //how many data points each process will take care
        public static int[][] PointsperThreadperProcess = null; // Number of data points in each process-thread

        //	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
        public static int[] PointsperThread = null;//how many data points each thread will take care
        public static int[] StartPointperThread = null;//the starting point that a thread will take care

        public static int ThreadCount = 1;  // maximum number of parallel threads in a process
        public static int NodeCount = 1;  // maximum number of separate nodes in run
        public static int MPIperNodeCount = 1;  // Number of MPI processes per node

        public static int checkerboard = 0;	// If 0 store full matrix, If 1 lower triangular, If 2 use load balanced checkerboard pattern
        public static int MPIIOStrategy = 0;  // MPI I/O Strategy

        public static Matrix<TDistance> PointDistances;

        // Timing Parameters
        public static DateTime startTime;
        public static DateTime endTime;
        public static Stopwatch PreciseTimer;   //    Hold Precise Timing
        public static int NumberofSubTimings = 0; // Number of subtimings
        public static Stopwatch[] SubTimers;   // Timing Objects
	    public static long TimerFrequency; // Frequency of stopwatches - ticks/sec
        public static double HPDuration = 0.0;    // Time with Precision
        public static double[] SubDurations;     // Hold partial timing
        public static string[] SubTimingNames;   //  Labels of partial timing
        public static bool[] SubTimingEnable;
        public static int MPIREDUCETiming = -1;
        public static int MPIREDUCETiming1 = -1;
        public static int MPISENDRECEIVETiming = -1;
        public static int MPISENDRECEIVEEigenTiming = -1;
        public static int MPIBROADCASTTiming = -1;
        public static int ThreadTiming = -1;

        //  These are general parameters for C# codes
        public static ArrayList CosmicOutput = new ArrayList(1000); // Monitoring Output
        public static bool ConsoleDebugOutput = true;               // If true send Monitoring output to console
        public static int DebugPrintOption = 2;                     // Control Printing (= 0 None, ==1 Summary, = 2 Full)

        //  Center Finding Parameters
        public static int NumberofCenters = 8;  // Number of centers to be found with each method
        public static double[] BucketFractions = new double[] { 0.15, 0.4, 0.75 };  // Fractions to be used in Bucket method -- centers are those with most neighbors in a radius determined by BucketFraction
        public static int NumberofBuckets = BucketFractions.Length; // Number of Buckets
        public static int addMDS = 1;   // Specify MDS versions of center finding; = 0 ignore Currently all nonzero values treated in same way
        public static string addMDSfile = "";   // File with MDS Information
        public static string Labelfile = "";    // File with Label and Length Information
        public static string ClusterNumberfile = "";    // File with Cluster Numbers; can be same as addMDSFile
	    public static int CenterPointsPerCenterTypeInOutput = 3; // number of center points to include for each center type in the output plot
	    public static string CenterPlotFile = ""; // output plot file with centers
        public static double MinimumDistanceCut = -0.1; // Insist Sequence distances bigger than this; negative values remove test
        public static int LinkCountinCenterFinding = 0;   // Minimum Number of links
        public static int LengthCut1 = -1;  // Insists first sequence has length bigger thanm this
        public static int LengthCut2 = -1;  // Insist comparison sequence has length bigger than this


	    public static Exception SALSAError(string message)
        {
            Console.WriteLine("SALSA Error " + message);
            Exception e = new Exception(message);
            return e;

        } // end SALSAError

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

        } // End SALSAPrint

        public static string PrintFixedInteger(int data, int digits)
        {
            string returned = data.ToString();
            for (int charloop = returned.Length; charloop < digits; charloop++)
                returned = " " + returned;
            return returned;

        }   // End PrintFixedInteger

        public static string LeftPrintFixedInteger(int data, int digits)
        {
            string returned = data.ToString();
            for (int charloop = returned.Length; charloop < digits; charloop++)
                returned = returned + " ";
            return returned;

        }   // End LeftPrintFixedInteger

        public static void InitializeTiming(int InputNumberofTimers)
        {
            NumberofSubTimings = InputNumberofTimers;
            SubTimers = new Stopwatch[NumberofSubTimings]; // Timing Objects
            SubDurations = new double[NumberofSubTimings];   // Hold partial timing
            SubTimingEnable = new bool[NumberofSubTimings];
            SubTimingNames = new string[NumberofSubTimings];

            for (int itimer = 0; itimer < NumberofSubTimings; itimer++)
            {
                SubTimers[itimer] = new Stopwatch();
                SubDurations[itimer] = 0.0;
                SubTimingEnable[itimer] = true;
            }
            PreciseTimer = Stopwatch.StartNew();

            startTime = DateTime.Now; // Set the initial time.
            TimerFrequency = Stopwatch.Frequency;

        }                             

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
        } // End SetUpSubTimer

        public static void SetUpMPISubTimers(int StartTimingIndex, string MPISetLabel)
        {
            if ((StartTimingIndex + 5) >= NumberofSubTimings)
            {
                SALSAError("Error in  MPI Timing Index " + (StartTimingIndex + 2).ToString() + " Max " + NumberofSubTimings.ToString());
                return;
            }
            SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce");
            SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI SendRec");
            SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI SendRec Eig");
            SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Bcast");
            SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "MPI Global Reductions");
            SetUpSubTimer(StartTimingIndex + 5, MPISetLabel + "Thread Global Reductions");

            if ((MPISetLabel != "") && (MPISetLabel != "Lib "))
                return;
            MPIREDUCETiming = StartTimingIndex;
            MPISENDRECEIVETiming = StartTimingIndex + 1;
            MPISENDRECEIVEEigenTiming = StartTimingIndex + 2;
            MPIBROADCASTTiming = StartTimingIndex + 3;
            MPIREDUCETiming1 = StartTimingIndex + 4;
            ThreadTiming = StartTimingIndex + 5;

        } // End SetUpMPISubTimers

        public static void InterimTiming()
        {
            PreciseTimer.Stop();
            HPDuration += (((double)(PreciseTimer.ElapsedTicks)) / TimerFrequency) * 1000000;
            PreciseTimer.Restart();

        } // end InterimTiming

        public static void EndTiming()
        {
            endTime = DateTime.Now; // Set the final time.
            PreciseTimer.Stop();
            HPDuration += (((double)(PreciseTimer.ElapsedTicks)) / TimerFrequency) * 1000000;

        } // end EndTiming

        public static void StartSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Start();
            }

        } // End StartSubTimer

        public static void StopSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Stop();
                SubDurations[TimingIndex] += (((double)(SubTimers[TimingIndex].ElapsedTicks)) / TimerFrequency) * 1000000;
                SubTimers[TimingIndex].Reset();
            }

        } // End StopSubTimer


        // Synchronize a boolean across MPI processes
        public static void SynchronizeMPIvariable(ref bool cosmicboolean)
        {
            if (PWCUtility.MPI_Size > 1)
            {
                PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                PWCUtility.MPI_communicator.Broadcast<bool>(ref cosmicboolean, 0);
                PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
            }
            return;
        } // End synchronizeboolean(bool cosmicboolean)

        // Synchronize a Double across MPI processes
        public static void SynchronizeMPIvariable(ref double cosmicdouble)
        {
            if (PWCUtility.MPI_Size > 1)
            {
                PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                PWCUtility.MPI_communicator.Broadcast<double>(ref cosmicdouble, 0);
                PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
            }
            return;
        } // End synchronizeboolean(double cosmicdouble)

        public static void SynchronizeMPIvariable(ref int cosmicint)
        {
            if (PWCUtility.MPI_Size > 1)
            {
                PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                PWCUtility.MPI_communicator.Broadcast<int>(ref cosmicint, 0);
                PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
            }
            return;
        } // End synchronizeboolean(int cosmicint)


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
            if (PWCUtility.MPI_Size > 1)
            {
                MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = PWCUtility.MPI_communicator.Allreduce<MPIReducePlusIndex>(LocalStructure, MinwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        } // End AllReduceMinWithIndex

        // Return Value and Index corresponding to maximum over all MPI Processes
        //  Replace Input values by maximum values
        public static void AllReduceMaxWithIndex(ref double ProcessValue, ref int ProcessIndex)
        {
            if (PWCUtility.MPI_Size > 1)
            {
                MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = PWCUtility.MPI_communicator.Allreduce<MPIReducePlusIndex>(LocalStructure, MaxwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        } // End AllReduceMaxWithIndex

        public static void WriteResults_Cluster(string ControlFileName, ArrayList linestooutput)
        {
            try
            {
                StreamWriter FileToWrite = null;

                if (!string.IsNullOrEmpty(ControlFileName))
                {
                    FileToWrite = new StreamWriter(ControlFileName, true, Encoding.UTF8); // Overwrite
                }


                if (FileToWrite != null)
                {
                    FileToWrite.WriteLine();

                    int LineLength = linestooutput.Count;

                    for (int countlines = 0; countlines < LineLength; countlines++)
                    {
                        FileToWrite.WriteLine(linestooutput[countlines]);
                    }
                }
                FileToWrite.Flush();
                FileToWrite.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing data {0} " + e, ControlFileName);
                throw (e);
            }
        }   // End WriteResults_Cluster

	}   // End PWCUtility

    [Serializable()]
    public class MPIReducePlusIndex
    {
        public int index;
        public double value;

        public MPIReducePlusIndex(int _index, double _value)
        {
            index = _index;
            value = _value;
        }
    } // End MPIMinPlusIndex

}   // End Namespace SALSALibrary
