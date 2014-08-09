using System;
using System.Collections;
using System.IO;
using System.Text;
using HPC.Utilities;
using MPI;
using Salsa.Core.Blas;
using Salsa.Core;


namespace SALSALibrary
{
	public class DAVectorUtility
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

        // Timing Parameters
        public static DateTime startTime;
        public static DateTime endTime;
        public static HiPerfTimer PreciseTimer;   //    Hold Precise Timing
        public static int NumberofSubTimings = 0; // Number of subtimings
        public static HiPerfTimer[] SubTimers;   // Timing Objects
        public static double HPDuration = 0.0;    // Time with Precision
        public static double[] SubDurations;     // Hold partial timing
        public static string[] SubTimingNames;   //  Labels of partial timing
        public static int[] SubTimingCalls;     // Number of calls
        public static bool[] SubTimingEnable;
        public static int[] TimingOutputOrder;  // Order of Output
        public static int MPIREDUCETiming1 = -1;
        public static int MPIREDUCETiming2 = -1;
        public static int MPIREDUCETiming3 = -1;
        public static int MPIREDUCETiming4 = -1;
        public static int MPIREDUCETiming5 = -1;
        public static int MPIREDUCETiming6 = -1;
        public static int MPISENDRECEIVETiming = -1;
        public static int MPIGATHERTiming = -1;
        public static int MPIDistributedREDUCETiming = -1;
        public static int MPIBROADCASTTiming = -1;
        public static int MPISynchTiming = -1;
        public static int ThreadTiming = -1;    // viewed as MPI timing

        //  These are general parameters for C# codes
        public static ArrayList CosmicOutput = new ArrayList(1000); // Monitoring Output
        public static bool ConsoleDebugOutput = true;               // If true send Monitoring output to console
        public static int DebugPrintOption = 2;                     // Control Printing (= 0 None, ==1 Summary, = 2 Full)

        public static ArrayList TemperatureValues = new ArrayList(10000);   // Record Temperatures per step
        public static ArrayList ClusterCountValues = new ArrayList(10000);  // Record Cluster Counts per Step

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

        public static void SALSAFullPrint(int PrintOption, string StufftoPrint)
        {
            if (DebugPrintOption < PrintOption)
                return;
            if (MPI_Rank == 0)
                CosmicOutput.Add(StufftoPrint);

            if (ConsoleDebugOutput)
                Console.WriteLine(" Node:" + MPI_Rank.ToString() + " " + StufftoPrint);
            return;

        } // End SALSAFullPrint

        public static void SALSASyncPrint(int PrintOption,string GlobalStufftoPrint, string StufftoPrint)
        {
            if (DebugPrintOption < PrintOption)
                return;
            string TotalStufftoPrint = "";
            if(StufftoPrint.Length > 0)
                TotalStufftoPrint = " Node:" + MPI_Rank.ToString() + " " + StufftoPrint;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            TotalStufftoPrint = DAVectorUtility.MPI_communicator.Allreduce<string>(TotalStufftoPrint, Operation<string>.Add);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            TotalStufftoPrint = GlobalStufftoPrint + TotalStufftoPrint;
            if (MPI_Rank != 0)
                return;
            CosmicOutput.Add(TotalStufftoPrint);

            if (ConsoleDebugOutput)
                Console.WriteLine(TotalStufftoPrint);
            return;

        } // End SALSASyncPrint

        public static string PrintFixedInteger(int data, int digits)
        {
            string returned = data.ToString();
            for (int charloop = returned.Length; charloop < digits; charloop++)
                returned = " " + returned;
            return returned;

        }   // End PrintFixedInteger

        public static string PadString(string start, int digits)
        {
            string returned = start;
            for (int charloop = start.Length; charloop < digits; charloop++)
                returned += " ";
            return returned;
        }   // End padstring

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
            SubTimers = new HiPerfTimer[NumberofSubTimings]; // Timing Objects
            SubDurations = new double[NumberofSubTimings];   // Hold partial timing
            SubTimingEnable = new bool[NumberofSubTimings];
            SubTimingNames = new string[NumberofSubTimings];
            SubTimingCalls = new int[NumberofSubTimings];
            TimingOutputOrder = new int[NumberofSubTimings];

            for (int itimer = 0; itimer < NumberofSubTimings; itimer++)
            {
                SubTimers[itimer] = new HiPerfTimer();
                SubDurations[itimer] = 0.0;
                SubTimingEnable[itimer] = true;
                SubTimingCalls[itimer] = 0;
                TimingOutputOrder[itimer] = itimer;
            }
            PreciseTimer = new HiPerfTimer();
            PreciseTimer.Start();

            startTime = DateTime.Now; // Set the initial time.

        }                             // End InitializeTiming

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
            if ((StartTimingIndex + 12) >= NumberofSubTimings)
            {
                SALSAError("Error in  MPI Timing Index " + (StartTimingIndex + 12).ToString() + " Max " + NumberofSubTimings.ToString());
                return;
            }
            SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce1(GlobalReductions)");
            SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI Reduce1(GlobalReductions while Dist)");
            SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI Reduce2(DAVectorClustering)");
            SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Reduce3(DAVectorDist)");
            SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "MPI Reduce4(Dist Reduce Initial)");
            SetUpSubTimer(StartTimingIndex + 5, MPISetLabel + "MPI Reduce5(Dist Exec Control)");
            SetUpSubTimer(StartTimingIndex + 6, MPISetLabel + "MPI Reduce6(Dist Reduce Transport)");
            SetUpSubTimer(StartTimingIndex + 7, MPISetLabel + "MPI DistributedReduce");
            SetUpSubTimer(StartTimingIndex + 8, MPISetLabel + "MPI Gather");
            SetUpSubTimer(StartTimingIndex + 9, MPISetLabel + "MPI Bcast");
            SetUpSubTimer(StartTimingIndex + 10, MPISetLabel + "MPI SendRecv");
            SetUpSubTimer(StartTimingIndex + 11, MPISetLabel + "MPI Sync");
            SetUpSubTimer(StartTimingIndex + 12, MPISetLabel + "Thread Timing");

            if ((MPISetLabel != "") && (MPISetLabel != "Lib "))
                return;
            MPIREDUCETiming1 = StartTimingIndex;    // StartTimingIndex+1 also used for different calls in Dist
            MPIREDUCETiming2 = StartTimingIndex + 2;
            MPIREDUCETiming3 = StartTimingIndex + 3;
            MPIREDUCETiming4 = StartTimingIndex + 4;
            MPIREDUCETiming5 = StartTimingIndex + 5;
            MPIREDUCETiming6 = StartTimingIndex + 6;
            MPIBROADCASTTiming = StartTimingIndex + 9;
            MPIGATHERTiming = StartTimingIndex + 8;
            MPIDistributedREDUCETiming = StartTimingIndex + 7;
            MPISENDRECEIVETiming = StartTimingIndex + 10;
            MPISynchTiming = StartTimingIndex + 11;
            ThreadTiming = StartTimingIndex + 12;

        } // End SetUpMPISubTimers

        public static int MPIaddviaGather(int thisnodevalue)
        {
            int[] Nodevalues = new int[DAVectorUtility.MPI_Size];
            Nodevalues = DAVectorUtility.MPI_communicator.Allgather<int>(thisnodevalue);
            int maxval = 0;
            for (int nodes = 0; nodes < DAVectorUtility.MPI_Size; nodes++)
                maxval += Nodevalues[nodes];
            return maxval;
        }


        public static void InterimTiming()
        {
            PreciseTimer.Stop();
            HPDuration += PreciseTimer.Duration;
            PreciseTimer.Start();

        } // end InterimTiming

        public static void EndTiming()
        {
            endTime = DateTime.Now; // Set the final time.
            PreciseTimer.Stop();
            HPDuration += PreciseTimer.Duration;

        } // end EndTiming

        public static void StartSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Start();
                ++SubTimingCalls[TimingIndex];
            }

        } // End StartSubTimer

        public static void StopSubTimer(int TimingIndex)
        {
            if (TimingIndex < 0)
                return;

            if (SubTimingEnable[TimingIndex])
            {
                SubTimers[TimingIndex].Stop();
                SubDurations[TimingIndex] += SubTimers[TimingIndex].Duration;
            }

        } // End StopSubTimer

        public static void CopyVector(int[] VectorC, int[] VectorA, int StartIndex, int TotalSize)
        {   // Can be parallelized
                for (int LongIndex = 0; LongIndex < TotalSize; LongIndex++)
                    VectorC[LongIndex + StartIndex] = VectorA[LongIndex];
                return;

        }

        // Synchronize a boolean across MPI processes
        public static void SynchronizeMPIvariable(ref bool cosmicboolean)
        {
            if (DAVectorUtility.MPI_Size > 1)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
                DAVectorUtility.MPI_communicator.Broadcast<bool>(ref cosmicboolean, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
            }
            return;
        } // End synchronizeboolean(bool cosmicboolean)

        // Synchronize a Double across MPI processes
        public static void SynchronizeMPIvariable(ref double cosmicdouble)
        {
            if (DAVectorUtility.MPI_Size > 1)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
                DAVectorUtility.MPI_communicator.Broadcast<double>(ref cosmicdouble, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
            }
            return;
        } // End synchronizeboolean(double cosmicdouble)

        public static void SynchronizeMPIvariable(ref int cosmicint)
        {
            if (DAVectorUtility.MPI_Size > 1)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPISynchTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref cosmicint, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPISynchTiming);
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
            if (DAVectorUtility.MPI_Size > 1)
            {
                MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = DAVectorUtility.MPI_communicator.Allreduce<MPIReducePlusIndex>(LocalStructure, MinwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        } // End AllReduceMinWithIndex

        // Return Value and Index corresponding to maximum over all MPI Processes
        //  Replace Input values by maximum values
        public static void AllReduceMaxWithIndex(ref double ProcessValue, ref int ProcessIndex)
        {
            if (DAVectorUtility.MPI_Size > 1)
            {
                MPIReducePlusIndex LocalStructure = new MPIReducePlusIndex(ProcessIndex, ProcessValue);
                MPIReducePlusIndex TotalStructure = DAVectorUtility.MPI_communicator.Allreduce<MPIReducePlusIndex>(LocalStructure, MaxwithIndex);
                ProcessValue = TotalStructure.value;
                ProcessIndex = TotalStructure.index;
            }
        } // End AllReduceMaxWithIndex

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

        public static void SetUpHistogramRange(int PointsinHistogram, ref double Histmin, ref double Histmax)
        {   // Choose good range to get rounded labels

            if (Histmax <= Histmin)
                return;
            if (Math.Abs(Histmin) < 0.000000001)
                Histmin = 0.0;
            else
            {
                double newminvalue = NewBinSize(Math.Abs(Histmin));
                if (Histmin < 0.0)
                    Histmin = -newminvalue;
                else
                    Histmin = newminvalue;
            }
            double binsize = (Histmax - Histmin) / PointsinHistogram;
            Histmax = Histmin + PointsinHistogram * NewBinSize(binsize);
            return;

        }   // End SetUpHistogramRange

        public static double NewBinSize(double oldBinSize)
        {   // Round Bin size up to a pretty value
            if (oldBinSize <= 0.000000001 || oldBinSize > 10000000000.0)
                return oldBinSize;

            double logvalue = Math.Log10(oldBinSize);
            int intlogvalue = Convert.ToInt32(Math.Floor(logvalue));
            double fudgepower = 1.0 - Convert.ToDouble(intlogvalue);
            double fudge = Math.Pow(10.0, fudgepower);
            double scaled = fudge * oldBinSize;
            scaled = Math.Min(scaled, 100.0);
            double Intversionofnewbinsize = Math.Ceiling(scaled) / fudge;
            //            SALSAUtility.SALSAPrint(1, "Hist " + oldBinSize.ToString("F4") + " " + intlogvalue + " " + Intversionofnewbinsize.ToString("F4") + " scaled " + scaled.ToString("F2"));
            return Intversionofnewbinsize;
        }

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

    }   // End DAVectorUtility

}   // End Namespace SALSALibrary
