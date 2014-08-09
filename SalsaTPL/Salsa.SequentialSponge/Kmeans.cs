using System;
using System.Collections;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using MPI;
using Salsa.Core;
using SALSALibrary;

namespace Salsa.DAVectorSponge
{

    public class Kmeans
    {
        public static ParallelOptions _parallelOptions; // Options for TPL Threading

        public static double[][] PointPosition; // Point Positions
        public static int[] InitialPointAssignment; // Initial Assignment of Points to Clusters

        public static int MaxNcent_Global = 50; // Maximum number of centers
        public static int ParameterVectorDimension = 2;        // Vector Dimension of Points
        public static int Ncent_Global = 1;        // Total Number of Clusters
        public static int Ncent_Global_Parallel = 50;    // Use Parallelism if total number of centers greater than this

        public static double[][] ClusterCenter;    // Cluster Centers stored in every node
        public static int[] ClusterSize;  // Size of Cluster (number of points)
        public static double[] ClusterRadius;  // Radius of Cluster
        public static double[] ClusterWidth;  // Width of Cluster

        //  Center Parallelism over threads and processes
        public static int CenterCount_Process = 0;	// Total number of Centers summed over all threads calculated in this process
        public static int CenterCount_Largest = 0;	// Largest number of points in all processes
        public static int CenterStart_Process = 0;	//	First data point in this process 

        //	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
        public static int[] FullParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
        public static int[] FullParallel_StartCenterperThread = null; //the starting point that a thread will take care for thread+process parallelism

        //  Center Parallelism over threads withon each process 
        public static int[] LocalParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
        public static int[] LocalParallel_StartCenterperThread = null;  //the starting point that a thread will take care for thread+process parallelism

        public static bool UseParallelismoverCenters = false;   // If true major center loops run in parallel over nodes

        public static int CountKmeansIterations = 0;  // Count Iterations in Kmeans
        public static int KmeansIterationCut = 1000;    // Cut on Iterations to stop
        public static double KmeansCenterChangeCut = 0.001; // Fractional change to stop
        public static double AverageCenterChange = 0.0; //  Average Center Change each iteration
        public static double AverageRadius = 0.0;   // Average Radius each iteration
        public static double AverageWidth = 0.0;   // Average Width each iteration

        public static void InitializeKmeans(double[][] PointPositionINPUT, ParallelOptions _parallelOptionsINPUT, string FileName, int ClusterPosition, int FirstClusterValue, int StartingPosition,
            int Ncent_GlobalINPUT, int MaxNcent_GlobalINPUT, int Ncent_Global_ParallelINPUT, int ParameterVectorDimensionINPUT, double CenterChangeINPUT, int IterationCutINPUT)
        {
            
            Kmeans.PointPosition = PointPositionINPUT;
            Kmeans._parallelOptions = _parallelOptionsINPUT;

            Kmeans.Ncent_Global = Ncent_GlobalINPUT;
            Kmeans.MaxNcent_Global = MaxNcent_GlobalINPUT;
            Kmeans.Ncent_Global_Parallel  = Ncent_Global_ParallelINPUT;
            Kmeans.KmeansCenterChangeCut = CenterChangeINPUT;
            Kmeans.KmeansIterationCut = IterationCutINPUT;
    
            Kmeans.ParameterVectorDimension = ParameterVectorDimensionINPUT;

            Kmeans.SetParallelCenterDecomposition();

            Kmeans.InitialPointAssignment = new int[DAVectorUtility.PointCount_Process];
            if (FileName.Length == 0)
                return;

            //  Read Initial Assignments
            DAVectorUtility.SALSAPrint(0, "Kmeans Read File " + FileName + " Points " + DAVectorUtility.PointCount_Global.ToString() + " Starting at position "
                + StartingPosition.ToString() + " Dimension " + Kmeans.ParameterVectorDimension.ToString() + " Cluster Position " + ClusterPosition.ToString() + " Initial value " + FirstClusterValue.ToString());
            Kmeans.ReadDataFromFile(FileName, ClusterPosition, FirstClusterValue,  StartingPosition);

        }   // End InitializeKmeans

        public static void SetupKmeans(int[] FullAssignment)
        {   // Set up initial assignment from existing array

            Parallel.For(0, Kmeans._parallelOptions.MaxDegreeOfParallelism, Kmeans._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    Kmeans.InitialPointAssignment[alpha] = FullAssignment[alpha + DAVectorUtility.PointStart_Process];

                }   // End loop over Points
            }); // End Sum over Threads

        }   // End SetupKmeans


        public static void RunKmeans(double[][] ClusterCenterINPUT, int[] ClusterSizeINPUT, double[] ClusteRadiusINPUT,
            out int Ncent_GlobalFINAL, out double AverageWidthFINAL)
        {

            ArrayList KeepPCfractions = new ArrayList(200);
            ArrayList KeepCCfractions = new ArrayList(200);

            //  Inherit Solution arrays
            Kmeans.ClusterCenter = ClusterCenterINPUT;
            Kmeans.ClusterSize = ClusterSizeINPUT;
            Kmeans.ClusterRadius = ClusteRadiusINPUT;
            Kmeans.ClusterWidth = new double[Kmeans.MaxNcent_Global];
            
            //  Set up TriangleInequality
            KmeansTriangleInequality.SetExternalFunctions(GetClusterRadius, GetClusterCenters, FindClusterCenters);
            KmeansTriangleInequality.InitializeTriangleInequality(Kmeans.PointPosition, Kmeans._parallelOptions, Kmeans.ClusterCenter,
                Kmeans.Ncent_Global, Kmeans.MaxNcent_Global, Kmeans.Ncent_Global_Parallel, Kmeans.ParameterVectorDimension);

            DAVectorUtility.SALSAPrint(0, "Start Kmeans ****** Number of Centers " + Kmeans.Ncent_Global.ToString() + " Max Number of Centers " + Kmeans.MaxNcent_Global.ToString()
                + " Center Limit for Parallelism " + Kmeans.Ncent_Global_Parallel.ToString() + " Vector Dimension " + Kmeans.ParameterVectorDimension.ToString());

            Kmeans.FindClusterCenters(true, Kmeans.InitialPointAssignment, null, null);
            Kmeans.CountKmeansIterations = 0;
            bool StartStop = false;
            int CountStops = 0;
            while (Kmeans.CountKmeansIterations < Kmeans.KmeansIterationCut)
            {
                double save1 = KmeansTriangleInequality.NumberFullDistancesCalculatedCC;
                double save2 = KmeansTriangleInequality.NumberFullDistancesCalculatedPC;

                KmeansTriangleInequality.NextIteration();
                ++Kmeans.CountKmeansIterations;
                bool WillStop = false;
                if (!StartStop)
                {
                    if (Kmeans.AverageCenterChange < Kmeans.AverageRadius * Kmeans.KmeansCenterChangeCut)
                        StartStop = true;
                }
                else
                {
                    ++CountStops;
                    if (CountStops > 10)
                        WillStop = true;
                }
                double tmp1 = (KmeansTriangleInequality.NumberFullDistancesCalculatedCC - save1) / (double)Kmeans.MaxNcent_Global;
                double tmp2 = (KmeansTriangleInequality.NumberFullDistancesCalculatedPC - save2) / ((double)Kmeans.MaxNcent_Global * (double)DAVectorUtility.PointCount_Global);
                double tmp3 = KmeansTriangleInequality.NumberFullDistancesCalculatedPC / ((double)Kmeans.MaxNcent_Global * (double)(DAVectorUtility.PointCount_Global * Kmeans.CountKmeansIterations));
                double tmp4 = (KmeansTriangleInequality.NumberFullDistancesCalculatedPC + KmeansTriangleInequality.NumberFullDistancesCalculatedCC) / ((double)Kmeans.MaxNcent_Global * (double)(DAVectorUtility.PointCount_Global * Kmeans.CountKmeansIterations));
                DAVectorUtility.SALSAPrint(0, "Iteration " + Kmeans.CountKmeansIterations.ToString() + " Average Center Change " + Kmeans.AverageCenterChange.ToString("E4")
                    + " Average Radius " + Kmeans.AverageRadius.ToString("E4") + " Average Width " + Kmeans.AverageWidth.ToString("E4")
                    + " CC calcs per C " + tmp1.ToString("F4") + " PC calcs per P&C " + tmp2.ToString("F6") 
                    + " Cumul PC / Max " + tmp3.ToString("F6") + " Cumul PC+CC / PC Max " + tmp4.ToString("F6"));
                KeepPCfractions.Add(tmp2);
                KeepCCfractions.Add(tmp1 / DAVectorUtility.PointCount_Global);
                if (((Kmeans.CountKmeansIterations % 10) == 1) || WillStop)
                {
                    string message = " Sizes";
                    for (int CenterIndex = 0; CenterIndex < Kmeans.Ncent_Global; CenterIndex++)
                        message += " " + Kmeans.ClusterSize[CenterIndex].ToString();
                    DAVectorUtility.SALSAPrint(0, message);
                }
                if(WillStop)
                    break;
            }
            DAVectorUtility.SALSAPrint(0, "End Kmeans Iterations " + Kmeans.CountKmeansIterations.ToString() + " Iteration Cut " + Kmeans.KmeansIterationCut.ToString() +
                " Average Center Change " + Kmeans.AverageCenterChange.ToString("E4") + " Average Radius " + Kmeans.AverageRadius.ToString("E4") +
                " Average Width " + Kmeans.AverageWidth.ToString("E4") + " Fractional Cut " + Kmeans.KmeansCenterChangeCut.ToString("F4"));
            KmeansTriangleInequality.PrintDiagnostics();
            string messagePC = "\nPC Calcs per Point iteration";
            string messageCC = "\nCC Calcs per Point iteration";
            int numPC = KeepPCfractions.Count;
            for (int linecount = 0; linecount < numPC; linecount++)
            {
                messagePC += " " + ((double) KeepPCfractions[linecount]).ToString("F4") + ",";
                messageCC += " " + ((double) KeepCCfractions[linecount]).ToString("F4") + ",";
            }
            DAVectorUtility.SALSAPrint(0, messagePC);
            DAVectorUtility.SALSAPrint(0, messageCC);
            Ncent_GlobalFINAL = Kmeans.Ncent_Global;
            AverageWidthFINAL = Kmeans.AverageWidth;

            //  Print Histograms
            if (KmeansTriangleInequality.UseTriangleInequality != 0)
            {
                KmeansTriangleInequality.PlotPointHistograms(Math.Sqrt(AverageWidthFINAL));
                KmeansTriangleInequality.PlotCenterHistograms(Math.Sqrt(AverageWidthFINAL));
            }
            return;

        }   // End RunKmeans()

        public static void SetParallelCenterDecomposition()
        {
            Kmeans.UseParallelismoverCenters = false;
            if (Kmeans.Ncent_Global <= Kmeans.Ncent_Global_Parallel)
            {
                Kmeans.CenterStart_Process = 0;
                Kmeans.CenterCount_Process = Kmeans.MaxNcent_Global;
                Kmeans.CenterCount_Largest = Kmeans.MaxNcent_Global;
                return;
            }
            Kmeans.UseParallelismoverCenters = true;

            //	First divide centers among processes
            Range[] processRanges = RangePartitioner.Partition(Kmeans.MaxNcent_Global, DAVectorUtility.MPI_Size);
            Range processRange = processRanges[DAVectorUtility.MPI_Rank];  // The answer for this process

            Kmeans.CenterStart_Process = processRange.StartIndex;
            Kmeans.CenterCount_Process = processRange.Length;

            Kmeans.CenterCount_Largest = int.MinValue;
            foreach (Range r in processRanges)
                Kmeans.CenterCount_Largest = Math.Max(r.Length, Kmeans.CenterCount_Largest);


            //	Now divide centers among threads for this process
            Range[] Full_ThreadRanges = RangePartitioner.Partition(processRange, DAVectorUtility.ThreadCount);
            Kmeans.FullParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
            Kmeans.FullParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

            for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
            {
                Kmeans.FullParallel_CentersperThread[j] = Full_ThreadRanges[j].Length;
                Kmeans.FullParallel_StartCenterperThread[j] = Full_ThreadRanges[j].StartIndex;
            }

            //  Process only Center Parallelism
            Range[] Local_ThreadRanges = RangePartitioner.Partition(Kmeans.MaxNcent_Global, DAVectorUtility.ThreadCount);
            Kmeans.LocalParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
            Kmeans.LocalParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

            for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
            {
                Kmeans.LocalParallel_CentersperThread[j] = Local_ThreadRanges[j].Length;
                Kmeans.LocalParallel_StartCenterperThread[j] = Local_ThreadRanges[j].StartIndex;
            }

        }   // End SetParallelCenterDecomposition()

        public static void GetClusterCenters(int CenterIndex, double[] IndividualClusterCenter)
        {   //  Return Cluster Center vector

            IndividualClusterCenter = Kmeans.ClusterCenter[CenterIndex];
            return;

        }   // End GetClusterCenters

        public static double GetClusterRadius(int CenterIndex)
        {   // Return Cluster Radius

            return Kmeans.ClusterRadius[CenterIndex];

        }   // End GetClusterRadius

        // Use LastClusterCenter if Size 0
        public static void FindClusterCenters(bool begin, int[] NearestCentertoPoint, double[] Distance_NearestCentertoPoint, double[][] LastClusterCenter)
        {   // Calculate Cluster Parameters

            GlobalReductions.FindVectorDoubleSum3 FindCenterVectorSums = new GlobalReductions.FindVectorDoubleSum3(DAVectorUtility.ThreadCount,
                Kmeans.ParameterVectorDimension, Kmeans.Ncent_Global);
            GlobalReductions.FindVectorIntSum FindCenterSizeSums = new GlobalReductions.FindVectorIntSum(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);
            GlobalReductions.FindVectorDoubleMax FindClusterRadius = new GlobalReductions.FindVectorDoubleMax(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);
            GlobalReductions.FindVectorDoubleSum FindClusterWidth = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, Kmeans.Ncent_Global);
;
            Parallel.For(0, Kmeans._parallelOptions.MaxDegreeOfParallelism, Kmeans._parallelOptions, (ThreadIndex) =>
            {
                FindCenterVectorSums.startthread(ThreadIndex);
                FindCenterSizeSums.startthread(ThreadIndex);
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int ClusterIndex = NearestCentertoPoint[alpha];
                    if ((ClusterIndex >= Kmeans.Ncent_Global) || (ClusterIndex < 0))
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal Cluster Index " + ClusterIndex.ToString() + " Number " + Kmeans.Ncent_Global.ToString()
                            + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString() );
                        throw (e);
                    }
                    FindCenterVectorSums.addapoint(ThreadIndex, PointPosition[alpha], ClusterIndex);
                    FindCenterSizeSums.addapoint(ThreadIndex, 1, ClusterIndex);
                    if (!begin)
                    {
                        FindClusterRadius.addapoint(ThreadIndex, Distance_NearestCentertoPoint[alpha], ClusterIndex);
                        FindClusterWidth.addapoint(ThreadIndex, Distance_NearestCentertoPoint[alpha] * Distance_NearestCentertoPoint[alpha], ClusterIndex);
                    }

                }   // End loop over Points
            }); // End Sum over Threads

            FindCenterVectorSums.sumoverthreadsandmpi();
            FindCenterSizeSums.sumoverthreadsandmpi();
            if (!begin)
            {
                FindClusterRadius.sumoverthreadsandmpi();
                FindClusterWidth.sumoverthreadsandmpi();
            }
            Kmeans.AverageRadius = 0.0;
            Kmeans.AverageWidth = 0.0;
            Kmeans.AverageCenterChange = 0.0;

            if (Kmeans.UseParallelismoverCenters)
            {   // Centers Parallel over Threads NOT nodes
                double[] AccumulateRadius = new double[DAVectorUtility.ThreadCount];
                double[] AccumulateWidth = new double[DAVectorUtility.ThreadCount];
                double[] AccumulateCenterChange = new double[DAVectorUtility.ThreadCount];
                for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                {
                    AccumulateRadius[ThreadIndex] = 0.0;
                    AccumulateWidth[ThreadIndex] = 0.0;
                    AccumulateCenterChange[ThreadIndex] = 0.0;
                }

                Parallel.For(0, Kmeans._parallelOptions.MaxDegreeOfParallelism, Kmeans._parallelOptions, (ThreadIndex) =>
                {
                    int indexlen = KmeansTriangleInequality.LocalParallel_CentersperThread[ThreadIndex];
                    int beginpoint = KmeansTriangleInequality.LocalParallel_StartCenterperThread[ThreadIndex];
                    for (int CenterIndex = beginpoint; CenterIndex < indexlen + beginpoint; CenterIndex++)
                    {
                        Kmeans.ClusterSize[CenterIndex] = FindCenterSizeSums.TotalVectorSum[CenterIndex];

                        if (Kmeans.ClusterSize[CenterIndex] > 0)
                        {
                            double divisor = 1.0 / (double)Kmeans.ClusterSize[CenterIndex];
                            for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++)
                            {
                                int totalindex = VectorIndex + CenterIndex * Kmeans.ParameterVectorDimension;
                                Kmeans.ClusterCenter[CenterIndex][VectorIndex] = FindCenterVectorSums.TotalVectorSum[totalindex] * divisor;
                            }

                            if (!begin)
                            {
                                Kmeans.ClusterRadius[CenterIndex] = FindClusterRadius.TotalVectorMax[CenterIndex];
                                Kmeans.ClusterWidth[CenterIndex] = FindClusterWidth.TotalVectorSum[CenterIndex];
                                AccumulateRadius[ThreadIndex] += Kmeans.ClusterRadius[CenterIndex];
                                AccumulateWidth[ThreadIndex] += Kmeans.ClusterWidth[CenterIndex];
                                AccumulateCenterChange[ThreadIndex] += DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(Kmeans.ClusterCenter[CenterIndex], LastClusterCenter[CenterIndex]);
                            }
                        }
                        else
                        {   // Zero Size Cluster
                            if(begin)
                            {
                                Exception e = DAVectorUtility.SALSAError("Empty Input Cluster " + CenterIndex.ToString() + " Number " + Kmeans.Ncent_Global.ToString()
                                    + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                                throw (e);
                            }
                            for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++)
                                Kmeans.ClusterCenter[CenterIndex][VectorIndex] = LastClusterCenter[CenterIndex][VectorIndex];
                            Kmeans.ClusterRadius[CenterIndex] = 0.0;
                            Kmeans.ClusterWidth[CenterIndex] = 0.0;
                        }

                    }
                }); // End Sum over Threads
                if (!begin)
                {
                    for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                    {
                        Kmeans.AverageRadius += AccumulateRadius[ThreadIndex];
                        Kmeans.AverageWidth += AccumulateWidth[ThreadIndex];
                        Kmeans.AverageCenterChange += AccumulateCenterChange[ThreadIndex];
                    }
                }
            }

            else
            {   // Centers Sequential

                for (int CenterIndex = 0; CenterIndex < Kmeans.Ncent_Global; CenterIndex++)
                {
                    Kmeans.ClusterSize[CenterIndex] = FindCenterSizeSums.TotalVectorSum[CenterIndex];

                    if (Kmeans.ClusterSize[CenterIndex] > 0)
                    {
                        double divisor = 1.0 / (double)Kmeans.ClusterSize[CenterIndex];
                        for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++)
                        {
                            int totalindex = VectorIndex + CenterIndex * Kmeans.ParameterVectorDimension;
                            Kmeans.ClusterCenter[CenterIndex][VectorIndex] = FindCenterVectorSums.TotalVectorSum[totalindex] * divisor;
                        }

                        if (!begin)
                        {
                            Kmeans.ClusterRadius[CenterIndex] = FindClusterRadius.TotalVectorMax[CenterIndex];
                            Kmeans.ClusterWidth[CenterIndex] = FindClusterWidth.TotalVectorSum[CenterIndex];
                            Kmeans.AverageRadius += Kmeans.ClusterRadius[CenterIndex];
                            Kmeans.AverageWidth += Kmeans.ClusterWidth[CenterIndex];
                            Kmeans.AverageCenterChange += DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(Kmeans.ClusterCenter[CenterIndex], LastClusterCenter[CenterIndex]);
                        }
                    }
                    else
                    {
                        if (begin)
                        {
                            Exception e = DAVectorUtility.SALSAError("Empty Input Cluster " + CenterIndex.ToString() + " Number " + Kmeans.Ncent_Global.ToString()
                                + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                            throw (e);
                        }
                        Kmeans.ClusterRadius[CenterIndex] = 0.0;
                        Kmeans.ClusterWidth[CenterIndex] = 0.0;
                        for (int VectorIndex = 0; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++)
                            Kmeans.ClusterCenter[CenterIndex][VectorIndex] = LastClusterCenter[CenterIndex][VectorIndex];
                    }

                }   // End Sequential Center Loop
            }

            if (begin)
                return;
            Kmeans.AverageCenterChange = Kmeans.AverageCenterChange / Kmeans.Ncent_Global;
            Kmeans.AverageRadius = Kmeans.AverageRadius / Kmeans.Ncent_Global;
            Kmeans.AverageWidth = Kmeans.AverageWidth / DAVectorUtility.PointCount_Global;
            return;

        }   // End FindClusterCenters(int[] NearestCentertoPoint, double[][] LastClusterCenter)

        public static void ReadDataFromFile(string fname, int ClusterPosition, int FirstClustervalue, int StartPointPosition)
        {
            char[] _sep = new[] { ' ', ',', '\t' };

            int FirstPointPosition = 0;
            int TotalNumberPointstoRead = 0;
            FirstPointPosition = DAVectorUtility.PointStart_Process;
            TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;
            Random RandomObject = new Random(10101010 + DAVectorUtility.MPI_Rank);
            if(ClusterPosition < 0)
                DAVectorUtility.SALSAPrint(0, "Random Start 10101010 plus rank ******************* Option " + ClusterPosition.ToString());
            int MinSplitSize = ClusterPosition+1;
            if (StartPointPosition >= 0)
                MinSplitSize = Math.Max(MinSplitSize, StartPointPosition + Kmeans.ParameterVectorDimension);
            else
            {
                Exception e = DAVectorUtility.SALSAError("Illegal Start Position on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank.ToString()
                    + " POsition " + StartPointPosition.ToString() + " Number to Read " + TotalNumberPointstoRead.ToString());
                throw (e);
            }
            bool success = false;
            string line = " Unset";
            int CountLinesinFile = 0;

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
                            string[] splits = line.Trim().Split(_sep, StringSplitOptions.RemoveEmptyEntries);
                            if (splits.Length < MinSplitSize)
                            {
                                DAVectorUtility.SALSAPrint(0, "Count " + CountLinesinFile.ToString() + " Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                continue;
                            }   // Skip header lines

                            double junk;
                            if (!Double.TryParse(splits[StartPointPosition], out junk))
                                continue;   // Skip header lines

                            if (CountLinesinFile < FirstPointPosition)
                            {
                                CountLinesinFile += 1;
                                continue;
                            }

                            int ActualPointPosition = CountLinesinFile - FirstPointPosition;
                            int label = 0;

                            Kmeans.PointPosition[ActualPointPosition][0] = double.Parse(splits[StartPointPosition]);
                            Kmeans.PointPosition[ActualPointPosition][1] = double.Parse(splits[StartPointPosition + 1]);
                            if (Kmeans.ParameterVectorDimension > 2)
                            {
                                for (int VectorIndex = 2; VectorIndex < Kmeans.ParameterVectorDimension; VectorIndex++)
                                    Kmeans.PointPosition[ActualPointPosition][VectorIndex] = double.Parse(splits[VectorIndex + StartPointPosition]);
                            }

                            if (ClusterPosition >= 0)
                            {
                                if (!Int32.TryParse(splits[ClusterPosition], out label))
                                    label = FirstClustervalue;
                                Kmeans.InitialPointAssignment[ActualPointPosition] = label - FirstClustervalue;
                            }
                            else
                            {
                                Kmeans.InitialPointAssignment[ActualPointPosition] = RandomObject.Next(Program.InitialNcent);
                                if (ClusterPosition == -2)
                                {   // Force each cluster to have one point
                                    if (CountLinesinFile < Program.InitialNcent)
                                        Kmeans.InitialPointAssignment[ActualPointPosition] = CountLinesinFile;
                                }
                                if (ClusterPosition == -3)
                                {
                                    int divisor = Program.NumberDataPoints / Program.InitialNcent;
                                    if (CountLinesinFile % divisor == 0)
                                        Kmeans.InitialPointAssignment[ActualPointPosition] = CountLinesinFile / divisor;
                                }
                                if (ClusterPosition == -4)
                                {
                                    int divisor = Program.NumberDataPoints / Program.InitialNcent;
                                    Kmeans.InitialPointAssignment[ActualPointPosition] = CountLinesinFile / divisor;
                                }
                            }
                            ++ActualPointPosition;
                            ++CountLinesinFile;
                            if (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead) )
                                break;

                        }
                    }
                    if ( CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead) )
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank.ToString()
                            + " Lines in File " + CountLinesinFile.ToString() + " Number to Read " + TotalNumberPointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Points data " + DAVectorUtility.MPI_Rank.ToString() + " " + CountLinesinFile.ToString() + " Start "
                    + FirstPointPosition.ToString() + " Number " + TotalNumberPointstoRead.ToString() + " " + line + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File read error " + fname);
                throw (e);
            }
        }   // End ReadDataFromFile  

    }   //  End Kmeans

}   // End Namespace Salsa.DAVectorSponge
