using System;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using MPI;
using Salsa.Core;
using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    public class KmeansTriangleInequality
    {
        public delegate double ClusterRadiusSignature(int CenterIndex);
        public static ClusterRadiusSignature GetRadiusperCenter;

        public delegate void GetIndividualClusterCenterSignature(int CenterIndex, double[] CenterPositions );
        public static GetIndividualClusterCenterSignature GetIndividualCenterPosition;

        public delegate void FindClusterCentersSignature(bool begin, int[] NearestCentertoPoint, double[] Distance_NearestCentertoPoint, double[][] LastClusterCenter);
        public static FindClusterCentersSignature FindClusterCenters;

        public static ParallelOptions _parallelOptions; // Options for TPL Threading

        public static double[][] PointPosition;    // Distributed Point Positions

        public static int[] NearestCentertoPoint;   // Nearest Cluster Center to each Point
        public static double[] Distance_NearestCentertoPoint;   // Distance from Point to Nearest Cluster Center

        public static PointData[] LB_Old_alpha_ClusterPointer;  //  "Old" Lower Bound data structure
        public static PointData[] LB_Last_alpha_ClusterPointer; //  "Last" Lower Bound data structure
        public static PointData[] LB_Current_alpha_ClusterPointer;  //  "Current" Lower Bound data structure

        public static PointData[] LB_Buffer1_alpha_ClusterPointer;  //  Buffer1 for "Current" and "Last" Lower Bound data structure
        public static PointData[] LB_Buffer2_alpha_ClusterPointer;  //  Buffer2 for "Current" and "Last" Lower Bound data structure

        public static double[][] CenterOld;    // This is only used in TriangleInequality and is "old" list of centers
        public static double[][] CenterLast;    // This is only used in TriangleInequality and is "last" list of centers
        public static double[][] CenterCurrent;    // Defined and stored outside TriangleInequality. This is a Pointer
        public static bool CenterCurrentAlreadyset = true;  // If true cluster centers set automatically
        public static CenterData[] CurrentInterCenter; // Center and Center-Center Data stored in all nodes
        public static double[][] DistributedCenterDifference;   // Center-Center data stored for "old centers"
        public static bool StartingOut = true;  // True 

        public static int Ncent_Global = 1;        // Total Number of Clusters

        public static int ParameterVectorDimension = 2;        // Vector Dimension of Points
        public static int MaxNcent_Global = 50; // Maximum number of centers
        public static int Ncent_Global_Parallel = 50;    // Use Parallelism if total number of centers greater than this
        public static int OldCenterOption = 0;  // -1 Don't use, 0 Use with incremental update, >0 Refresh every OldCenterOption iterations
        public static bool DoBackwardFacingTests = true;

        public static int CenterCount_Process = 0;	// Total number of Centers summed over all threads calculated in this process
        public static int CenterCount_Largest = 0;	// Largest number of points in all processes
        public static int CenterStart_Process = 0;	//	First data point in this process 

        //	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
        public static int[] FullParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
        public static int[] FullParallel_StartCenterperThread = null; //the starting point that a thread will take care for thread+process parallelism

        public static int[] LocalParallel_CentersperThread = null; //how many data points each thread will take care for thread+process parallelism
        public static int[] LocalParallel_StartCenterperThread = null;  //the starting point that a thread will take care for thread+process parallelism

        public static bool UseParallelismoverCenters = false;   // If true major center loops run in parallel
        //  Even if centers calculated in parallel, they are stored globally per process

        //  Specify parameters for controlling use of triangle inequality
        public static int UseTriangleInequality = 0;    // if 0 do NOT use trianglew inequality; > 0 use it in a way specified by integer
        public static double TriangleInequality_Delta1_old = 0.1;  // Test for center change and old lower bounds (normalized by radius)
        public static double TriangleInequality_Delta1_current = 0.1;  // Test for Center change and current lower bounds (normalized by radius)
        public static int MaxClusterLBsperPoint = 50;   // Maximum number of Lower Bound values
        public static int MaxCentersperCenter = 49;   // Maximum number of Centers in Center Difference Array

        public static GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Center;  // Accumulate Diagnostic Sums in a parallel Environment
        public static int NumberDiagnosticSums_Center = 0; //   Number of Center Diagnostic Sums
        public static int SaveMPI_Size = 0; //  Save DAVectorUtility.MPI_Size
        public static GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Points;  // Accumulate Diagnostic Sums in a parallel Environment
        public static int NumberDiagnosticSums_Points = 0; //   Number of Diagnostic Sums

        //  Specify parameters for monitoring use of triangle inequality
        public static double NumberFullDistancesCalculatedCC = 0.0;   // Count Number of Actual Full Distance Computations for center-center cases
        public static double NumberFullDistancesCalculatedPC = 0.0;   // Count Number of Actual Full Distance Computations for point-center cases
        public static double CalculationCenterIterations = 0.0;  // Count Center times Iterations
        public static double CountToobigAShift = 0.0;  // Count Centers with large shift between iteration
        public static double AccumulateRadius = 0.0;    // Accumulate Radius
        public static double AccumulateDmax = 0.0;  // Accumulate Dmax
        public static double CountCCEntries = 0.0;  // Accumulate Number of Center-Center distances
        public static double AccumulateCCvalues = 0.0;  // Average value for Center-Center distances
        public static double AccumulateLastCenterchanges = 0.0;  // Average value for Last Center-Center changes
        public static double AccumulateOldCenterchanges = 0.0;  // Average value for Old Center-Center changes

        public static double CalculationPointIterations = 0.0;   // Accumulate Product of Points and Iterations
        public static double AccumulateOLD_LB = 0.0;   // OLD Lower Bound used successfully
        public static double AccumulateCURRENT_LB = 0.0;   // CURRENT Lower Bound used successfully
        public static double AccumulateCURRENT_LBLimit = 0.0;   // CURRENT Lower Bound Limit used successfully
        public static double AccumulateOLD_LBLimit = 0.0;   // OLD Lower Bound Limit used successfully
        public static double AccumulateOLD_LBFail = 0.0;   // OLD Lower Bound used unsuccessfully
        public static double AccumulateCURRENT_LBFail = 0.0;   // CURRENT Lower Bound used unsuccessfully
        public static double CountLBEntries = 0.0;  //  Number in Lower Bound Array
        public static double AccumulateBestScore = 0.0;  // Best Score
        public static double AccumulateBreakLimit = 0.0;  // Break Score - Best Score
        public static double AccumulateMAX_LB = 0.0;  // Maximum Lower Bound
        public static double AccumulateBreakReason7 = 0.0;  // Break reason 7 -- while ends
        public static double AccumulateBreakReason8 = 0.0;  // Break reason 8 -- Center Scan Dmax
        public static double AccumulateBreakReason9 = 0.0;  // Break reason 9 -- Center Scan Middle
        public static double AccumulateBreakReason10 = 0.0;  // Break reason 10 -- Center Scan Start
        public static double NumberFullDistancesCalculatedPC_Old = 0.0;   // Count Number of Actual Full Distance Computations for point-center cases to refresh Old Centers
        public static double NumberFullDistancesCalculatedPC_Last = 0.0;   // Count Number of Actual Full Distance Computations for point-center cases to refresh Last Centers
        public static double NumberFullDistancesCalculatedLoop = 0.0;   // Number of Full Distances calculated in main lowwer bound loop
        public static double PossibleNumberFullDistancesCalculated = 0.0;   // Possible Number of Full Distances calculated in main lower bound loop

        public static double NumberCentersinBestSet = 0.0;  // Number of centers examined in best loop Method = -2
        public static double NumberCentersinSimpleSet = 0.0;    // Number of centers examined in Method = -1 loop
        public static double NumberCentersinCenterbasedset = 0.0;   // Number of centers examined in method >=0 loop
        public static double AverageCentersLookedAt = 0.0;  // Average number of centers NOT looked at
        public static double NumberCentersRemovedBackwards = 0.0; // Centers Removed in Backwards test
        public static double NumberBackwardTests = 0.0; // Number of Backwards test
        public static double NumberForwardTests = 0.0; // Number of Forwards test

        public static int CalculationIterations = 0;    // Divide NumberFullDistancesCalculated to get Number per position and per iteration
        public static bool RecalculateOldCenters = false;

        public static int HistogramSize = 202;  // Size of final histograms including underflow and overflow

        public class PointData
        {
            public double[] DistceLowerBound;    // Lower Bound on Point Center Distance
            public int[] CenterIndices;  // Centers whose Point-Center distances stored in DistceLowerBound
            public int NumCenters; // Number of Centers in list
            public int AssociatedClusterIndex;  // Center associated with this point
            public double NearestDistance;  // Distance to center
            public double LimitforRest;     // Limit on Centers not in CenterIndices

            public PointData(int ArraySize)
            {
                DistceLowerBound = new double[ArraySize];
                CenterIndices = new int[ArraySize];
                NumCenters = 0;
                AssociatedClusterIndex = -1;
                NearestDistance = 0.0;
            }

        }   // End PointData for TriangleInequality

        public class CenterData
        {   // Only gathered for current Centers i.e. Clusters
            // Center positions available for current and last iteration

            public static int NumberPackingInts = 2;
            public static int NumberPackingDoubles = 4;
            public int RecalculateStatus; // If 1 Recalculate all distances for this center at start of each x loop
            public double ClusterRadius;    // Radius of Cluster
            public int NearbyClusters;  // Number of Nearby Clusters Stored
            public double Dmax; // Minimum of distances NOT in CennterDifference
            public double LastCenterChange; // Value of distance between current and last center position
            public double OldCenterChange;  // Value of distance between current and old center position
            public double[] CenterDifference;   // Distances between centers
            public int[] CenterDifferenceIndices; // Indices of Centers in CenterDifference

            public CenterData()
            {
                CenterDifference = new double[MaxCentersperCenter];
                CenterDifferenceIndices = new int[MaxCentersperCenter];
            }   // End constructor CenterData()

            public void serialize(int[] intarray, double[] doublearray, int intShift, int doubleShift)
            {
                intarray[intShift] = this.NearbyClusters;
                intarray[1 + intShift] = this.RecalculateStatus;
                doublearray[doubleShift] = this.ClusterRadius;
                doublearray[1 + doubleShift] = this.Dmax;
                doublearray[2 + doubleShift] = this.LastCenterChange;
                doublearray[3 + doubleShift] = this.OldCenterChange;

                for (int loopindex = 0; loopindex < this.NearbyClusters; loopindex++)
                {
                    intarray[NumberPackingInts + loopindex + intShift] = this.CenterDifferenceIndices[loopindex];
                    doublearray[NumberPackingDoubles + loopindex + doubleShift] = this.CenterDifference[loopindex];
                }

            }   // End serialize

            public void deserialize(int[] intarray, double[] doublearray, int intShift, int doubleShift)
            {
                this.NearbyClusters = intarray[intShift];
                this.RecalculateStatus = intarray[1 + intShift];
                this.ClusterRadius = doublearray[doubleShift];
                this.Dmax = doublearray[1 + doubleShift];
                this.LastCenterChange = doublearray[2 + doubleShift];
                this.OldCenterChange = doublearray[3 + doubleShift];

                for (int loopindex = 0; loopindex < this.NearbyClusters; loopindex++)
                {
                    this.CenterDifferenceIndices[loopindex] = intarray[NumberPackingInts + loopindex + intShift];
                    this.CenterDifference[loopindex] = doublearray[NumberPackingDoubles + loopindex + doubleShift];
                }

            }   // End deserialize

        }   // End CenterData for TriangleInequality

        public static void SetExternalFunctions(ClusterRadiusSignature GetRadiusperCenterINPUT, GetIndividualClusterCenterSignature GetCenterPositionsINPUT,
            FindClusterCentersSignature FindClusterCentersINPUT)
        {
            GetRadiusperCenter = GetRadiusperCenterINPUT;
            GetIndividualCenterPosition = GetCenterPositionsINPUT;
            FindClusterCenters = FindClusterCentersINPUT;

        }   // End SetExternalFunctions

        public static void SetTriangleInequalityParameters(int UseTriangleInequalityINPUT, int MaxClusterLBsperPointINPUT, int MaxCentersperCenterINPUT,
            double TriangleInequality_Delta1_oldINPUT, double TriangleInequality_Delta1_currentINPUT, int OldCenterOptionINPUT, bool DoBackwardFacingTestsINPUT)
        {
            KmeansTriangleInequality.UseTriangleInequality = UseTriangleInequalityINPUT;    // if 0 do NOT use trianglew inequality; > 0 use it in a way specified by integer
            KmeansTriangleInequality.MaxClusterLBsperPoint = MaxClusterLBsperPointINPUT;   // Maximum number of Lower Bound values
            KmeansTriangleInequality.MaxCentersperCenter = MaxCentersperCenterINPUT;   // Maximum number of Centers in Center Difference Array
            KmeansTriangleInequality.TriangleInequality_Delta1_old = TriangleInequality_Delta1_oldINPUT;  // Test for center change and old lower bounds (normalized by radius)
            KmeansTriangleInequality.TriangleInequality_Delta1_current = TriangleInequality_Delta1_currentINPUT;  // Test for Center change and current lower bounds (normalized by radius)
            KmeansTriangleInequality.OldCenterOption = OldCenterOptionINPUT;
            KmeansTriangleInequality.DoBackwardFacingTests = DoBackwardFacingTestsINPUT;

        }   // End SetTriangleInequalityParameters

        public static void InitializeTriangleInequality(double[][] PointPositionINPUT, ParallelOptions _parallelOptionsINPUT, double[][] CenterCurrentINPUT,
            int Ncent_GlobalINPUT, int MaxNcent_GlobalINPUT, int Ncent_Global_ParallelINPUT, int ParameterVectorDimensionINPUT)
        {
            KmeansTriangleInequality.PointPosition = PointPositionINPUT;
            KmeansTriangleInequality._parallelOptions = _parallelOptionsINPUT;
            KmeansTriangleInequality.SaveMPI_Size = DAVectorUtility.MPI_Size;

            KmeansTriangleInequality.CenterCurrent = CenterCurrentINPUT;  // Array for current centers

            KmeansTriangleInequality.Ncent_Global = Ncent_GlobalINPUT;
            KmeansTriangleInequality.MaxNcent_Global = MaxNcent_GlobalINPUT;
            KmeansTriangleInequality.Ncent_Global_Parallel = Ncent_Global_ParallelINPUT;
            KmeansTriangleInequality.ParameterVectorDimension = ParameterVectorDimensionINPUT;

            KmeansTriangleInequality.SetParallelCenterDecomposition();

            //  Some Sizes
            double CenterSize = KmeansTriangleInequality.ParameterVectorDimension * KmeansTriangleInequality.MaxNcent_Global * 8;
            double CenterSizeNode = CenterSize * 3 * DAVectorUtility.MPIperNodeCount;
            DAVectorUtility.SALSAPrint(0, " One Center Array " + CenterSize.ToString("F0") + " Per node " + CenterSizeNode.ToString("F0") + " " + CenterSizeNode.ToString("E4"));
            double CenterCenterSize = KmeansTriangleInequality.MaxNcent_Global * MaxCentersperCenter * 12;
            double CenterCenterSizeNode = CenterCenterSize * DAVectorUtility.MPIperNodeCount;
            DAVectorUtility.SALSAPrint(0, " Center Center Array " + CenterCenterSize.ToString("F0") + " Per node " + CenterCenterSizeNode.ToString("F0") + " " + CenterCenterSizeNode.ToString("E4"));
            double PointSize = DAVectorUtility.PointCount_Largest * KmeansTriangleInequality.ParameterVectorDimension * 8;
            double PointSizeNode = PointSize * DAVectorUtility.MPIperNodeCount;
            DAVectorUtility.SALSAPrint(0, " Point Vector Array " + PointSize.ToString("F0") + " Per node " + PointSizeNode.ToString("F0") + " " + PointSizeNode.ToString("E4"));
            double PointLBSize = DAVectorUtility.PointCount_Largest * KmeansTriangleInequality.MaxClusterLBsperPoint * 12;
            if (KmeansTriangleInequality.OldCenterOption >= 0)
                PointLBSize *= 3;
            else
                PointLBSize *= 2;
            double PointLBSizeNode = PointLBSize * DAVectorUtility.MPIperNodeCount;
            DAVectorUtility.SALSAPrint(0, " Point LB Array " + PointLBSize.ToString("F0") + " Per node " + PointLBSizeNode.ToString("F0") + " " + PointLBSizeNode.ToString("E4"));
            double FindClusterCenterSize = CenterSize * DAVectorUtility.ThreadCount;
            double FindClusterCenterSizeNode = FindClusterCenterSize * DAVectorUtility.MPIperNodeCount;
            DAVectorUtility.SALSAPrint(0, " Find Cluster Center Size " + FindClusterCenterSize.ToString("F0") + " Per node " + FindClusterCenterSize.ToString("F0") + " " + FindClusterCenterSize.ToString("E4"));
            double TotalSize = CenterSizeNode + CenterCenterSizeNode + PointSizeNode + PointLBSizeNode + FindClusterCenterSizeNode;

            DAVectorUtility.SALSAPrint(0, " Total per node " + TotalSize.ToString("F0") + " " + TotalSize.ToString("E4"));

            //  Center Related Data
            KmeansTriangleInequality.CenterOld = new double[KmeansTriangleInequality.MaxNcent_Global][];
            KmeansTriangleInequality.CenterLast = new double[KmeansTriangleInequality.MaxNcent_Global][];
            KmeansTriangleInequality.CurrentInterCenter = new CenterData[KmeansTriangleInequality.MaxNcent_Global];
            KmeansTriangleInequality.DistributedCenterDifference = new double[KmeansTriangleInequality.CenterCount_Process][];
            for (int LocalCenterIndex = 0; LocalCenterIndex < KmeansTriangleInequality.CenterCount_Process; LocalCenterIndex++)
                KmeansTriangleInequality.DistributedCenterDifference[LocalCenterIndex] = new double[KmeansTriangleInequality.MaxNcent_Global];

            if (KmeansTriangleInequality.UseParallelismoverCenters)
            {   // Centers Parallel over Threads NOT nodes
                Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
                {
                    int NumberCentersperThread = KmeansTriangleInequality.LocalParallel_CentersperThread[ThreadIndex];
                    int BeginCenter = KmeansTriangleInequality.LocalParallel_StartCenterperThread[ThreadIndex];
                    for (int CenterIndex = BeginCenter; CenterIndex < NumberCentersperThread + BeginCenter; CenterIndex++)
                    {
                        KmeansTriangleInequality.CenterOld[CenterIndex] = new double[KmeansTriangleInequality.ParameterVectorDimension];
                        KmeansTriangleInequality.CenterLast[CenterIndex] = new double[KmeansTriangleInequality.ParameterVectorDimension];
                        KmeansTriangleInequality.CurrentInterCenter[CenterIndex] = new CenterData();
                    }
                }); // End Sum over Threads
            }
            else
            {   // Centers Sequential
                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.MaxNcent_Global; CenterIndex++)
                {
                    KmeansTriangleInequality.CenterOld[CenterIndex] = new double[KmeansTriangleInequality.ParameterVectorDimension];
                    KmeansTriangleInequality.CenterLast[CenterIndex] = new double[KmeansTriangleInequality.ParameterVectorDimension];
                    KmeansTriangleInequality.CurrentInterCenter[CenterIndex] = new CenterData();
                }
            }

            // Point Related Data
            KmeansTriangleInequality.NearestCentertoPoint = new int[DAVectorUtility.PointCount_Process];
            KmeansTriangleInequality.Distance_NearestCentertoPoint = new double[DAVectorUtility.PointCount_Process];
            if (KmeansTriangleInequality.OldCenterOption >= 0)
                KmeansTriangleInequality.LB_Old_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];
            KmeansTriangleInequality.LB_Buffer1_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];
            KmeansTriangleInequality.LB_Buffer2_alpha_ClusterPointer = new PointData[DAVectorUtility.PointCount_Process];

            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha] = new PointData(KmeansTriangleInequality.MaxClusterLBsperPoint);
                    KmeansTriangleInequality.LB_Buffer1_alpha_ClusterPointer[alpha] = new PointData(KmeansTriangleInequality.MaxClusterLBsperPoint);
                    KmeansTriangleInequality.LB_Buffer2_alpha_ClusterPointer[alpha] = new PointData(KmeansTriangleInequality.MaxClusterLBsperPoint);   
                }
            }); // End Loop over Threads for points

            KmeansTriangleInequality.LB_Current_alpha_ClusterPointer = KmeansTriangleInequality.LB_Buffer1_alpha_ClusterPointer;
            KmeansTriangleInequality.LB_Last_alpha_ClusterPointer = null;

        }   // End Initialize()

        public static void SetParallelCenterDecomposition()
        {
            KmeansTriangleInequality.UseParallelismoverCenters = false;
            if (KmeansTriangleInequality.Ncent_Global <= KmeansTriangleInequality.Ncent_Global_Parallel)
            {
                KmeansTriangleInequality.CenterStart_Process = 0;
                KmeansTriangleInequality.CenterCount_Process = KmeansTriangleInequality.MaxNcent_Global;
                KmeansTriangleInequality.CenterCount_Largest = KmeansTriangleInequality.MaxNcent_Global;
                return;
            }
            KmeansTriangleInequality.UseParallelismoverCenters = true;
            DAVectorUtility.SALSAPrint(0, "Center Parallelism used as well as Point Parallelism in Triangle Inequality");

            //	First divide centers among processes
            Range[] processRanges = RangePartitioner.Partition(KmeansTriangleInequality.MaxNcent_Global, DAVectorUtility.MPI_Size);
            Range processRange = processRanges[DAVectorUtility.MPI_Rank];  // The answer for this process

            KmeansTriangleInequality.CenterStart_Process = processRange.StartIndex;
            KmeansTriangleInequality.CenterCount_Process = processRange.Length;

            KmeansTriangleInequality.CenterCount_Largest = int.MinValue;
            foreach (Range r in processRanges)
                KmeansTriangleInequality.CenterCount_Largest = Math.Max(r.Length, KmeansTriangleInequality.CenterCount_Largest);


            //	Now divide centers among threads for this process
            Range[] Full_ThreadRanges = RangePartitioner.Partition(processRange, DAVectorUtility.ThreadCount);
            KmeansTriangleInequality.FullParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
            KmeansTriangleInequality.FullParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

            for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
            {
                KmeansTriangleInequality.FullParallel_CentersperThread[j] = Full_ThreadRanges[j].Length;
                KmeansTriangleInequality.FullParallel_StartCenterperThread[j] = Full_ThreadRanges[j].StartIndex;
            }

            //  Process only Center Parallelism
            Range[] Local_ThreadRanges = RangePartitioner.Partition(KmeansTriangleInequality.MaxNcent_Global, DAVectorUtility.ThreadCount);
            KmeansTriangleInequality.LocalParallel_CentersperThread = new int[DAVectorUtility.ThreadCount];
            KmeansTriangleInequality.LocalParallel_StartCenterperThread = new int[DAVectorUtility.ThreadCount];

            for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
            {
                KmeansTriangleInequality.LocalParallel_CentersperThread[j] = Local_ThreadRanges[j].Length;
                KmeansTriangleInequality.LocalParallel_StartCenterperThread[j] = Local_ThreadRanges[j].StartIndex;
            }

        }   // End SetParallelCenterDecomposition()

        public static void NextIteration()
        {   //  Move to next iteration

            //  Center Related
            if (!KmeansTriangleInequality.CenterCurrentAlreadyset)
            {
                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    KmeansTriangleInequality.GetIndividualCenterPosition(CenterIndex, KmeansTriangleInequality.CenterCurrent[CenterIndex]);   // Set Pointer to Current Centers set outside
            }

            if (KmeansTriangleInequality.UseTriangleInequality != 0)
            {
                KmeansTriangleInequality.LB_Last_alpha_ClusterPointer = KmeansTriangleInequality.LB_Current_alpha_ClusterPointer;
                if (KmeansTriangleInequality.LB_Last_alpha_ClusterPointer.Equals(KmeansTriangleInequality.LB_Buffer1_alpha_ClusterPointer))
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer = KmeansTriangleInequality.LB_Buffer2_alpha_ClusterPointer;
                else
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer = KmeansTriangleInequality.LB_Buffer1_alpha_ClusterPointer;

            }

            //  Set up Diagnostics
            KmeansTriangleInequality.NumberDiagnosticSums_Center = 9; 
            KmeansTriangleInequality.FindDiagnosticSums_Center = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.NumberDiagnosticSums_Center);
            KmeansTriangleInequality.NumberDiagnosticSums_Points = 27; 
            KmeansTriangleInequality.FindDiagnosticSums_Points = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.NumberDiagnosticSums_Points);

            if (KmeansTriangleInequality.UseTriangleInequality != 0)
            {
                KmeansTriangleInequality.RecalculateOldCenters = false;
                if( (KmeansTriangleInequality.OldCenterOption > 2) && (KmeansTriangleInequality.CalculationIterations > 5) && (KmeansTriangleInequality.CalculationIterations%KmeansTriangleInequality.OldCenterOption == 0) )
                    KmeansTriangleInequality.RecalculateOldCenters = true;

                // Analyze over centers
                DAVectorUtility.StartSubTimer(2);
                KmeansTriangleInequality.CenterFacingAnalysis();
                DAVectorUtility.StopSubTimer(2);
                // Analyze over points
                DAVectorUtility.StartSubTimer(3);
                if (KmeansTriangleInequality.StartingOut)
                    KmeansTriangleInequality.InitialPointAnalysis();
                else
                    KmeansTriangleInequality.PointDependentAnalysis();
                DAVectorUtility.StopSubTimer(3);
            }
            else
            {
                DAVectorUtility.StartSubTimer(0);
                KmeansTriangleInequality.PureKmeans();
                DAVectorUtility.StopSubTimer(0);
            }

            //  Process Diagnostics
            if (!KmeansTriangleInequality.UseParallelismoverCenters)
                DAVectorUtility.MPI_Size = 1;
            KmeansTriangleInequality.FindDiagnosticSums_Center.sumoverthreadsandmpi();
            if (!KmeansTriangleInequality.UseParallelismoverCenters)
                DAVectorUtility.MPI_Size = KmeansTriangleInequality.SaveMPI_Size;
            KmeansTriangleInequality.FindDiagnosticSums_Points.sumoverthreadsandmpi();

            KmeansTriangleInequality.NumberFullDistancesCalculatedCC += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[0];
            KmeansTriangleInequality.CalculationCenterIterations += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[1];
            KmeansTriangleInequality.CountToobigAShift += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[2];
            KmeansTriangleInequality.AccumulateRadius += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[3];
            KmeansTriangleInequality.AccumulateDmax += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[4];
            KmeansTriangleInequality.CountCCEntries += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[5];
            KmeansTriangleInequality.AccumulateCCvalues += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[6];
            KmeansTriangleInequality.AccumulateLastCenterchanges += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[7];
            KmeansTriangleInequality.AccumulateOldCenterchanges += KmeansTriangleInequality.FindDiagnosticSums_Center.TotalVectorSum[8];

            KmeansTriangleInequality.NumberFullDistancesCalculatedPC += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[0];   // Total
            KmeansTriangleInequality.CalculationPointIterations += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[1];   // Accumulate Product of Points and Iterations
            KmeansTriangleInequality.AccumulateOLD_LB += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[2];   // OLD Lower Bound used successfully
            KmeansTriangleInequality.AccumulateCURRENT_LB += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[3];   // CURRENT Lower Bound used successfully
            KmeansTriangleInequality.CountLBEntries += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[4];  //  Number in Lower Bound Array
            KmeansTriangleInequality.AccumulateBestScore += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[5];  // Best Score
            KmeansTriangleInequality.AccumulateMAX_LB += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[6];  // Maximum Lower Bound
            KmeansTriangleInequality.AccumulateBreakReason7 += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[7];  // Break in LB while as all centers done
            KmeansTriangleInequality.AccumulateBreakReason8 += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[8];  // Break in LB while as Center Scan Dmax
            KmeansTriangleInequality.AccumulateBreakReason9 += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[9];  // Break in LB while as Center Scan middle
            KmeansTriangleInequality.AccumulateBreakReason10 += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[10];  // Break in LB while Center Scan start
            KmeansTriangleInequality.NumberFullDistancesCalculatedPC_Old += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[11];   // Total PC calculation in Old Lower Bound Update
            KmeansTriangleInequality.NumberFullDistancesCalculatedPC_Last += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[12];   // Total PC calculation in Last Lower Bound Update
            KmeansTriangleInequality.AccumulateOLD_LBFail += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[13];   // OLD Lower Bound used UNsuccessfully
            KmeansTriangleInequality.AccumulateCURRENT_LBFail += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[14];   // CURRENT Lower Bound used UNsuccessfully
            KmeansTriangleInequality.NumberFullDistancesCalculatedLoop += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[15];   // Number of Full Distances calculated in main lowwer bound loop
            KmeansTriangleInequality.NumberCentersinBestSet += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[16];  // Number of centers examined in best loop Method = -2
            KmeansTriangleInequality.NumberCentersinSimpleSet += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[17];    // Number of centers examined in Method = -1 loop
            KmeansTriangleInequality.NumberCentersinCenterbasedset += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[18];   // Number of centers examined in method >=0 loop
            KmeansTriangleInequality.AverageCentersLookedAt += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[19];   // Average number of centers NOT looked at
            KmeansTriangleInequality.NumberCentersRemovedBackwards += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[20]; // Centers Removed in Backwards test
            KmeansTriangleInequality.NumberBackwardTests += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[21]; // Number of backward tests
            KmeansTriangleInequality.NumberForwardTests += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[22]; // Number of forward tests
            KmeansTriangleInequality.PossibleNumberFullDistancesCalculated += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[23];   // Possible Number of Full Distances calculated in main lowwer bound loop
            KmeansTriangleInequality.AccumulateCURRENT_LBLimit += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[24];    // Center veto used Limit on Rest
            KmeansTriangleInequality.AccumulateBreakLimit += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[25];    // Break Limit Value
            KmeansTriangleInequality.AccumulateOLD_LBLimit += KmeansTriangleInequality.FindDiagnosticSums_Points.TotalVectorSum[26];    // Center veto used Limit on Rest

            //  End Iteration
            DAVectorUtility.StartSubTimer(1);
            KmeansTriangleInequality.EndIteration();
            DAVectorUtility.StopSubTimer(1);

        }   // End NextIteration()

        public static void EndIteration()
        {
            if (KmeansTriangleInequality.UseParallelismoverCenters)
            {   // Centers Parallel over Threads NOT nodes
                Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
                {
                    int NumberCentersperThread = KmeansTriangleInequality.LocalParallel_CentersperThread[ThreadIndex];
                    int BeginCenter = KmeansTriangleInequality.LocalParallel_StartCenterperThread[ThreadIndex];
                    for (int CenterIndex = BeginCenter; CenterIndex < NumberCentersperThread + BeginCenter; CenterIndex++)
                    {
                        for (int VectorIndex = 0; VectorIndex < KmeansTriangleInequality.ParameterVectorDimension; VectorIndex++)
                            KmeansTriangleInequality.CenterLast[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                    }
                }); // End Sum over Threads
            }
            else
            {   // Centers Sequential
                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                {
                    for (int VectorIndex = 0; VectorIndex < KmeansTriangleInequality.ParameterVectorDimension; VectorIndex++)
                        KmeansTriangleInequality.CenterLast[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                }
            }
            KmeansTriangleInequality.StartingOut = false;
            //  Iteration Count
            KmeansTriangleInequality.CalculationIterations++;

            // Reset Cluster Centers
            DAVectorUtility.StartSubTimer(4);
            KmeansTriangleInequality.FindClusterCenters(false, KmeansTriangleInequality.NearestCentertoPoint, KmeansTriangleInequality.Distance_NearestCentertoPoint, KmeansTriangleInequality.CenterLast);
            DAVectorUtility.StopSubTimer(4);

        }   // End EndIteration()

        public static void CenterFacingAnalysis()
        {   // This does NOT get factor of 2 speed up from symmetry of distance computation

            //  First set O(Center) changes using thread not distributed parallelism
            // Distributed version could be used
            if (KmeansTriangleInequality.UseParallelismoverCenters)
            {   // Centers Parallel over Threads NOT nodes

                DAVectorUtility.StartSubTimer(5);
                double[] AccumulateCC = new double[DAVectorUtility.ThreadCount];
                double[] AccumulateTooBig = new double[DAVectorUtility.ThreadCount];
                double[] AccumulateRadius = new double[DAVectorUtility.ThreadCount];
                for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                {
                    AccumulateCC[ThreadIndex] = 0.0;
                    AccumulateTooBig[ThreadIndex] = 0.0;
                    AccumulateRadius[ThreadIndex] = 0.0;
                }

                Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
                {
                    int NumberCentersperThread = KmeansTriangleInequality.LocalParallel_CentersperThread[ThreadIndex];
                    int BeginCenter = KmeansTriangleInequality.LocalParallel_StartCenterperThread[ThreadIndex];
                    for (int CenterIndex = BeginCenter; CenterIndex < NumberCentersperThread + BeginCenter; CenterIndex++)
                    {
                        KmeansTriangleInequality.CFA_InitialCenterCalculation(CenterIndex, ref AccumulateCC[ThreadIndex], ref AccumulateTooBig[ThreadIndex],
                            ref AccumulateRadius[ThreadIndex]);
                    }
                }); // End Sum over Threads

                for (int ThreadIndex = 0; ThreadIndex < DAVectorUtility.ThreadCount; ThreadIndex++)
                {
                    KmeansTriangleInequality.NumberFullDistancesCalculatedCC += AccumulateCC[ThreadIndex];
                    KmeansTriangleInequality.CountToobigAShift += AccumulateTooBig[ThreadIndex];
                    KmeansTriangleInequality.AccumulateRadius += AccumulateRadius[ThreadIndex];
                }
                DAVectorUtility.StopSubTimer(5);

            }   // End Center Parallel over threads

            else
            {   // Centers Sequential
                DAVectorUtility.StartSubTimer(5);
                for (int CenterIndex = 0; CenterIndex < Ncent_Global; CenterIndex++)
                    KmeansTriangleInequality.CFA_InitialCenterCalculation(CenterIndex, ref KmeansTriangleInequality.NumberFullDistancesCalculatedCC, ref KmeansTriangleInequality.CountToobigAShift,
                        ref KmeansTriangleInequality.AccumulateRadius);
                DAVectorUtility.StopSubTimer(5);

            }   // End Sequential center analysis doing things O(Number of Centers)

            //  Reconcile TriangleInequality.CountToobigAShift and TriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus if MPI Parallel
            if (DAVectorUtility.MPI_Size > 1)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<double>(ref KmeansTriangleInequality.CountToobigAShift, 0);
                int[] TransmitRecalculateStatus = new int[KmeansTriangleInequality.Ncent_Global];
                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    TransmitRecalculateStatus[CenterIndex] = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus;
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref TransmitRecalculateStatus, 0);
                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus = TransmitRecalculateStatus[CenterIndex];
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            }

            //  Now do full analysis with O(Number of Centers ^2)
            DAVectorUtility.StartSubTimer(6);
            int NumCentersinOrderedList = Math.Min(KmeansTriangleInequality.Ncent_Global - 1, KmeansTriangleInequality.MaxCentersperCenter);
            if (KmeansTriangleInequality.UseParallelismoverCenters)
            {   // Centers Parallel over Threads and processes

                //  Broadcast Data
                int IntSize = 3 + (CenterData.NumberPackingInts + NumCentersinOrderedList) * KmeansTriangleInequality.CenterCount_Process;
                int DoubleSize = (CenterData.NumberPackingDoubles + NumCentersinOrderedList) * KmeansTriangleInequality.CenterCount_Process;
                int IntSize1 = 3 + (CenterData.NumberPackingInts + NumCentersinOrderedList) * KmeansTriangleInequality.CenterCount_Largest;
                int DoubleSize1 = (CenterData.NumberPackingDoubles + NumCentersinOrderedList) * KmeansTriangleInequality.CenterCount_Largest;
                int[] IntSendTransport = new int[IntSize1];
                double[] DoubleSendTransport = new double[DoubleSize1];

                //  Parallel Center Processing
                Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
                {
                    int NumberCentersToProcess = KmeansTriangleInequality.FullParallel_CentersperThread[ThreadIndex];
                    int BeginCenterThread = KmeansTriangleInequality.FullParallel_StartCenterperThread[ThreadIndex];
                    for (int CenterIndex = BeginCenterThread; CenterIndex < NumberCentersToProcess + BeginCenterThread; CenterIndex++)
                    {
                        KmeansTriangleInequality.CFA_CenterCenterCalculation(ThreadIndex, KmeansTriangleInequality.CenterStart_Process, CenterIndex, NumCentersinOrderedList);
                        int intposition = 3 + (CenterData.NumberPackingInts + NumCentersinOrderedList) * (CenterIndex - KmeansTriangleInequality.CenterStart_Process);
                        int doubleposition = (CenterData.NumberPackingDoubles + NumCentersinOrderedList) * (CenterIndex - KmeansTriangleInequality.CenterStart_Process);
                        KmeansTriangleInequality.CurrentInterCenter[CenterIndex].serialize(IntSendTransport, DoubleSendTransport, intposition, doubleposition);
                    }

                }); // End Loop over Threads
                IntSendTransport[0] = KmeansTriangleInequality.CenterStart_Process;
                IntSendTransport[1] = KmeansTriangleInequality.CenterCount_Process;
                IntSendTransport[2] = DAVectorUtility.MPI_Rank;


                int[] IntReceiveTransport = new int[IntSize1];
                double[] DoubleReceiveTransport = new double[DoubleSize1];
                
                //  Do a broadcast from each process
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                for (int mpiloop = 0; mpiloop < DAVectorUtility.MPI_Size; mpiloop++)
                {
                    int[] IntSend = IntSendTransport;
                    if (mpiloop != DAVectorUtility.MPI_Rank)
                        IntSend = IntReceiveTransport;
                    //DAVectorUtility.SALSAFullPrint(0, mpiloop + " " + IntSend.Length + " " + IntSize1 + " " + DoubleSize1 + " " + IntSize + " " + DoubleSize + " " + KmeansTriangleInequality.CenterCount_Process + " " + KmeansTriangleInequality.CenterCount_Largest);

                    DAVectorUtility.MPI_communicator.Broadcast<int>(ref IntSend, (int) mpiloop);
                    if(IntSend[2] != mpiloop)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Inconsistent Broadcast Data from rank " + IntSend[2].ToString() + " " + mpiloop.ToString());
                        throw (e);
                    }

                    double[] DoubleSend = DoubleSendTransport;
                    if (mpiloop != DAVectorUtility.MPI_Rank)
                        DoubleSend = DoubleReceiveTransport;
                    DAVectorUtility.MPI_communicator.Broadcast<double>(ref DoubleSend, (int) mpiloop);

                    if (mpiloop == DAVectorUtility.MPI_Rank)
                        continue;

                    int BeginCenter = IntReceiveTransport[0];
                    int NumCentersTransported = IntReceiveTransport[1];
                    for (int CenterIndex = BeginCenter; CenterIndex < NumCentersTransported + BeginCenter; CenterIndex++)
                    {
                        int intposition = 3 + (CenterData.NumberPackingInts + NumCentersinOrderedList) * (CenterIndex - BeginCenter);
                        int doubleposition = (CenterData.NumberPackingDoubles + NumCentersinOrderedList) * (CenterIndex - BeginCenter);
                        KmeansTriangleInequality.CurrentInterCenter[CenterIndex].deserialize(IntReceiveTransport, DoubleReceiveTransport, intposition, doubleposition);
                    }
                }
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            }   // End Parallel Centers analysis of O(Number of Centers ^2) calculation

            else
            {   // Centers Sequential

                for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    KmeansTriangleInequality.CFA_CenterCenterCalculation(0, 0, CenterIndex, NumCentersinOrderedList);

            }   // End Sequential Centers analysis of O(Number of Centers ^2) calculation
            DAVectorUtility.StopSubTimer(6);

            //  Loop over all centers addressing TriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus == 1
            if (KmeansTriangleInequality.StartingOut)
                return;
            DAVectorUtility.StartSubTimer(7);
            for (int CenterIndex = 0; CenterIndex < Ncent_Global; CenterIndex++)
            {
                if (KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus == 0)
                    continue;
                KmeansTriangleInequality.CurrentInterCenter[CenterIndex].LastCenterChange = 0.0;
                KmeansTriangleInequality.CurrentInterCenter[CenterIndex].OldCenterChange = 0.0;
                for (int VectorIndex = 0; VectorIndex < KmeansTriangleInequality.ParameterVectorDimension; VectorIndex++)
                {
                    KmeansTriangleInequality.CenterOld[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                    KmeansTriangleInequality.CenterLast[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                }
            }
            DAVectorUtility.StopSubTimer(7);

        }   // End CenterFacingAnalysis()

        //  Initial O(Number of Centers) Computation
        public static void CFA_InitialCenterCalculation(int CenterIndex, ref double ThreadCCAccumulate, ref double ThreadTooBigAccumulate, ref double ThreadRadiusAccumulate)
        {
            double CurrentLastDifference = 0.0;
            double CurrentOldDifference = 0.0;
            double Radius = 0.0;
            bool TooBigaShift = false;

            if (KmeansTriangleInequality.StartingOut)
            {   // Set Old Centers
                for (int VectorIndex = 0; VectorIndex < KmeansTriangleInequality.ParameterVectorDimension; VectorIndex++)
                {
                    KmeansTriangleInequality.CenterOld[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                    KmeansTriangleInequality.CenterLast[CenterIndex][VectorIndex] = KmeansTriangleInequality.CenterCurrent[CenterIndex][VectorIndex];
                }
            }

            else
            {
                Radius = KmeansTriangleInequality.GetRadiusperCenter(CenterIndex);
                CurrentLastDifference = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.CenterLast[CenterIndex]);
                CurrentOldDifference = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.CenterOld[CenterIndex]);
                ThreadCCAccumulate += 2.0;

                if ((CurrentLastDifference > KmeansTriangleInequality.TriangleInequality_Delta1_current * Radius) ||
                    (CurrentOldDifference > KmeansTriangleInequality.TriangleInequality_Delta1_old * Radius))
                {
                    TooBigaShift = true;
                    ThreadTooBigAccumulate += 1.0;
                }
            }

            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus = 0;
            if (TooBigaShift)
                KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus = 1;
            if (KmeansTriangleInequality.RecalculateOldCenters )
               KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus = 1;

            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].ClusterRadius = Radius;
            ThreadRadiusAccumulate += Radius;
            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].LastCenterChange = CurrentLastDifference;
            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].OldCenterChange = CurrentOldDifference;

        }   // End CFA_IndividualCenterCalculation

        public static void CFA_CenterCenterCalculation(int ThreadIndex, int CenterStart, int CenterIndex, int NumCentersinOrderedList)
        {   // Do calculation for one center in O(Number of Centers ^2)  part of CenterFacingAnalysis()

            KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, KmeansTriangleInequality.CurrentInterCenter[CenterIndex].LastCenterChange, 7);
            KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 1);
            bool MainChanged = false;
            if (KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0)
                MainChanged = true;
            double CurrentOldDifference = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].OldCenterChange;
            KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, CurrentOldDifference, 8);

            //  Set up ordered matrix of inter-center distances
            int currentcut = -1;
            double Minrejected = -1.0;
            int ActualNumber = 0;
            double[] SmallValues = new double[NumCentersinOrderedList];
            int[] SmallIndices = new int[NumCentersinOrderedList];

            //  Prepare Dynamic Looping
            bool[] LookedAt = new bool[KmeansTriangleInequality.Ncent_Global];
            for (int CenterIndex1 = 0; CenterIndex1 < KmeansTriangleInequality.Ncent_Global; CenterIndex1++)
                LookedAt[CenterIndex1] = false;
            int NumCentersToGo = KmeansTriangleInequality.Ncent_Global;  // Number of Centers still to look at -- end when <= 0

            int CurrentCandidateCenter = -1;    // Current Center being looked at in Method = -1 approach
            int Method = -2;    // Method = -2 Start with Previous list, Method = -1 scan List in order of centers
            int Methodminus2Size = 0;
            if (KmeansTriangleInequality.StartingOut)
                Method = -1;
            else
                Methodminus2Size = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].NearbyClusters;
            int MethodPosition = -1; // Position in Center List of Center Method = -2
            int CenterIndexPrime = -1;
            while(NumCentersToGo > 0)
            {
                if (Method == -1)
                {
                    ++CurrentCandidateCenter;
                    CenterIndexPrime = CurrentCandidateCenter;
                }
                else // Method = -2
                {
                    ++MethodPosition;
                    if (MethodPosition >= Methodminus2Size)
                    {
                        Method = -1;
                        continue;
                    }
                    CenterIndexPrime = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].CenterDifferenceIndices[MethodPosition];
                }

                if (LookedAt[CenterIndexPrime])
                    continue;
                LookedAt[CenterIndexPrime] = true;
                NumCentersToGo--;

                if (CenterIndexPrime == CenterIndex)
                {
                    KmeansTriangleInequality.DistributedCenterDifference[CenterIndex - CenterStart][CenterIndexPrime] = 0.0;
                    continue;
                }
                bool PrimeChanged = KmeansTriangleInequality.CurrentInterCenter[CenterIndexPrime].RecalculateStatus > 0;
                if ((!KmeansTriangleInequality.StartingOut) && (!MainChanged) && (!PrimeChanged) && (currentcut >= 0) && (ActualNumber == NumCentersinOrderedList))
                {
                    double AtLeast = KmeansTriangleInequality.DistributedCenterDifference[CenterIndex - CenterStart][CenterIndexPrime]  - CurrentOldDifference
                        - KmeansTriangleInequality.CurrentInterCenter[CenterIndexPrime].OldCenterChange;
                    if (AtLeast >= SmallValues[currentcut])
                    {
                        if (Minrejected < 0.0)
                            Minrejected = AtLeast;
                        else
                            Minrejected = Math.Min(Minrejected, AtLeast);
                        continue;
                    }
                }
                double Centerdistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.CenterCurrent[CenterIndexPrime]);
                KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 0);
                KmeansTriangleInequality.FindMinimumSetwithRemainder(Centerdistce, CenterIndexPrime, ref currentcut, SmallValues, SmallIndices, NumCentersinOrderedList, ref ActualNumber, ref Minrejected);

                if (KmeansTriangleInequality.StartingOut || (MainChanged && PrimeChanged))
                    KmeansTriangleInequality.DistributedCenterDifference[CenterIndex - CenterStart][CenterIndexPrime] = Centerdistce;
                else
                {
                    if ((!MainChanged) && PrimeChanged)
                    {
                        KmeansTriangleInequality.DistributedCenterDifference[CenterIndex - CenterStart][CenterIndexPrime] =
                            DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterOld[CenterIndex], KmeansTriangleInequality.CenterCurrent[CenterIndexPrime]);
                        KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 0);
                    }
                    if ( MainChanged && (!PrimeChanged) )
                    {
                        KmeansTriangleInequality.DistributedCenterDifference[CenterIndex - CenterStart][CenterIndexPrime] =
                            DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.CenterOld[CenterIndexPrime]);
                        KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, 1.0, 0);
                    }
                }

            }   // End While loop over centers
            Array.Sort(SmallValues, SmallIndices, 0, ActualNumber);
               
            //  Set Data Structure
            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].NearbyClusters = ActualNumber;
            KmeansTriangleInequality.CurrentInterCenter[CenterIndex].Dmax = Minrejected;
            KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, Minrejected, 4);
            KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, (double) ActualNumber, 5);

            for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
            {
                KmeansTriangleInequality.CurrentInterCenter[CenterIndex].CenterDifference[loopindex] = SmallValues[loopindex];
                KmeansTriangleInequality.FindDiagnosticSums_Center.addapoint(ThreadIndex, SmallValues[loopindex], 6);
                KmeansTriangleInequality.CurrentInterCenter[CenterIndex].CenterDifferenceIndices[loopindex] = SmallIndices[loopindex];
            }

        }   // End CFA_CenterCenterCalculation

        public static void InitialPointAnalysis()
        {
            if (KmeansTriangleInequality.OldCenterOption > 0)
            {
                InitialPointAnalysis_FullOld();
                return;
            }
            
            //  Set Associated Cluster and initial lower bounds before any previous point data available
            //  It loops through centers from beginning to end
            //  When a new best center is found, it changes to loop over associated centers
            //  If best center loop exhausted it reverts to ploughing through center list

            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 1); 
                    
                    int NumCentersToGo = KmeansTriangleInequality.Ncent_Global;  // Number of Centers still to look at -- end when <= 0
                    bool[] LookedAt = new bool[KmeansTriangleInequality.Ncent_Global];
                    for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                        LookedAt[CenterIndex] = false;

                    int CurrentCandidateCenter = -1;    // Current Center being looked at in Method = -1 (scan centers in numerical order) approach

                    int Method = -1;    // Method = -1 scan Center List; Method > 0 Scan centers in list of center Method
                    int MethodPosition = -1; // Position in Center List of Center Method

                    int BestCenter = -2;    // Index of nearest center
                    double BestScore = -1.0;    // Distance of nearest center
                    int ExaminedCenter = -1;    // in while loop, number of center being examined
                    int BreakReason = 7;    // Reason for Break -- 7 is hit continue

                    int Num_LBValuesforPoint = KmeansTriangleInequality.MaxClusterLBsperPoint;    // Max Size of ordered list preserved in LB_Current_alpha_ClusterPointer
                    int ActualNumber = 0;   // For creating ordered list preserved in LB_Current_alpha_ClusterPointer
                    int currentcut = -1;    // For creating ordered list preserved in LB_Current_alpha_ClusterPointer
                    double Minrejected = -1.0;  // For creating ordered list preserved in LB_Current_alpha_ClusterPointer

                    int ActuallyUse = KmeansTriangleInequality.Ncent_Global;
                    double[] SmallValues = new double[ActuallyUse];    // For creating ordered list preserved in LB_Current_alpha_ClusterPointer
                    int[] SmallIndices = new int[ActuallyUse]; // For creating ordered list preserved in LB_Current_alpha_ClusterPointer

                    Double BreakLimit = -1.0;    // Lower Bound on Stopping

                    //  Loop over candidates
                    while (NumCentersToGo > 0)
                    {
                        if (Method == -1)
                        {
                            CurrentCandidateCenter++;
                            if (LookedAt[CurrentCandidateCenter])
                                continue;
                            ExaminedCenter = CurrentCandidateCenter;
                        }
                        else
                        {   // Examining centers near current best center
                            ++MethodPosition;
                            if (MethodPosition >= KmeansTriangleInequality.CurrentInterCenter[Method].NearbyClusters)
                            {
                                if (KmeansTriangleInequality.CurrentInterCenter[Method].Dmax > 2.0 * BestScore)
                                {
                                    BreakReason = 8;
                                    BreakLimit = KmeansTriangleInequality.CurrentInterCenter[Method].Dmax - BestScore;
                                    break;
                                }
                                Method = -1;
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 17);
                                continue;
                            }
                            if (KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] >= 2.0 * BestScore)
                            {
                                BreakReason = 9;
                                BreakLimit = KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] - BestScore;
                                break;
                            }
                            ExaminedCenter = KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifferenceIndices[MethodPosition];
                            if (LookedAt[ExaminedCenter])
                                continue;
                            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 18);
                        }

                        //  Need to look at explicit distance
                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterCurrent[ExaminedCenter], KmeansTriangleInequality.PointPosition[alpha]);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 15);

                        KmeansTriangleInequality.FindMinimumSetwithRemainder(PointCenterDistce, ExaminedCenter, ref currentcut, SmallValues, SmallIndices, ActuallyUse, 
                            ref ActualNumber, ref Minrejected);
                        if( (BestCenter >= 0) && (PointCenterDistce > BestScore) )
                        {
                            --NumCentersToGo;
                            LookedAt[ExaminedCenter] = true;
                            continue;
                        }

                        //  New best Center discovered
                        BestScore = PointCenterDistce;
                        BestCenter = ExaminedCenter;
                        if (KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] >= 2.0 * PointCenterDistce)
                        {
                            BreakReason = 10;
                            BreakLimit = KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] - BestScore;
                            break;
                        }
                        Method = ExaminedCenter;
                        MethodPosition = -1;
                        --NumCentersToGo;
                        LookedAt[ExaminedCenter] = true;
                        continue;

                    }   // End while over centers
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, BreakReason);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double)NumCentersToGo, 19);

                    BreakLimit = Math.Max(BestScore, BreakLimit);
                    Array.Sort(SmallValues, SmallIndices, 0, ActualNumber);
                    if (Num_LBValuesforPoint < ActualNumber)
                    {
                        Minrejected = SmallValues[Num_LBValuesforPoint];
                        if (BreakLimit < 0.0)
                            BreakLimit = Minrejected;
                        else
                            BreakLimit = Math.Min(BreakLimit, Minrejected);
                        ActualNumber = Num_LBValuesforPoint;
                    }
                    

                    if (BestCenter != SmallIndices[0])
                    {
                        if (Math.Abs(BestScore - SmallValues[0]) < 0.001 * BestScore)
                        {
                            BestCenter = SmallIndices[0];
                            BestScore = SmallValues[0];
                        }
                        else
                        {
                            Exception e = DAVectorUtility.SALSAError("Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " in Process " + DAVectorUtility.MPI_Rank.ToString()
                                + " Inconsistent Best Centers 1 " + SmallIndices[0].ToString() + " " + BestCenter.ToString()
                                + " " + SmallValues[0].ToString("F4") + " " + BestScore.ToString("F4"));
                            throw (e);
                        }
                    }
                    if (Math.Abs(BestScore - SmallValues[0]) > 0.001 * BestScore)
                    {
                        Exception e = DAVectorUtility.SALSAError("Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " in Process " + DAVectorUtility.MPI_Rank.ToString()
                            + " Incorrect Score " + BestScore.ToString("F4") + " " + SmallValues[0].ToString("F4"));
                        throw (e);
                    }
                    if (BestCenter < 0)
                    {
                            Exception e = DAVectorUtility.SALSAError("Illegal Best Cluster in Initial Analysis " + BestCenter.ToString() + " Cluster Count " + Kmeans.Ncent_Global.ToString()
                                + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                            throw (e);
                    }
                    if (ActualNumber > KmeansTriangleInequality.Ncent_Global)
                    {
                        Exception e = DAVectorUtility.SALSAError("Too many centers in Initial Analysis " + ActualNumber.ToString() + " " + KmeansTriangleInequality.Ncent_Global.ToString()
                            + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                        throw (e);
                    }
                    
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].LimitforRest = BreakLimit;
                    for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                    {
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                    }

                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                    {
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest = BreakLimit;
                        for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                        {
                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                        }
                    }

                    KmeansTriangleInequality.NearestCentertoPoint[alpha] = BestCenter;
                    KmeansTriangleInequality.Distance_NearestCentertoPoint[alpha] = BestScore;

                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double)ActualNumber, 4);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BestScore, 5);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BreakLimit - BestScore, 25);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, SmallValues[ActualNumber-1], 6);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, KmeansTriangleInequality.Ncent_Global, 23);

                }   // End loop over Points
            }); // End Sum over Threads

        }   // End InitialPointAnalysis()

        public static void InitialPointAnalysis_FullOld()
        {
            //  Set Associated Cluster and initial lower bounds before any previous point data available
            //  This version cycles through all centers and does not try to be quick
            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 1);

                    int Num_LBValuesforPoint = KmeansTriangleInequality.MaxClusterLBsperPoint;    // Max Size of ordered list preserved in LB_Current_alpha_ClusterPointer
                    int ActualNumber = 0;   // For creating ordered list preserved in LB_Current_alpha_ClusterPointer
                    int currentcut = -1;    // For creating ordered list preserved in LB_Current_alpha_ClusterPointer
                    double Minrejected = -1.0;  // For creating ordered list preserved in LB_Current_alpha_ClusterPointer

                    int ActuallyUse = KmeansTriangleInequality.Ncent_Global;
                    double[] SmallValues = new double[ActuallyUse];    // For creating ordered list preserved in LB_Current/Old_alpha_ClusterPointer
                    int[] SmallIndices = new int[ActuallyUse]; // For creating ordered list preserved in LB_Current/Old_alpha_ClusterPointer

                    for(int CenterIndex =0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    {
                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterOld[CenterIndex], KmeansTriangleInequality.PointPosition[alpha]);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 15);

                        KmeansTriangleInequality.FindMinimumSetwithRemainder(PointCenterDistce, CenterIndex, ref currentcut, SmallValues, SmallIndices, ActuallyUse,
                            ref ActualNumber, ref Minrejected);

                    }   // End while over centers
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 7);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 0.0 , 19);

                    Array.Sort(SmallValues, SmallIndices, 0, ActualNumber);
                    if (Num_LBValuesforPoint < ActualNumber)
                    {
                        Minrejected = SmallValues[Num_LBValuesforPoint];
                        ActualNumber = Num_LBValuesforPoint;
                    }
                    if (ActualNumber > KmeansTriangleInequality.Ncent_Global)
                    {
                        Exception e = DAVectorUtility.SALSAError("Too many centers in Initial Analysis " + ActualNumber.ToString() + " " + KmeansTriangleInequality.Ncent_Global.ToString()
                            + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                        throw (e);
                    }
                    int BestCenter = SmallIndices[0];
                    double BestScore = SmallValues[0];
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].LimitforRest = Minrejected;
                    for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                    {
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                    }

                    KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                    KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                    KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                    KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest = Minrejected;
                    for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                    {
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                        KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                    }

                    KmeansTriangleInequality.NearestCentertoPoint[alpha] = BestCenter;
                    KmeansTriangleInequality.Distance_NearestCentertoPoint[alpha] = BestScore;

                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double)ActualNumber, 4);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BestScore, 5);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 0.0, 25);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, SmallValues[ActualNumber - 1], 6);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, KmeansTriangleInequality.Ncent_Global, 23);

                }   // End loop over Points
            }); // End Sum over Threads

        }   // End InitialPointAnalysis_FullOld()

        public static void PointDependentAnalysis()
        {
            if ( KmeansTriangleInequality.RecalculateOldCenters )
            {
                InitialPointAnalysis_FullOld();
                return;
            }
            
            //  Set Associated Cluster and reset lower bounds following previous iterations
            //  It loops through centers in order of centers near its previous center
            //  When a new best center is found, it changes to loop over its associated centers
            //  If best center loop exhausted it reverts to ploughing through center list

            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 1);

                    // Set up Look up arrays
                    bool[] LookedAt = new bool[KmeansTriangleInequality.Ncent_Global];
                    int[] LastBoundLookup = new int[KmeansTriangleInequality.Ncent_Global];
                    int Oldsize = 1;
                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                        Oldsize = KmeansTriangleInequality.Ncent_Global;
                    int[] OldBoundLookup = new int[Oldsize];

                    for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    {
                        LookedAt[CenterIndex] = false;
                        if (KmeansTriangleInequality.OldCenterOption >= 0)
                            OldBoundLookup[CenterIndex] = -1;
                        LastBoundLookup[CenterIndex] = -1;
                    }

                    //  Arrays for finding minimum set
                    int Num_LBValuesforPoint = KmeansTriangleInequality.MaxClusterLBsperPoint;
                    int ActuallyUse = KmeansTriangleInequality.Ncent_Global;
                    double[] SmallValues = new double[ActuallyUse];
                    int[] SmallIndices = new int[ActuallyUse];

                    //  Set up parameters of Lower Bound Finding for use by old
                    int ActualNumber = 0;
                    int currentcut = -1;
                    double Minrejected = -1.0;

                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                    {
                        for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                            OldBoundLookup[KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex]] = loopindex;
                    }
                    for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                        LastBoundLookup[KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex]] = loopindex;

                    //  Update bounds for recalculated centers in both Old and Last arrays
                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                    {
                        bool changed = false;
                        for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                        {
                            int CenterIndex = KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                            if (KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0)
                            {
                                changed = true;
                                double lowerbound = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.PointPosition[alpha]);
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 11);
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                                KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = lowerbound;
                                if (LastBoundLookup[CenterIndex] >= 0)
                                    KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[LastBoundLookup[CenterIndex]] = lowerbound;
                            }
                            KmeansTriangleInequality.FindMinimumSetwithRemainder(KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex], CenterIndex, ref currentcut, SmallValues, SmallIndices, ActuallyUse,
                                ref ActualNumber, ref Minrejected);
                        }
                        //  Add in Updated Centers NOT in magic set
                        for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                        {
                            if (OldBoundLookup[CenterIndex] >= 0)
                                continue;
                            if (KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0)
                            {
                                changed = true;
                                double lowerbound = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.PointPosition[alpha]);
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 11);
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                                KmeansTriangleInequality.FindMinimumSetwithRemainder(lowerbound, CenterIndex, ref currentcut, SmallValues, SmallIndices, ActuallyUse,
                                ref ActualNumber, ref Minrejected);
                            }
                        }
                        
                        if (changed)
                        {
                            //  Set has all centers that were in original list plus those that were changed -- whether or not in list

                            Array.Sort(SmallValues, SmallIndices, 0, ActualNumber);
                            if (Num_LBValuesforPoint < ActualNumber)
                            {
                                Minrejected = SmallValues[Num_LBValuesforPoint];
                                ActualNumber = Num_LBValuesforPoint;
                            }
                            else
                                Minrejected = KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest;
                            Minrejected = Math.Min(Minrejected, KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest);

                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].AssociatedClusterIndex = SmallIndices[0];
                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NearestDistance = SmallValues[0];
                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                            KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest = Minrejected;
                            for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                            {
                                KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                                KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                            }
                        }

                    }   // End Analysis of Old Centers

                    for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                    {
                        double lowerbound;
                        int CenterIndex = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                        if (KmeansTriangleInequality.OldCenterOption >= 0)
                        {
                            if (OldBoundLookup[CenterIndex] >= 0)
                                continue;
                        }
                        if (KmeansTriangleInequality.CurrentInterCenter[CenterIndex].RecalculateStatus > 0)
                        {
                            lowerbound = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], KmeansTriangleInequality.PointPosition[alpha]);
                            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 12);
                            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                            KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = lowerbound;
                        }
                    }

                    //  Reset OldBoundLookup
                    if (KmeansTriangleInequality.OldCenterOption >= 0)
                    {
                        for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                                OldBoundLookup[CenterIndex] = -1;
                        for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                            OldBoundLookup[KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].CenterIndices[loopindex]] = loopindex;
                    }
                    //  Reset up parameters of Lower Bound Finding for use by current
                    ActualNumber = 0;
                    currentcut = -1;
                    Minrejected = -1.0;

                    //  Set up Scan for best center
                    int NumCentersToGo = KmeansTriangleInequality.Ncent_Global;  // Number of Centers still to look at -- end when <= 0

                    int CurrentCandidateCenter = -1;    // Current Center being looked at in Method = -1 approach

                    int Method = -2;    // Method = -2 Start with Previous best, Method = -1 scan List in order of centers; Method > 0 Scan centers in list of center Method
                    int Bestposition = -2;
                    int MethodPosition = -1; // Position in Center List of Center Method

                    int BestCenter = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].AssociatedClusterIndex;    // Previous nearest center

                    double BestScore = -1.0;    // Distance of nearest center
                    int ExaminedCenter = BestCenter;    // Current center being looked at
                    int BreakReason = 7;    // Reason for Break; 7 is continue

                    double BreakLimit = -1.0;

                    //  Loop over candidates
                    while (NumCentersToGo > 0)
                    {
                        if( Method == -2)
                        {
                            Bestposition++;
                            if (Bestposition == -1)
                            {
                                ExaminedCenter = BestCenter;
                                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 16);
                            }
                            else
                            {
                                if( Bestposition < KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters)
                                {
                                    ExaminedCenter = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[Bestposition];
                                    if (LookedAt[ExaminedCenter])
                                        continue;
                                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 16);
                                }
                                else
                                    Method = -1;
                            }
                        }
                        if (Method == -1)
                        {
                            CurrentCandidateCenter++;
                            if (LookedAt[CurrentCandidateCenter])
                                continue;
                            ExaminedCenter = CurrentCandidateCenter;
                            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 17);
                        }
                        else if (Method >= 0)
                        {   // Examining centers near current best center in position "Method"
                            ++MethodPosition;
                            if (MethodPosition >= KmeansTriangleInequality.CurrentInterCenter[Method].NearbyClusters)
                            {   // No more centers left in center list
                                if (KmeansTriangleInequality.CurrentInterCenter[Method].Dmax > 2.0 * BestScore)
                                {
                                    BreakReason = 8;
                                    if (BreakLimit < -0.5)
                                        BreakLimit = KmeansTriangleInequality.CurrentInterCenter[Method].Dmax - BestScore;
                                    else
                                        BreakLimit = Math.Min(BreakLimit, KmeansTriangleInequality.CurrentInterCenter[Method].Dmax - BestScore);
                                    
                                    break;
                                }
                                Method = -2;    // Go back to linear scan
                                continue;
                            }
                            if (KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] >= 2.0 * BestScore)
                            {
                                BreakReason = 9;
                                if (BreakLimit < -0.5)
                                    BreakLimit = KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] - BestScore;
                                else
                                    BreakLimit = Math.Min(BreakLimit, KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifference[MethodPosition] - BestScore);
                                break;
                            }
                            ExaminedCenter = KmeansTriangleInequality.CurrentInterCenter[Method].CenterDifferenceIndices[MethodPosition];
                            if (LookedAt[ExaminedCenter])
                                continue;
                            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 18);
                        }

                        //  See if bound exists
                        if( BestScore > - 0.5 )
                        {
                            int status_last = LastBoundLookup[ExaminedCenter];
                            if (status_last >= 0)
                            {
                                if ((KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[status_last]
                                - CurrentInterCenter[ExaminedCenter].LastCenterChange) > BestScore)
                                {
                                    --NumCentersToGo;
                                    LookedAt[ExaminedCenter] = true;
                                    LastBoundLookup[ExaminedCenter] = -2;
                                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 3);
                                    continue;
                                }
                                else
                                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 14);
                            }

                            else
                            {
                                double TestLimit = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].LimitforRest - CurrentInterCenter[ExaminedCenter].LastCenterChange;
                                if (TestLimit > BestScore)
                                {
                                    --NumCentersToGo;
                                    LookedAt[ExaminedCenter] = true;
                                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 24);
                                    if (BreakLimit < -0.5)
                                        BreakLimit = TestLimit;
                                    else
                                        BreakLimit = Math.Min(BreakLimit, TestLimit);
                                    continue;
                                }
                            }

                            if (KmeansTriangleInequality.OldCenterOption >= 0)
                            {
                                int status_old = OldBoundLookup[ExaminedCenter];
                                if (status_old >= 0)
                                {
                                    double TestLimit = KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].DistceLowerBound[status_old]
                                        - CurrentInterCenter[ExaminedCenter].OldCenterChange;
                                    if ( TestLimit > BestScore)
                                    {
                                        --NumCentersToGo;
                                        LookedAt[ExaminedCenter] = true;
                                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 2);
                                        if (BreakLimit < -0.5)
                                            BreakLimit = TestLimit;
                                        else
                                            BreakLimit = Math.Min(BreakLimit, TestLimit);
                                        continue;
                                    }
                                    else
                                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 13);
                                }
                                else
                                {
                                    double TestLimit = KmeansTriangleInequality.LB_Old_alpha_ClusterPointer[alpha].LimitforRest - CurrentInterCenter[ExaminedCenter].OldCenterChange;
                                    if ( TestLimit > BestScore)
                                    {
                                        --NumCentersToGo;
                                        LookedAt[ExaminedCenter] = true;
                                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 26);
                                        if (BreakLimit < -0.5)
                                            BreakLimit = TestLimit;
                                        else
                                            BreakLimit = Math.Min(BreakLimit, TestLimit);
                                        continue;
                                    }
                                }
                            }
                        }

                        //  Calculate Explicit Distance
                        if (LastBoundLookup[ExaminedCenter] >= 0)
                            LastBoundLookup[ExaminedCenter] = -1;

                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterCurrent[ExaminedCenter], KmeansTriangleInequality.PointPosition[alpha]);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0 , 15);

                        KmeansTriangleInequality.FindMinimumSetwithRemainder(PointCenterDistce, ExaminedCenter, ref currentcut, SmallValues, SmallIndices, ActuallyUse, 
                            ref ActualNumber, ref Minrejected);

                        if( (BestScore > -0.5) && (PointCenterDistce > BestScore))
                        {
                            if (Method < 0)
                            {
                                Method = ExaminedCenter;
                                MethodPosition = -1;
                            }
                            --NumCentersToGo;
                            LookedAt[ExaminedCenter] = true;
                            if( KmeansTriangleInequality.DoBackwardFacingTests)
                                KmeansTriangleInequality.BackwardTest(ThreadIndex, KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter], ref NumCentersToGo, LookedAt, PointCenterDistce, BestScore, BreakLimit);
                            continue;
                        }

                        //  New best Center discovered
                        if(BestScore > -0.5)
                        {
                            if( KmeansTriangleInequality.DoBackwardFacingTests )
                                KmeansTriangleInequality.BackwardTest(ThreadIndex, KmeansTriangleInequality.CurrentInterCenter[BestCenter], ref NumCentersToGo, LookedAt, PointCenterDistce, BestScore, BreakLimit);
                        }

                        BestScore = PointCenterDistce;
                        BestCenter = ExaminedCenter;
                        if (KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] >= 2.0 * PointCenterDistce)
                        {
                            BreakReason = 10;
                            if (BreakLimit < -0.5)
                                BreakLimit = KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] - BestScore;
                            else
                                BreakLimit = Math.Min(BreakLimit, KmeansTriangleInequality.CurrentInterCenter[ExaminedCenter].CenterDifference[0] - BestScore);
                            break;
                        }
                        Method = ExaminedCenter;
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 22);
                        MethodPosition = -1;
                        --NumCentersToGo;
                        LookedAt[ExaminedCenter] = true;
                        continue;

                    }   // End while over centers

                    if (BestCenter < 0)
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal Best Cluster in Ongoing Analysis " + BestCenter.ToString() + " Cluster Count " + Kmeans.Ncent_Global.ToString()
                            + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                        throw (e);
                    }

                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, BreakReason);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double) NumCentersToGo, 19);

                    //  Combine Last and Current Lower Bounds (stored above in Small values/indices and store back into Current
                    BreakLimit = Math.Max(BreakLimit, BestScore);

                    for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                    {
                        int CenterIndex = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                        if (LastBoundLookup[CenterIndex] == -1)
                            continue;   // This center is aleady in SmallValues/Indices
                        double lowerbound = KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] - CurrentInterCenter[CenterIndex].LastCenterChange;
                        lowerbound = Math.Max(BestScore, lowerbound);
                        if (LastBoundLookup[CenterIndex] >= 0)
                        {   // This lower bound was not used. So at least BreakLimit
                            if(lowerbound < BreakLimit)
                                continue;
                        }
                        KmeansTriangleInequality.FindMinimumSetwithRemainder(lowerbound, CenterIndex, ref currentcut, SmallValues, SmallIndices, ActuallyUse,
                            ref ActualNumber, ref Minrejected);
                    }

                    if (ActualNumber > KmeansTriangleInequality.Ncent_Global)
                    {
                        Exception e = DAVectorUtility.SALSAError("Too many centers in Ongoing Analysis " + ActualNumber.ToString() + " " + KmeansTriangleInequality.Ncent_Global.ToString()
                            + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Rank " + DAVectorUtility.MPI_Rank.ToString());
                        throw (e);
                    }
  
                    Array.Sort(SmallValues, SmallIndices, 0, ActualNumber);
                    if (Num_LBValuesforPoint < ActualNumber)
                    {
                        Minrejected = SmallValues[Num_LBValuesforPoint];
                        BreakLimit = Math.Min(BreakLimit, Minrejected);

                        ActualNumber = Num_LBValuesforPoint;
                    }

                    //  Store LB_Current_alpha_ClusterPointer
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].AssociatedClusterIndex = BestCenter;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NearestDistance = BestScore;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].NumCenters = ActualNumber;
                    KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].LimitforRest = BreakLimit;
                    for (int loopindex = 0; loopindex < ActualNumber; loopindex++)
                    {
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex] = SmallValues[loopindex];
                        KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex] = SmallIndices[loopindex];
                    }

                    KmeansTriangleInequality.NearestCentertoPoint[alpha] = BestCenter;
                    KmeansTriangleInequality.Distance_NearestCentertoPoint[alpha] = BestScore;
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double)ActualNumber, 4);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BestScore, 5);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BreakLimit - BestScore, 25);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, SmallValues[ActualNumber-1], 6);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double) KmeansTriangleInequality.Ncent_Global, 23);

                }   // End loop over Points

            }); // End Sum over Threads

        }   // End PointDependentAnalysis()

        public static void BackwardTest(int ThreadIndex, CenterData CurrentCenterData, ref int NumberCentersLeft, bool[] LookedAtArray, double Distance, double BestScore, double BreakLimit)
        {   // Test backwards in Center-Center Array
            KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 21);

            for (int CenterIndexPointer = CurrentCenterData.NearbyClusters - 1; CenterIndexPointer >= 0; CenterIndexPointer--)
            {
                double TestLimit = CurrentCenterData.CenterDifference[CenterIndexPointer] - Distance;
                if( TestLimit < BestScore )
                    return;
                int CenterIndex = CurrentCenterData.CenterDifferenceIndices[CenterIndexPointer];
                if (LookedAtArray[CenterIndex])
                    continue;
                LookedAtArray[CenterIndex] = true;
                --NumberCentersLeft;
                
                if (BreakLimit < -0.5)
                    BreakLimit = TestLimit;
                else
                    BreakLimit = Math.Min(BreakLimit, TestLimit);

                KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 20);
            }
        }   // End BackwardTest


        public static void PlotCenterHistograms(double PlotPeakSize)
        {
            double Histmin = 0.0;
            double Histmax = 10.0 * PlotPeakSize;
            DAVectorUtility.SetUpHistogramRange(KmeansTriangleInequality.HistogramSize - 2, ref Histmin, ref Histmax);

            if( !KmeansTriangleInequality.UseParallelismoverCenters)
                return;

            GlobalReductions.FindVectorDoubleSum FindInCCBounds = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindActualCCBounds = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindOutCCBounds = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);

            //  Parallel Center Processing
            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int NumberCentersToProcess = KmeansTriangleInequality.FullParallel_CentersperThread[ThreadIndex];
                int BeginCenterThread = KmeansTriangleInequality.FullParallel_StartCenterperThread[ThreadIndex];
                for (int CenterIndex = BeginCenterThread; CenterIndex < NumberCentersToProcess + BeginCenterThread; CenterIndex++)
                {
                    bool[] LookedAt = new bool[KmeansTriangleInequality.Ncent_Global];
                    for (int CenterIndexPrime = 0; CenterIndexPrime < KmeansTriangleInequality.Ncent_Global; CenterIndexPrime++)
                        LookedAt[CenterIndexPrime] = false;
                    int NumberPrimedCenters = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].NearbyClusters;
                    for (int CenterPrimeLoop = 0; CenterPrimeLoop < NumberPrimedCenters; CenterPrimeLoop++)
                    {
                        int CenterIndexPrime = KmeansTriangleInequality.CurrentInterCenter[CenterIndex].CenterDifferenceIndices[CenterPrimeLoop];
                        if (CenterIndexPrime == CenterIndex)
                            continue;
                        if (LookedAt[CenterIndexPrime])
                            continue;
                        LookedAt[CenterIndexPrime] = true;
                        double CenterCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], 
                            KmeansTriangleInequality.CenterCurrent[CenterIndexPrime]);
                        int HistPos = GetHistogramPosition(CenterCenterDistce, Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        FindInCCBounds.addapoint(ThreadIndex, 1.0, HistPos);
                        HistPos = GetHistogramPosition(KmeansTriangleInequality.CurrentInterCenter[CenterIndex].CenterDifference[CenterPrimeLoop], 
                            Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        FindActualCCBounds.addapoint(ThreadIndex, 1.0, HistPos);
                    }
                    for (int CenterIndexPrime = 0; CenterIndexPrime < KmeansTriangleInequality.Ncent_Global; CenterIndexPrime++)
                    {
                        if (CenterIndexPrime == CenterIndex)
                            continue;
                        if (LookedAt[CenterIndexPrime])
                            continue;
                        double CenterCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex],
                            KmeansTriangleInequality.CenterCurrent[CenterIndexPrime]);
                        int HistPos = GetHistogramPosition(CenterCenterDistce, Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        FindOutCCBounds.addapoint(ThreadIndex, 1.0, HistPos);
                    }
                }

            }); // End Loop over Threads

            FindInCCBounds.sumoverthreadsandmpi();
            FindActualCCBounds.sumoverthreadsandmpi();
            FindOutCCBounds.sumoverthreadsandmpi();

            double[] DistcesinCCBounds = new double[KmeansTriangleInequality.HistogramSize];
            double[] ActualCCBounds = new double[KmeansTriangleInequality.HistogramSize];
            double[] DistcesNotinCCBounds = new double[KmeansTriangleInequality.HistogramSize];
            for (int HistLoop = 0; HistLoop < KmeansTriangleInequality.HistogramSize; HistLoop++)
            {
                DistcesinCCBounds[HistLoop] = FindInCCBounds.TotalVectorSum[HistLoop];
                ActualCCBounds[HistLoop] = FindActualCCBounds.TotalVectorSum[HistLoop];
                DistcesNotinCCBounds[HistLoop] = FindOutCCBounds.TotalVectorSum[HistLoop];
            }

            DAVectorUtility.SALSAPrint(0, "\nCenter Histograms with Minimum " + Histmin.ToString("F3") + " Maximum " + Histmax.ToString("E4"));
            KmeansTriangleInequality.PrintPointHistogram("In Center-Center Bounds    ", DistcesinCCBounds, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Actual Center-Center Bounds", ActualCCBounds, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Not in Center-Center Bounds", DistcesNotinCCBounds, KmeansTriangleInequality.HistogramSize);

        }   // End PlotCenterHistograms

        public static void PlotPointHistograms(double PlotPeakSize)
        {
            double Histmin = 0.0;
            double Histmax = 10.0 * PlotPeakSize;
            DAVectorUtility.SetUpHistogramRange(KmeansTriangleInequality.HistogramSize - 2, ref Histmin, ref Histmax);

            GlobalReductions.FindVectorDoubleSum FindInsideCenter = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindOutsideCenterbutinLB = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindOutsideCenterActualLB = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindCenterActualLB = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);
            GlobalReductions.FindVectorDoubleSum FindOutofLB = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, KmeansTriangleInequality.HistogramSize);

            //   Loop over points
            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    // Set up Look up array
                    bool[] LookedAt = new bool[KmeansTriangleInequality.Ncent_Global];
                    for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                        LookedAt[CenterIndex] = false;

                    for (int loopindex = 0; loopindex < KmeansTriangleInequality.LB_Last_alpha_ClusterPointer[alpha].NumCenters; loopindex++)
                    {
                        int CenterIndex = KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].CenterIndices[loopindex];
                        if (LookedAt[CenterIndex])
                            continue;
                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex], 
                            KmeansTriangleInequality.PointPosition[alpha]);
                        int HistPos = GetHistogramPosition(PointCenterDistce, Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        if (CenterIndex == KmeansTriangleInequality.NearestCentertoPoint[alpha])
                            FindInsideCenter.addapoint(ThreadIndex, 1.0, HistPos);
                        else
                            FindOutsideCenterbutinLB.addapoint(ThreadIndex, 1.0, HistPos);
                        HistPos = GetHistogramPosition(KmeansTriangleInequality.LB_Current_alpha_ClusterPointer[alpha].DistceLowerBound[loopindex], 
                            Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        if (CenterIndex == KmeansTriangleInequality.NearestCentertoPoint[alpha])
                            FindCenterActualLB.addapoint(ThreadIndex, 1.0, HistPos);
                        else
                            FindOutsideCenterActualLB.addapoint(ThreadIndex, 1.0, HistPos);
                        LookedAt[CenterIndex] = true;
                    }
                    for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    {
                        if (LookedAt[CenterIndex])
                            continue;
                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(KmeansTriangleInequality.CenterCurrent[CenterIndex],
                            KmeansTriangleInequality.PointPosition[alpha]);
                        int HistPos = GetHistogramPosition(PointCenterDistce, Histmin, Histmax, KmeansTriangleInequality.HistogramSize);
                        if (CenterIndex == KmeansTriangleInequality.NearestCentertoPoint[alpha])
                            FindInsideCenter.addapoint(ThreadIndex, 1.0, HistPos);
                        else
                            FindOutofLB.addapoint(ThreadIndex, 1.0, HistPos);
                    }

                }   // End loop over Points
            }); // End Sum over Threads

            FindInsideCenter.sumoverthreadsandmpi();
            FindOutsideCenterbutinLB.sumoverthreadsandmpi();
            FindOutsideCenterActualLB.sumoverthreadsandmpi();
            FindCenterActualLB.sumoverthreadsandmpi();
            FindOutofLB.sumoverthreadsandmpi();

            double[] PointsInsideCenter = new double[KmeansTriangleInequality.HistogramSize];
            double[] PointsOutsideCenterbutinLB = new double[KmeansTriangleInequality.HistogramSize];
            double[] PointsOutsideCenterActualLB = new double[KmeansTriangleInequality.HistogramSize];
            double[] PointsCenterActualLB = new double[KmeansTriangleInequality.HistogramSize];
            double[] PointsOutofLB = new double[KmeansTriangleInequality.HistogramSize];
            for (int HistLoop = 0; HistLoop < KmeansTriangleInequality.HistogramSize; HistLoop++)
            {
                PointsInsideCenter[HistLoop] = FindInsideCenter.TotalVectorSum[HistLoop];
                PointsOutsideCenterbutinLB[HistLoop] = FindOutsideCenterbutinLB.TotalVectorSum[HistLoop];
                PointsOutsideCenterActualLB[HistLoop] = FindOutsideCenterActualLB.TotalVectorSum[HistLoop];
                PointsCenterActualLB[HistLoop] = FindCenterActualLB.TotalVectorSum[HistLoop];
                PointsOutofLB[HistLoop] = FindOutofLB.TotalVectorSum[HistLoop];
            }

            DAVectorUtility.SALSAPrint(0, "\nPoint-Center Histograms with Minimum " + Histmin.ToString("F3") + " Maximum " + Histmax.ToString("E4"));
            KmeansTriangleInequality.PrintPointHistogram("Inside Center       ", PointsInsideCenter, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Inside Center LB    ", PointsCenterActualLB, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Outside Center in LB", PointsOutsideCenterbutinLB, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Not Center Actual LB", PointsOutsideCenterActualLB, KmeansTriangleInequality.HistogramSize);
            KmeansTriangleInequality.PrintPointHistogram("Outside Lower Bounds", PointsOutofLB, KmeansTriangleInequality.HistogramSize);

        }   //  End PlotPointHistograms

        public static int GetHistogramPosition( double value, double Histmin, double Histmax, int HistSize)
        {
            int HistPos = 1 + (int)Math.Floor( ((double) (HistSize - 2)) * (value - Histmin) / (Histmax - Histmin) );
            if(HistPos < 0 )
                HistPos = 0;
            HistPos = Math.Min( HistPos, HistSize - 1);
            return HistPos;

        }   // End GetHistogramPosition( double value, double Histmin, double Histmax, int HistSize)

        public static void PrintPointHistogram(string Label, double[] HistCounts, int NumberBinsplus2)
        {
            if (DAVectorUtility.MPI_Rank != 0)
                return;

            double Total = 0.0;
            for (int HistLoop = 0; HistLoop < NumberBinsplus2; HistLoop++)
                Total += HistCounts[HistLoop];
            if (Total < 0.5)
                return;

            for (int HistLoop = 0; HistLoop < NumberBinsplus2; HistLoop++)
                HistCounts[HistLoop] *= (100.0 / Total);
            int ilabel = 0;
            string message = Label + " Total " + Total.ToString("E4");
            for (int HistLoop = 0; HistLoop < NumberBinsplus2; HistLoop++)
            {
                if (HistLoop % 10 == 1)
                {
                    message += " *" + ilabel.ToString() + "*";
                    ++ilabel;
                }
                if (HistLoop == 0)
                    message += " Underflow";
                if (HistLoop == (NumberBinsplus2 - 1))
                    message += " Overflow";
                message += " " + HistCounts[HistLoop].ToString("F3") + ","; ;
            }
            DAVectorUtility.SALSAPrint(0, "\n" + message);

        }   // End PrintPointHistogram

        public static void PrintDiagnostics()
        {
            DAVectorUtility.SALSAPrint(0, "\nTriangle Inequality Diagnostics ***** Number of Iterations " + KmeansTriangleInequality.CalculationIterations.ToString());

            DAVectorUtility.SALSAPrint(0, "Triangle Inequality Option " + UseTriangleInequality.ToString()  + " Old Center " +  KmeansTriangleInequality.OldCenterOption.ToString() 
                + " Use Backward Tests? " + KmeansTriangleInequality.DoBackwardFacingTests.ToString() );
            DAVectorUtility.SALSAPrint(0, "Test for center change and old lower bounds (normalized by radius) " + KmeansTriangleInequality.TriangleInequality_Delta1_old.ToString("F4") );
            DAVectorUtility.SALSAPrint(0, "Test for center change and current lower bounds (normalized by radius) " + KmeansTriangleInequality.TriangleInequality_Delta1_current.ToString("F4") );
            DAVectorUtility.SALSAPrint(0, "Maximum Allowed # Centers " + KmeansTriangleInequality.MaxNcent_Global.ToString() + " Cut for Center Parallelism "
                + KmeansTriangleInequality.Ncent_Global_Parallel.ToString() + " Parallelism used? " + KmeansTriangleInequality.UseParallelismoverCenters);
            DAVectorUtility.SALSAPrint(0, "Vector Dimension " + KmeansTriangleInequality.ParameterVectorDimension.ToString());
            DAVectorUtility.SALSAPrint(0, "Maximum number of Lower Bound values for each Point " + KmeansTriangleInequality.MaxClusterLBsperPoint.ToString());
            DAVectorUtility.SALSAPrint(0, "Maximum number of Centers in Center Difference Array " + KmeansTriangleInequality.MaxCentersperCenter.ToString());

            double tmp1 = 1.0;
            if (KmeansTriangleInequality.CalculationIterations > 0 )
                tmp1 = 1.0 / (double)KmeansTriangleInequality.CalculationIterations;
            double tmp2 = 1.0;
            if (KmeansTriangleInequality.CalculationPointIterations > 0)
                tmp2 = 1.0 / KmeansTriangleInequality.CalculationPointIterations;
            double tmp3 = 1.0;
            if (KmeansTriangleInequality.CalculationCenterIterations > 0)
                tmp3 = 1.0 / KmeansTriangleInequality.CalculationCenterIterations;
            double tmp4 = 1.0;
            if (KmeansTriangleInequality.CountCCEntries > 0)
                tmp4 = 1.0 / KmeansTriangleInequality.CountCCEntries;

            DAVectorUtility.SALSAPrint(0, "Points per Iteration " + (KmeansTriangleInequality.CalculationPointIterations * tmp1).ToString("F2"));
            DAVectorUtility.SALSAPrint(0, "Centers per Iteration " + (KmeansTriangleInequality.CalculationCenterIterations * tmp1).ToString("F2"));

            DAVectorUtility.SALSAPrint(0, "Center-Center Distances per Center-Iteration " + (KmeansTriangleInequality.NumberFullDistancesCalculatedCC * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Center Moving Too much per Iteration " + (KmeansTriangleInequality.CountToobigAShift * tmp1).ToString("F6"));

            DAVectorUtility.SALSAPrint(0, "\nAverage Center Radius " + (KmeansTriangleInequality.AccumulateRadius * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Largest Center Center Distance NOT Recorded (-1 if all done) " + (KmeansTriangleInequality.AccumulateDmax * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of recorded Center-Center distances " + (KmeansTriangleInequality.CountCCEntries * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average value  of recorded Center-Center distances " + (KmeansTriangleInequality.AccumulateCCvalues * tmp4).ToString("F6"));

            DAVectorUtility.SALSAPrint(0, "\nAverage number of OLD Lower Bounds tested Good " + (KmeansTriangleInequality.AccumulateOLD_LB * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of LAST Lower Bounds tested Good " + (KmeansTriangleInequality.AccumulateCURRENT_LB * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of LAST Lower Bound Limits used to Reject " + (KmeansTriangleInequality.AccumulateCURRENT_LBLimit * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of OLD Lower Bound Limits used to Reject " + (KmeansTriangleInequality.AccumulateOLD_LBLimit * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of OLD Lower Bounds tested Bad " + (KmeansTriangleInequality.AccumulateOLD_LBFail * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of LAST Lower Bounds tested Bad " + (KmeansTriangleInequality.AccumulateCURRENT_LBFail * tmp2).ToString("F6"));

            DAVectorUtility.SALSAPrint(0, "\nAverage number of CURRENT Lower Bound  distance entries " + (KmeansTriangleInequality.CountLBEntries * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of CURRENT Best Score " + (KmeansTriangleInequality.AccumulateBestScore * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of CURRENT Break Limit - Best Score " + (KmeansTriangleInequality.AccumulateBreakLimit * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average of Center change distances (Last) Changes  " + (KmeansTriangleInequality.AccumulateLastCenterchanges * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average of Center change distances (Old) Changes  " + (KmeansTriangleInequality.AccumulateOldCenterchanges * tmp3).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of CURRENT Maximum value of Lower Bound " + (KmeansTriangleInequality.AccumulateMAX_LB * tmp2).ToString("F6"));

            DAVectorUtility.SALSAPrint(0, "\nPoint-Centers Distances per Point-Iteration " + (KmeansTriangleInequality.PossibleNumberFullDistancesCalculated * tmp2).ToString("F6") + " Possible");
            DAVectorUtility.SALSAPrint(0, "\nPoint-Centers Distances per Point-Iteration " + (KmeansTriangleInequality.NumberFullDistancesCalculatedPC * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + (KmeansTriangleInequality.NumberFullDistancesCalculatedPC_Old * tmp2).ToString("F6") + " OLD refresh");
            DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + (KmeansTriangleInequality.NumberFullDistancesCalculatedPC_Last * tmp2).ToString("F6") + " LAST refresh");
            DAVectorUtility.SALSAPrint(0, "Point-Centers Distances per Point-Iteration " + (KmeansTriangleInequality.NumberFullDistancesCalculatedLoop * tmp2).ToString("F6") + " In while loop");

            DAVectorUtility.SALSAPrint(0, "\nAverage number of Breaks in LB/Center distance test as while limit hit " + (KmeansTriangleInequality.AccumulateBreakReason7 * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Dmax   " + (KmeansTriangleInequality.AccumulateBreakReason8 * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Middle " + (KmeansTriangleInequality.AccumulateBreakReason9 * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Average number of Breaks in LB/Center distance test as Center Scan Start  " + (KmeansTriangleInequality.AccumulateBreakReason10 * tmp2).ToString("F6"));
            DAVectorUtility.SALSAPrint(0, "Number of Centers NOT looked at " + (KmeansTriangleInequality.AverageCentersLookedAt * tmp2).ToString("F6") + " After while loop");
            DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  -2 (Best) " + (KmeansTriangleInequality.NumberCentersinBestSet * tmp2).ToString("F6") + " In while loop");
            DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  -1 (Linear) " + (KmeansTriangleInequality.NumberCentersinSimpleSet * tmp2).ToString("F6") + " In while loop");
            DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Method  >=0 (Center facing) " + (KmeansTriangleInequality.NumberCentersinCenterbasedset * tmp2).ToString("F6") + " In while loop");
            DAVectorUtility.SALSAPrint(0, "Number of such Center facing forward Tests " + (KmeansTriangleInequality.NumberForwardTests * tmp2).ToString("F6") + " In while loop"); 
            DAVectorUtility.SALSAPrint(0, "Number of Centers removed with Backward Test " + (KmeansTriangleInequality.NumberCentersRemovedBackwards * tmp2).ToString("F6") + " In while loop");
            DAVectorUtility.SALSAPrint(0, "Number of Center facing Backward Tests " + (KmeansTriangleInequality.NumberBackwardTests * tmp2).ToString("F6") + " In while loop"); 

        }   // End PrintDiagnostics()

        public static void PureKmeans()
        {
            //  Set Associated Cluster and reset lower bounds following previous iterations
            //  It loops through centers in order of centers near its previous center
            //  When a new best center is found, it changes to loop over its associated centers
            //  If best center loop exhausted it reverts to ploughing through center list
            Parallel.For(0, KmeansTriangleInequality._parallelOptions.MaxDegreeOfParallelism, KmeansTriangleInequality._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    // Set up Look up arrays
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 1);
                    int BestCenter = -1;
                    double BestScore = 0.0;
                    for (int CenterIndex = 0; CenterIndex < KmeansTriangleInequality.Ncent_Global; CenterIndex++)
                    {
                        double PointCenterDistce = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterCurrent[CenterIndex], KmeansTriangleInequality.PointPosition[alpha]);
                        KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, 1.0, 0);
                        if ((BestCenter < 0) || (PointCenterDistce < BestScore))
                        {
                            BestScore = PointCenterDistce;
                            BestCenter = CenterIndex;
                        }

                    }   // End loop over centers
                    KmeansTriangleInequality.NearestCentertoPoint[alpha] = BestCenter;
                    KmeansTriangleInequality.Distance_NearestCentertoPoint[alpha] = BestScore;
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, (double)KmeansTriangleInequality.Ncent_Global, 4);
                    KmeansTriangleInequality.FindDiagnosticSums_Points.addapoint(ThreadIndex, BestScore, 5);

                }   // End loop over Points
            }); // End Sum over Threads

        }   // End PureKmeans()

        //  Support finding list of minimum values by inserting new point with value newvalue and index newindex into lists SmallValues SmallIndices
        //  In SmallIndices negative values correspond to unset values
        //  NumberSmallOnes is total number wanted
        //  Currentcut is position in SmallValues, SmallIndices of largest min value
        //  Minrejected is Minimum of values not included
        //  Minrejected must be set to -1.0 and currentcut to -1
        public static void FindMinimumSetwithRemainder(double newvalue, int newindex, ref int currentcut, double[] SmallValues, int[] SmallIndices, int NumberSmallones,
            ref int ActualNumber, ref double Minrejected)
        {
            if (currentcut < 0)
            {
                currentcut = 0;
                SmallValues[0] = newvalue;
                SmallIndices[0] = newindex;
                ActualNumber = 1;
                return;
            }
            if (ActualNumber < NumberSmallones)
            {   // Not all positions are filled so add at next available
                // Reset currentcut if worst
                SmallValues[ActualNumber] = newvalue;
                SmallIndices[ActualNumber] = newindex;
                if (SmallValues[ActualNumber] > SmallValues[currentcut])
                    currentcut = ActualNumber;
                ++ActualNumber;
                return;
            }
            if (newvalue >= SmallValues[currentcut])
            {
                if (Minrejected < 0.0)
                {
                    Minrejected = newvalue;
                    return;
                }
                Minrejected = Math.Min(Minrejected, newvalue);
                return;
            }

            // Replace currentcut position with new values and Reset new worst position
            Minrejected = SmallValues[currentcut];
            SmallValues[currentcut] = newvalue;
            SmallIndices[currentcut] = newindex;
            double maxvalue = -1.0;
            for (int ivalue = 0; ivalue < NumberSmallones; ivalue++)
            {
                if (SmallIndices[ivalue] < 0)
                    continue;
                if (SmallValues[ivalue] > maxvalue)
                {
                    currentcut = ivalue;
                    maxvalue = SmallValues[ivalue];
                }
            }
            return;

        }   // End FindMinimumSetwithRemainder


    }   // end class TriangleInequality

}   // End Namespace Salsa.DAVectorSponge