using System;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using MPI;
using Salsa.Core;
using SALSALibrary;

namespace Salsa.DAVectorSponge
{

    // For fast look up of cluster locations
    //  IterationSet allows one to designate which information is really up to date
    public class ClusterIndirection
    {
        public int IterationSet;   // Iteration on which this cluster information set
        public int Availability;       // >0 1 + Local Cluster Storage location (Real Cluster Not Active Cluster)  =0 NOT available; < 0 -1 - TransportedPosition
        public int GlobalClusterNumber; // Occasionally set to Global Cluster Number
        public int PackedHost; // =0 unless Distributed when Packed Host index in form H + MULTIPLIER (H1 + MULTIPLIER * H2 ) where H must be current host labelled by ActiveCluster Index
        public double Ymapping; // Value of mapped Y last major iteration

        public ClusterIndirection(int IterationSetINPUT, int AvailabilityINPUT)
        {
            IterationSet = IterationSetINPUT;
            Availability = AvailabilityINPUT;
            GlobalClusterNumber = -1;
            PackedHost = 0;
            Ymapping = 0.0;
        }

    }   // End ClusterIndirection

   

    public class FullSolution
    {
        public int SpongeCluster;   // Sponge Cluster Location
        public double[][] CenterPosition;   // Center Position
        public int[] CreatedIndex;  // Created Index of Cluster
        public int[] CurrentNode;   // Current Node holding Cluster
        public int[] OccupationCount;   // Number of points in Cluster
        public int NumberofCenters; // Number of Centers
        public int MaximumNumberofCenters;    // maximum Number of Centers
        public int IterationSetAt;  // Iteration where this instance set at

        public static int PointsinHistogram = 50;   // Number of Points in Histograms

        public FullSolution(int MaxCenters)
        {
            SpongeCluster = -1;
            CenterPosition = new double[MaxCenters][];
            for(int CenterIndex = 0; CenterIndex < MaxCenters; CenterIndex++)
                CenterPosition[CenterIndex] = new double[Program.ParameterVectorDimension];
            CreatedIndex = new int[MaxCenters];
            CurrentNode = new int[MaxCenters];
            OccupationCount = new int[MaxCenters];
            NumberofCenters = 0;
            MaximumNumberofCenters = MaxCenters;
            IterationSetAt = -1;
        }


        public void HistogramClusterProperties()
        {   // Histogram Occupation Counts

            if (DAVectorUtility.MPI_Rank != 0)
                return;
            int maxcount = 0;
            double InterCenterDistcemax = 0.0;
            for (int FullCenterIndex1 = 0; FullCenterIndex1 < this.NumberofCenters; FullCenterIndex1++)
            {
                if( FullCenterIndex1 == SpongeCluster)
                    continue;
                maxcount = Math.Max(maxcount, this.OccupationCount[FullCenterIndex1]);
                for (int FullCenterIndex2 = FullCenterIndex1; FullCenterIndex2 < this.NumberofCenters; FullCenterIndex2++)
                {
                    if( FullCenterIndex2 == SpongeCluster)
                        continue;
                    InterCenterDistcemax = Math.Max(InterCenterDistcemax, DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterPosition[FullCenterIndex1], CenterPosition[FullCenterIndex2]));
                }
            }
            double hist1min = 0.0;
            double hist1max = maxcount;
            double hist2min = 0.0;
            double hist2max = InterCenterDistcemax;
            DAVectorUtility.SetUpHistogramRange(FullSolution.PointsinHistogram, ref hist1min, ref hist1max);
            DAVectorUtility.SetUpHistogramRange(FullSolution.PointsinHistogram, ref hist2min, ref hist2max);
            int [] Hist1 = new int[2+FullSolution.PointsinHistogram];
            int [] Hist2 = new int[2+FullSolution.PointsinHistogram];
            for (int FullCenterIndex1 = 0; FullCenterIndex1 < this.NumberofCenters; FullCenterIndex1++)
            {
                if( FullCenterIndex1 == SpongeCluster)
                    continue;
                int thiscount = FullSolution.getHistposition( (double) this.OccupationCount[FullCenterIndex1], hist1min, hist1max, FullSolution.PointsinHistogram);
                ++Hist1[thiscount];
                for (int FullCenterIndex2 = FullCenterIndex1; FullCenterIndex2 < this.NumberofCenters; FullCenterIndex2++)
                {
                    if( FullCenterIndex2 == SpongeCluster)
                        continue;
                    double tmp = DAVectorParallelism.getNOTSquaredUNScaledDistancebetweenVectors(CenterPosition[FullCenterIndex1], CenterPosition[FullCenterIndex2]);
                    int CorrelCount = FullSolution.getHistposition( tmp, hist2min, hist2max, FullSolution.PointsinHistogram);
                    ++Hist2[CorrelCount];
                }
            }
            DAVectorUtility.SALSAPrint(0, "\nOccupation Count Histogram " + HistOutput(hist1min, hist1max, FullSolution.PointsinHistogram, Hist1));
            DAVectorUtility.SALSAPrint(0, "\nCenter-Center Distance Histogram " + HistOutput(hist2min, hist2max, FullSolution.PointsinHistogram, Hist2));
        }

        public static int getHistposition(double value, double histmin, double histmax, int NumberofPoints)
        {
            if (value < histmin)
                return 0;
            if (value > histmax)
                return NumberofPoints + 1;
            int index = 1 + Convert.ToInt32( Math.Floor( (value - histmin)*NumberofPoints/(histmax-histmin) ) );
            if (index > NumberofPoints)
                index = NumberofPoints;
            return index;
        }

        public static string HistOutput(double Histmin, double Histmax, int NumberofPoints, int[] Count)
        {
            string message = "Min " + Histmin.ToString("F3") + " Max " + Histmax.ToString("F3") + " Bin Size " + ((Histmax - Histmin) / NumberofPoints).ToString("F3") + " Num " + NumberofPoints.ToString() + "\n";
            for (int histloop = 0; histloop < (NumberofPoints + 2); histloop++)
                message += Count[histloop].ToString() + " ";
            return message;
        }

    }   // End FullSolution

      public class ClusteringSolution
    {
        public bool DistributedExecutionMode = false;   // If True, run in distributed mode; if false only Global Clusters
        public double[][] M_alpha_kpointer_;                // Probability that point in Cluster pointed to by kpointer
        public int[][] Map_alpha_PointertoCreatedIndex;      // Map Pointer to Cluster Created Index for each point 
        public int[] NumClusters_alpha_;           // Number of Clusters associated with point alpha

        public double[] DiffMsummed_k_;              // Change in Malpha_k_ summed over alpha
        public double[] C_k_;                       // Summation of C(k) over parallel threads
        public double[] ClusterScaledSquaredWidth_k_;     // Mean Scaled (by Sigma_k_i_) Squared Width

        public double[] FreezingMeasure_k_;         // Freezing Measure for Cluster
        public double[] P_k_;                       // P(k) used in Continuous Clustering
        public int[] Splittable_k_;                 // If -1 used already; 0 ignored, =1 Try to find eigenvalue; =2 Eigenvalue > 0; = 3 Eigenvalue < 0
        public double[] Eigenvalue_k;               // Eigenvalue in Y_k_i_/SQRT(Sigma_k_i_) space
        public int[] OccupationCounts_k_;           // Occupation Counts for Clusters
        public int[] SplitPriority_k_;                 // -1 Infinite priority; 0 zero priority; 1 child just removed otherwise priority increases as value increases

        public double[][] Y_k_i_;                   // Y(k,i) is center of Cluster k where i runs over vector coordinates (NOT used for Sponge)
        public double[][] YPrevious_k_i_;           // YPrevious(k,i) is center of Cluster k where i runs over vector coordinates at previous iteration
        public bool[] YPreviousActuallySet;         // If true this Y has YPrevious set
        public double[][,] Correlation_k_i_j;       // Correlation that is second term in second derivation -- Symmetric -- Diagonal Sum is Square of Cluster Width
        public double[][] DC_k_DY_k_i_;             // Derivative of C[k] used in splitting
        public double[][] Sigma_k_i_;               // Squared Standard Deviation of Cluster k calculated as controlled by DAVectorSponge.SigmaMethod
        public double[][] Eigenvector_k_i;          // Splitting vector for cluster k in Y_k_/SQRT(Sigma_k_i_) -- only set for clusters that are candidates

        public int ClustertoSplit = -1;             // Current Cluster being split

        public int SpongeCluster = -1;                      // Label of Sponge
        public double PairwiseHammy = 0.0;                  // Value of Hamiltonian
        public double OldHammy = 0.0;                       //  Previous value of Hamiltonian
        public int Ncent_ThisNode = 1;                      //  the current number of clusters stored in this node
        public int Ncent_Global = 1;                        //  the current number of clusters summed over all nodes
        public double ActualCoolingFactor = Program.InitialCoolingFactor;          // Actual Cooling Factor
        public double Temperature;                          // The current temperature
        public int IterationSetAt = 0;                      // Iteration set at
        public bool Eigenvectorset = false;                 // If True Eigenvector set
        public int YPreviousSet = -1;                   // If -1 Y not set, = 0 Y but not YPrevious set, = 1 Y and YPrevious set
        public bool SolutionSet = false;                    // If True Solution Set
        public bool CorrelationsSet = false;                    // If true Correlation Set
        public int DiffMalpha_k_Set = -1;                   // -1 initial, 0 set Previous_Malpha_k_, 1 Set 
        public double[] AverageWidth;                       // Average widths over Clusters in each dimension
        public double TotaloverVectorIndicesAverageWidth;                    // Average Width summed over all dimensions

        public int[] LocalCreatedIndex;   // Created Index of Local Cluster
        public int[] LocalStatus;   // Status -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here 3 Distributed Stored Elsewhere (Not Possible)
        public int[] LocalSplitCreatedIndex;    // 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED

        //  Note quantities below are static and NOT saved in backup. However they change during job
        //  Backup must be followed by a master synchronization setting these variables
        public static int NumberLocalActiveClusters;     // Number of Active Clusters stored here -- arrays later on are defined for this range NOT NumberAvailableActiveClusters
        public static int NumberAvailableActiveClusters;     // Number of Active Clusters stored here or transported
        public static int NumberGlobalClusters; // Number of status 0 and 1 Clusters
        public static int[] RealClusterIndices;         // Maps Active Cluster Indices into Real cluster list 
        public static int[] ActiveClusterIndices;   // Cluster Numbers of Active Clusters -- maps  index of main arrays into ActiveCluster index
        public static int[] LocalNodeAccPosition;    // Current Node Accumulation Positions for this cluster labelled by ActiveCluster Index
        public static int[][] LocalThreadAccPosition;    // Current Thread Accumulation Position for this cluster labelled by ActiveCluster Index
        public static int[] LocalHost;    // =0 unless Distributed when Packed Host index in form H + MULTIPLIER (H1 + MULTIPLIER * H2 ) where H must be current host labelled by ActiveCluster Index

        //  These quantities are really static. Set at Start and do not change
        public static int NumberofPointsinProcess = -1;  // Number of Points in Process
        public static int MaximumCreatedClusters;   // Number of allowed created clusters
        public static int MaximumNumberClusterspernode = 0;     // Maximum Number of Centers per node summed over points
        public static int MaximumNumberClustersTotal = 0;     // Maximum Number of Centers per node summed over points
        public static int TargetClustersperPoint = 0;   // Target maximum number of clusters per point
        public static int TargetMinimumClustersperPoint = 0;   // Target minimum number of clusters per point
        public static int MaximumClustersperPoint = 0;   // Actual maximum number of clusters per point
        public static double ExponentArgumentCut = 20.0;     // Include all clusters with (Y(cluster)-X(point))^2 / (2 * Temperature) < ExpArgumentCut in Point Collection
        public static int cachelinesize = 0;             // Cacheline increment
       
        public static int CurrentIteration; // Current Iteration Number to monitor which clusters stale
        public static ClusterIndirection[] UniversalMapping;    // Maps CreatedIndex into location -- in distributed model any one node only has partial information 
        public static int CenterMaxforCreatedIndex = 0; // Center count to use for CreatedIndex
        public static FullSolution TotalClusterSummary; // Summary of Clusters

        public static int PACKINGMASK;          // 2^PACKINGSHIFT - 1
        public static int PACKINGSHIFT;
        public static int PACKINGMULTIPLIER;    // 2^PACKINGSHIFT

        public static int ClustersDeleted = 0;  // Clusters Deleted this iteration
        public static int ClustersMoved = 0;    // Clusters Moved this iteration
        public static int ClustersSplit = 0;    // Clusters Split this iteration

        // public static int DEBUGPrint1 = 0;

        public static void SetParameters(int NumberofPointsinProcessINPUT, int MaximumCreatedClustersINPUT, int MaximumNumberClustersTotalINPUT,
            int MaximumNumberClusterspernodeINPUT,  int cachelinesizeINPUT,
            int TargetClustersperPointINPUT, int MinimumNumberClusterspernodeINPUT, int MaximumClustersperPointINPUT, double ExponentArgumentCutINPUT)
        {
            if (NumberofPointsinProcess > 0)
                return;
            NumberofPointsinProcess = NumberofPointsinProcessINPUT;
            MaximumCreatedClusters = MaximumCreatedClustersINPUT;
            MaximumNumberClusterspernode = MaximumNumberClusterspernodeINPUT;
            TargetMinimumClustersperPoint = MinimumNumberClusterspernodeINPUT;
            MaximumNumberClustersTotal = MaximumNumberClustersTotalINPUT;
            cachelinesize = cachelinesizeINPUT;
            TargetClustersperPoint = TargetClustersperPointINPUT;
            MaximumClustersperPoint = MaximumClustersperPointINPUT;
            ExponentArgumentCut = ExponentArgumentCutINPUT;

            LocalNodeAccPosition = new int[MaximumNumberClusterspernode];
            LocalThreadAccPosition = new int[MaximumNumberClusterspernode][];
            LocalHost = new int[MaximumNumberClusterspernode];
            for (int ClusterIndex = 0; ClusterIndex < MaximumNumberClusterspernode; ClusterIndex++) 
                LocalThreadAccPosition[ClusterIndex] = new int[DAVectorUtility.ThreadCount];

            CurrentIteration = 0;

            PACKINGSHIFT = 0;
            PACKINGMASK = 1;
            while (PACKINGMASK <= DAVectorUtility.MPI_Size)
            {
                PACKINGMASK = 2 * PACKINGMASK;
                PACKINGSHIFT++;
            }
            PACKINGMULTIPLIER = PACKINGMASK;
            PACKINGMASK--;


            // Set Arrays that Map Created Indices to Local Storage
            int UniversalSize = (1 + MaximumCreatedClusters) * PACKINGMULTIPLIER;
            DAVectorUtility.SALSAPrint(0, "Create Universal Mapping " + UniversalSize.ToString());
            UniversalMapping = new ClusterIndirection[UniversalSize];
            DAVectorUtility.SALSAPrint(0, "Set up Full Solution " + MaximumNumberClustersTotal);
            TotalClusterSummary = new FullSolution(MaximumNumberClustersTotal);
            DAVectorUtility.SALSAPrint(0, "End ClusteringSolution");
        }

        public static int SetCreatedIndex(int RealLocalClusterIndex)
        {
            int host = 0;
            if (ParallelClustering.RunningSolution.DistributedExecutionMode && ((ParallelClustering.RunningSolution.SpongeCluster < 0) || (ParallelClustering.RunningSolution.SpongeCluster != RealLocalClusterIndex)) )
                host = DAVectorUtility.MPI_Rank;
            ++CenterMaxforCreatedIndex;
            if(CenterMaxforCreatedIndex > ClusteringSolution.MaximumCreatedClusters)
            {
                Exception e = DAVectorUtility.SALSAError(" Center Index too large " + ClusteringSolution.CenterMaxforCreatedIndex.ToString() + " "
                    + ParallelClustering.RunningSolution.Ncent_ThisNode.ToString() + " " + ClusteringSolution.MaximumNumberClustersTotal.ToString() + " " + ClusteringSolution.MaximumCreatedClusters.ToString() );
                throw (e);
            }
            int CreatedIndex = host + 1 + (CenterMaxforCreatedIndex << ClusteringSolution.PACKINGSHIFT);
            UniversalMapping[CreatedIndex] = new ClusterIndirection(CurrentIteration, RealLocalClusterIndex + 1);
            ParallelClustering.RunningSolution.LocalCreatedIndex[RealLocalClusterIndex] = CreatedIndex;
            return CreatedIndex;

        }   // End SetCreatedIndex()

        public static int MapActivetoCreatedIndex(int ActiveClusterIndex, ClusteringSolution Solution)
        {
            if (Solution.DistributedExecutionMode && (ActiveClusterIndex >= NumberLocalActiveClusters))
                return DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedCreatedIndex[ActiveClusterIndex - NumberLocalActiveClusters];
            else
                return Solution.LocalCreatedIndex[RealClusterIndices[ActiveClusterIndex]];
                

        }   // End MapActivetoCreatedIndex

        public ClusteringSolution(bool thereisasponge)
        {
            if (NumberofPointsinProcess < 0)
            {
                Exception e = DAVectorUtility.SALSAError("NumberofPointsinProcess Unset");
                throw (e);
            }
            this.Ncent_ThisNode = 1;
            this.SpongeCluster = -1;
            int FirstRealCluster = 0;
            this.DistributedExecutionMode = false;

            if (thereisasponge)
            {
                this.SpongeCluster = 0;
                this.Ncent_ThisNode = 2;
                FirstRealCluster = 1;
            }
            this.Ncent_Global = this.Ncent_ThisNode;

            //  Point Arrays
            DAVectorUtility.SALSAPrint(0, "Large Arrays Started " + NumberofPointsinProcess.ToString() + " Clusters per Point " + MaximumClustersperPoint.ToString());
            M_alpha_kpointer_ = new double[NumberofPointsinProcess][];
            Map_alpha_PointertoCreatedIndex = new int[NumberofPointsinProcess][];
            NumClusters_alpha_ = new int[NumberofPointsinProcess];

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    M_alpha_kpointer_[alpha] = new double[MaximumClustersperPoint];
                    Map_alpha_PointertoCreatedIndex[alpha] = new int[MaximumClustersperPoint];
                    Map_alpha_PointertoCreatedIndex[alpha][0] = 0;
                    Map_alpha_PointertoCreatedIndex[alpha][FirstRealCluster] = FirstRealCluster;
                    NumClusters_alpha_[alpha] = 1+FirstRealCluster;
                }


            }); // End loop initialing Point dependent quantities
            DAVectorUtility.SALSAPrint(0, "Large Arrays Created " + NumberofPointsinProcess.ToString());

            DiffMsummed_k_ = new double[MaximumNumberClusterspernode + cachelinesize]; 
            C_k_ = new double[MaximumNumberClusterspernode + cachelinesize]; 
            ClusterScaledSquaredWidth_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
            P_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
            FreezingMeasure_k_ = new double[MaximumNumberClusterspernode + cachelinesize];
            OccupationCounts_k_ = new int[MaximumNumberClusterspernode + cachelinesize]; 
            Splittable_k_ = new int[MaximumNumberClusterspernode + cachelinesize];
            SplitPriority_k_ = new int[MaximumNumberClusterspernode + cachelinesize];
            Eigenvalue_k = new double[MaximumNumberClusterspernode + cachelinesize]; 
            ClustertoSplit = -1;
            
            Y_k_i_ = new double[MaximumNumberClusterspernode][];                   // Y(k,i) is center of Cluster k where i runs over vector coordinates
            YPrevious_k_i_ = new double[MaximumNumberClusterspernode][];
            YPreviousActuallySet = new bool[MaximumNumberClusterspernode];
            Eigenvector_k_i = new double[MaximumNumberClusterspernode][];
            DC_k_DY_k_i_ = new double[MaximumNumberClusterspernode][];
            Sigma_k_i_ = new double[MaximumNumberClusterspernode][];                   // Standard Deviation Squared of Cluster k calculated as controlled by DAVectorSponge.SigmaMethod
            Correlation_k_i_j = new double[MaximumNumberClusterspernode][,];       // Correlation that is second term in second derivation -- Symmetric -- Diagonal Sum is Square of Cluster Width

            RealClusterIndices = new int[MaximumNumberClusterspernode];
            ActiveClusterIndices = new int[MaximumNumberClusterspernode];
            LocalCreatedIndex = new int[MaximumNumberClusterspernode];
            LocalStatus = new int[MaximumNumberClusterspernode];
            LocalSplitCreatedIndex = new int[MaximumNumberClusterspernode];

            for (int ClusterIndex = 0; ClusterIndex < MaximumNumberClusterspernode; ClusterIndex++)
            {
                YPreviousActuallySet[ClusterIndex] = false;
                Y_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
                YPrevious_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
                Eigenvector_k_i[ClusterIndex] = new double[Program.ParameterVectorDimension];
                DC_k_DY_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
                Sigma_k_i_[ClusterIndex] = new double[Program.ParameterVectorDimension];
                int correlsize = 1;
                if (Program.CalculateCorrelationMatrix)
                    correlsize = Program.ParameterVectorDimension;
                Correlation_k_i_j[ClusterIndex] = new double[correlsize, correlsize];
            }

            AverageWidth = new double[Program.ParameterVectorDimension];
            Eigenvectorset = false;
            SolutionSet = false;
            CorrelationsSet = false;
            DiffMalpha_k_Set = -1;


            if (thereisasponge)
            {
                SplitPriority_k_[this.SpongeCluster] = 0;
                Splittable_k_[this.SpongeCluster] = 0;
                LocalStatus[this.SpongeCluster] = 0;
            }
            SetActiveClusters();
            DAVectorUtility.SALSAPrint(0, "Clusters set up");

        }   // End ClusteringSolution
        
        // Set Active Cluster Array and count
        public void SetActiveClusters()
        {   // Update Ncent_Global and basic static arrays for locally controlled clusters
            //  Need to set accumulation positions and Host index separately
            //  Need to set remote clusters transported locally separately in 
            int ActiveCount = 0;
            int ActiveCount2 = 0;
            for (int RealClusterIndex = 0; RealClusterIndex < this.Ncent_ThisNode; RealClusterIndex++)
            {
                ActiveClusterIndices[RealClusterIndex] = -1;
                if (this.LocalStatus[RealClusterIndex] < 0 || this.LocalStatus[RealClusterIndex] > 2)
                    continue;
                if( this.LocalStatus[RealClusterIndex] == 2 )
                    ++ActiveCount2;
                ActiveClusterIndices[RealClusterIndex] = ActiveCount;
                RealClusterIndices[ActiveCount] = RealClusterIndex;
                ++ActiveCount;
            }
            NumberLocalActiveClusters = ActiveCount;
            NumberAvailableActiveClusters = ActiveCount;    // This is incremented on Major Synchronizations
            int ActiveCount01 = ActiveCount - ActiveCount2;
            if (this.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
                ActiveCount2 = DAVectorUtility.MPI_communicator.Allreduce<int>(ActiveCount2, Operation<int>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
            }

            //   Check agreement on Global Count
            NumberGlobalClusters = ActiveCount01;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
            NumberGlobalClusters = DAVectorUtility.MPI_communicator.Allreduce<int>(NumberGlobalClusters, Operation<int>.Max);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
            if(ActiveCount01 != NumberGlobalClusters)
            {
                Exception e = DAVectorUtility.SALSAError("Error in Global Count Max " + NumberGlobalClusters.ToString() + " Actual " + ActiveCount01.ToString());
                throw e;
            }
            this.Ncent_Global = ActiveCount2 + ActiveCount01;
            return;

        }   // End SetActiveClusters()

        //  These widths are divided by Sigma_k_i_ and NOT Square Rooted
        public void SetClusterWidths()
        {
            DAVectorUtility.StartSubTimer(4);
            GlobalReductions.FindIndirectVectorDoubleSum FindWidthsGlobal = null;
            DistributedReductions.FindIndirectMultiVectorDoubleSum FindWidthsinDistributedMode = null;
            int BeginFindWidthsinDistributedMode = -1;

            if (this.DistributedExecutionMode)
            {
                FindWidthsinDistributedMode = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginFindWidthsinDistributedMode = FindWidthsinDistributedMode.AddComponents(1);
                FindWidthsinDistributedMode.NodeInitialize();
            }
            else
                FindWidthsGlobal = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);

            GlobalReductions.FindVectorDoubleSum FindAverageWidth = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, Program.ParameterVectorDimension);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                if (this.DistributedExecutionMode)
                    FindWidthsinDistributedMode.ThreadInitialize((int)ThreadNo);
                else
                    FindWidthsGlobal.startthread((int)ThreadNo);
                int ThreadStorePosition = -1;

                double[] Sigma_Pointer;
                double[] Y_Pointer;

                double[] WorkSpace_Avg = new double[Program.ParameterVectorDimension];
                int ArraySize = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] WorkSpace_k_ = new double[ArraySize];

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                        if (RemoteIndex < 0)
                        {
                            Y_Pointer = this.Y_k_i_[RealClusterIndex];
                            Sigma_Pointer = this.Sigma_k_i_[RealClusterIndex];
                            if (this.DistributedExecutionMode)
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
                        }
                        else
                        {   // Only possible for distributed execution
                            Y_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex];
                            Sigma_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex];
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo];
                        }
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;

                        double wgt = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        if ((SpongeCluster < 0) || (RealClusterIndex != SpongeCluster))
                        {
                            WorkSpace_k_[IndirectClusterIndex] = wgt * DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(alpha, ActiveClusterIndex, this);
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            {
                                double tmp = Program.PointPosition[alpha][VectorIndex] - Y_Pointer[VectorIndex];
                                WorkSpace_Avg[VectorIndex] = wgt * tmp * tmp / Sigma_Pointer[VectorIndex];
                            }
                            FindAverageWidth.addapoint(ThreadNo, WorkSpace_Avg);
                        }
                        else
                            WorkSpace_k_[IndirectClusterIndex] = wgt * Program.SpongeFactor * Program.SpongeFactor;
                        if (this.DistributedExecutionMode)
                            FindWidthsinDistributedMode.addapoint((int)ThreadNo, ThreadStorePosition, BeginFindWidthsinDistributedMode, WorkSpace_k_[IndirectClusterIndex]);
                    }
                    if (!this.DistributedExecutionMode)
                        FindWidthsGlobal.addapoint(ThreadNo, IndirectSize, ActiveClustersperPoint, WorkSpace_k_);

                }
            }); // End loop initialing Point dependent quantities

            if (this.DistributedExecutionMode)
                FindWidthsinDistributedMode.sumoverthreadsandmpi();
            else
                FindWidthsGlobal.sumoverthreadsandmpi();

            // Broadcast zero size information for global clusters
            bool[] ZeroSizeClusters = new bool[ClusteringSolution.NumberGlobalClusters];
            int CountGlobalClusters = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if ((ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 1))
                    continue;
                bool zerosizecluster = false;
                if (RealClusterIndex != ParallelClustering.RunningSolution.SpongeCluster)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
                ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
                ++CountGlobalClusters;
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            DAVectorUtility.MPI_communicator.Broadcast<bool>(ref ZeroSizeClusters, 0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            CountGlobalClusters = 0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                double Result;
                if (this.DistributedExecutionMode)
                    Result = FindWidthsinDistributedMode.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindWidthsinDistributedMode];
                else
                    Result = FindWidthsGlobal.TotalVectorSum[LocalActiveClusterIndex];

                bool zerosizecluster = false;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                {
                    zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                    ++CountGlobalClusters;
                }
                else if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;

                if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                {
                    this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = Program.SpongeFactor * Program.SpongeFactor;
                    continue;
                }

                if (zerosizecluster)
                    this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = 0.0;
                else
                    this.ClusterScaledSquaredWidth_k_[RealClusterIndex] = Result / this.C_k_[RealClusterIndex];
            }

            double C_k_Sum1not0 = 0.0;
            double C_k_Sum2 = 0.0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (this.LocalStatus[RealClusterIndex] == 2)
                    C_k_Sum2 += this.C_k_[RealClusterIndex];
                else if (this.LocalStatus[RealClusterIndex] == 1)
                    C_k_Sum1not0 += this.C_k_[RealClusterIndex];
            }
            if (this.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
                C_k_Sum2 = DAVectorUtility.MPI_communicator.Allreduce<double>(C_k_Sum2, Operation<double>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);
            }
            double C_Sum = C_k_Sum1not0 + C_k_Sum2;

            //   No contribution from Sponge here
            FindAverageWidth.sumoverthreadsandmpi();
            this.TotaloverVectorIndicesAverageWidth = 0.0;
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                this.AverageWidth[VectorIndex] = FindAverageWidth.TotalVectorSum[VectorIndex] / C_Sum;
                this.TotaloverVectorIndicesAverageWidth += this.AverageWidth[VectorIndex];
            }
            DAVectorUtility.StopSubTimer(4);

        }   // End SetClusterWidths

        public double SetClusterRadius(int RealClusterIndex)
        {   // Note this definition -- Sqrt(squared width) is different from Maximum distance of Points in Cluster from Center as for DA we don't assign rigourously points to clusters

            return Math.Sqrt(this.ClusterScaledSquaredWidth_k_[RealClusterIndex]);
        }

        //  Set the C_k_ for Clusters in this Solution
        public void SetClusterSizes()
        {
            DAVectorUtility.StartSubTimer(4);
            GlobalReductions.FindIndirectVectorDoubleSum FindSizes = null;
            GlobalReductions.FindIndirectVectorDoubleSum FindFreezing = null;
            DistributedReductions.FindIndirectMultiVectorDoubleSum Find2Components = null;
            int BeginFindSizes = -1;
            int BeginFindFreezing = -1;

            if (this.DistributedExecutionMode)
            {
                Find2Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginFindSizes = Find2Components.AddComponents(1);
                BeginFindFreezing = Find2Components.AddComponents(1);
                Find2Components.NodeInitialize();
            }
            else
            {
                FindSizes = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
                FindFreezing = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                if (this.DistributedExecutionMode)
                    Find2Components.ThreadInitialize((int)ThreadNo);
                else
                {
                    FindSizes.startthread((int)ThreadNo);
                    FindFreezing.startthread((int)ThreadNo);
                }
                int ThreadStorePosition = -1;

                int ArraySize = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[] WorkSpace_k_ = new double[ArraySize];

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                        if (RemoteIndex < 0)
                        {
                            if (this.DistributedExecutionMode)
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
                        }
                        else
                        {   // Only possible for distributed execution
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo];
                        }
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        WorkSpace_k_[IndirectClusterIndex] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] * (1.0 - this.M_alpha_kpointer_[alpha][IndirectClusterIndex]);

                        if (this.DistributedExecutionMode)
                        {
                            if ((ThreadStorePosition < 0) || (ThreadStorePosition >= DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[ThreadNo]))
                            {
                                Exception e = DAVectorUtility.SALSAError(" Find Sizes Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString()
                                    + " Illegal Position " + ThreadStorePosition.ToString() + " Max " + DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[ThreadNo].ToString()
                                    + " CreatedIndex " + ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex].ToString()
                                    + " Node Total " + DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode.ToString() );
                                throw (e);
                            }
                            Find2Components.addapoint((int)ThreadNo, ThreadStorePosition, BeginFindSizes, this.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                            Find2Components.addapoint((int)ThreadNo, ThreadStorePosition, BeginFindFreezing, WorkSpace_k_[IndirectClusterIndex]);
                        }
                    }
                    if (!this.DistributedExecutionMode)
                    {
                        FindSizes.addapoint(ThreadNo, IndirectSize, ActiveClustersperPoint, this.M_alpha_kpointer_[alpha]);
                        FindFreezing.addapoint(ThreadNo, IndirectSize, ActiveClustersperPoint, WorkSpace_k_);
                    }
                }
            }); // End loop initialing Point dependent quantities

            if (this.DistributedExecutionMode)
                Find2Components.sumoverthreadsandmpi();
            else
            {
                FindSizes.sumoverthreadsandmpi();
                FindFreezing.sumoverthreadsandmpi();
            }

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (this.DistributedExecutionMode)
                    this.C_k_[RealClusterIndex] = Find2Components.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindSizes];
                else
                    this.C_k_[RealClusterIndex] = FindSizes.TotalVectorSum[LocalActiveClusterIndex];
            }

            // Broadcast zero size information for global clusters
            bool[] ZeroSizeClusters = new bool[ClusteringSolution.NumberGlobalClusters];
            int CountGlobalClusters = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if ((ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 1))
                    continue;
                bool zerosizecluster = false;
                if (RealClusterIndex != ParallelClustering.RunningSolution.SpongeCluster)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
                ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
                ++CountGlobalClusters;
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            DAVectorUtility.MPI_communicator.Broadcast<bool>(ref ZeroSizeClusters, 0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            CountGlobalClusters = 0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (this.DistributedExecutionMode)
                    this.FreezingMeasure_k_[RealClusterIndex] = Find2Components.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][BeginFindFreezing];
                else
                    this.FreezingMeasure_k_[RealClusterIndex] = FindFreezing.TotalVectorSum[LocalActiveClusterIndex];

                bool zerosizecluster = false;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                {
                    zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                    ++CountGlobalClusters;
                }
                else if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;

                if (!zerosizecluster)
                    this.FreezingMeasure_k_[RealClusterIndex] = this.FreezingMeasure_k_[RealClusterIndex] / this.C_k_[RealClusterIndex];
            }

            DAVectorUtility.StopSubTimer(4);

        }   // End SetClusterSizes


        //  Set the Correlation_k_i_j for Clusters in this Solution
        // Note NOT Divided by C_k_ for a true correlation as division not wanted in eigencalculation
        public void SetClusterCorrelations()
        {
            DAVectorUtility.StartSubTimer(4);
            int NumberofCorrelationPositions = 0;
            if(Program.CalculateCorrelationMatrix)
                NumberofCorrelationPositions = (Program.ParameterVectorDimension * (Program.ParameterVectorDimension + 1)) / 2;

            GlobalReductions.FindIndirectVectorDoubleSum[] FindDC_k_DY_k_ = null;
            GlobalReductions.FindIndirectVectorDoubleSum[] FindCorrelations = null;
            DistributedReductions.FindIndirectMultiVectorDoubleSum Find2Components = null;
            int BeginDC_k_DY_k_ = -1;
            int BeginCorrelations = -1;

            if (this.DistributedExecutionMode)
            {
                Find2Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginDC_k_DY_k_ = Find2Components.AddComponents(Program.ParameterVectorDimension);
                if (Program.CalculateCorrelationMatrix)
                    BeginCorrelations = Find2Components.AddComponents(NumberofCorrelationPositions);
                Find2Components.NodeInitialize();
            }
            else
            {
                FindDC_k_DY_k_ = new GlobalReductions.FindIndirectVectorDoubleSum[Program.ParameterVectorDimension];
                if (Program.CalculateCorrelationMatrix)
                    FindCorrelations = new GlobalReductions.FindIndirectVectorDoubleSum[NumberofCorrelationPositions];

                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    FindDC_k_DY_k_[VectorIndex] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
                if (Program.CalculateCorrelationMatrix)
                {
                    for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                        FindCorrelations[CorrelationPositions] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);
                }
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                if (this.DistributedExecutionMode)
                    Find2Components.ThreadInitialize((int)ThreadNo);
                else
                {
                    if (Program.CalculateCorrelationMatrix)
                    {
                        for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                            FindCorrelations[CorrelationPositions].startthread((int)ThreadNo);
                    }
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        FindDC_k_DY_k_[VectorIndex].startthread((int)ThreadNo);
                }

                double[] Sigma_Pointer;
                double[] Y_Pointer;

                int ThreadStorePosition = -1;
                int ArraySize = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];

                double[] DC_k_DY_k_Values = new double[Program.ParameterVectorDimension];
                double[][] DC_k_DY_k_Values_ClusterIndex = new double[Program.ParameterVectorDimension][];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    DC_k_DY_k_Values_ClusterIndex[VectorIndex] = new double[ArraySize];

                int Maynotneed = Math.Max(1, NumberofCorrelationPositions);
                double[][] WorkSpace_Sym_i_j_ClusterIndex = new double[Maynotneed][];
                double[] WorkSpace_Sym_i_j_ = new double[Maynotneed];
                if (Program.CalculateCorrelationMatrix)
                {
                    for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                        WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions] = new double[ArraySize];
                }

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;

                        if (RemoteIndex < 0)
                        {
                            Y_Pointer = this.Y_k_i_[RealClusterIndex];
                            Sigma_Pointer = this.Sigma_k_i_[RealClusterIndex];
                            if (this.DistributedExecutionMode)
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
                        }
                        else
                        {   // Only possible for distributed execution
                            Y_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex];
                            Sigma_Pointer = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex];
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo];
                        }

                        if ((RealClusterIndex == SpongeCluster) && (SpongeCluster >= 0))
                        {
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                                DC_k_DY_k_Values_ClusterIndex[VectorIndex][IndirectClusterIndex] = 0.0;
                            if (Program.CalculateCorrelationMatrix)
                            {
                                for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                                    WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions][IndirectClusterIndex] = 0.0;
                            }
                            continue;
                        }
                        //	0.5 as cluster halved
                        double CurrentMvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] * 0.5;
                        int countdoubleindex = 0;

                        for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                        {
                            double fudge = 1.0 / (Math.Sqrt(Sigma_Pointer[VectorIndex1]) * this.Temperature);
                            DC_k_DY_k_Values[VectorIndex1] = -(Y_Pointer[VectorIndex1] - Program.PointPosition[alpha][VectorIndex1]) * fudge * CurrentMvalue * (1.0 - CurrentMvalue);
                            DC_k_DY_k_Values_ClusterIndex[VectorIndex1][IndirectClusterIndex] = DC_k_DY_k_Values[VectorIndex1];

                            if (Program.CalculateCorrelationMatrix)
                            {
                                double vecdiff1 = Program.PointPosition[alpha][VectorIndex1] - Y_Pointer[VectorIndex1];
                                for (int VectorIndex2 = VectorIndex1; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                                {
                                    double vecdiff2 = Program.PointPosition[alpha][VectorIndex2] - Y_Pointer[VectorIndex2];
                                    WorkSpace_Sym_i_j_[countdoubleindex] = CurrentMvalue * vecdiff1 * vecdiff2;
                                    WorkSpace_Sym_i_j_ClusterIndex[countdoubleindex][IndirectClusterIndex] = WorkSpace_Sym_i_j_[countdoubleindex];
                                    ++countdoubleindex;
                                }
                            }
                        }
                        if (this.DistributedExecutionMode)
                        {
                            if ((ThreadStorePosition < 0) || (ThreadStorePosition >= DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[ThreadNo]))
                            {
                                Exception e = DAVectorUtility.SALSAError(" Find Correlations Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString()
                                    + " Illegal Position " + ThreadStorePosition.ToString() + " Max " + DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[ThreadNo].ToString()
                                    + " CreatedIndex " + ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex].ToString()
                                    + " Node Total " + DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperNode.ToString());
                                throw (e);
                            }
                            if (Program.CalculateCorrelationMatrix)
                                Find2Components.addapoint((int)ThreadNo, ThreadStorePosition, BeginCorrelations, NumberofCorrelationPositions, WorkSpace_Sym_i_j_);
                            Find2Components.addapoint((int)ThreadNo, ThreadStorePosition, BeginDC_k_DY_k_, Program.ParameterVectorDimension, DC_k_DY_k_Values);
                        }
                    }
                    if (!this.DistributedExecutionMode)
                    {
                        if (Program.CalculateCorrelationMatrix)
                        {
                            for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                                FindCorrelations[CorrelationPositions].addapoint((int)ThreadNo, IndirectSize, ActiveClustersperPoint, WorkSpace_Sym_i_j_ClusterIndex[CorrelationPositions]);
                        }
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            FindDC_k_DY_k_[VectorIndex].addapoint((int)ThreadNo, IndirectSize, ActiveClustersperPoint, DC_k_DY_k_Values_ClusterIndex[VectorIndex]);
                    }
                }
            }); // End loop initialing Point dependent quantities

            if (this.DistributedExecutionMode)
                Find2Components.sumoverthreadsandmpi();
            else
            {
                if (Program.CalculateCorrelationMatrix)
                {
                    for (int CorrelationPositions = 0; CorrelationPositions < NumberofCorrelationPositions; CorrelationPositions++)
                        FindCorrelations[CorrelationPositions].sumoverthreadsandmpi();
                }
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    FindDC_k_DY_k_[VectorIndex].sumoverthreadsandmpi();
            }


            int AccumulationPosition = -1;
            double result = 0.0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (RealClusterIndex == SpongeCluster)
                    continue;
                if (this.DistributedExecutionMode)
                    AccumulationPosition = ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex];

                int countdoubleindex = 0;
                for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                {
                    if (this.DistributedExecutionMode)
                        this.DC_k_DY_k_i_[RealClusterIndex][VectorIndex1] = Find2Components.TotalVectorSum[AccumulationPosition][BeginDC_k_DY_k_ + VectorIndex1];
                    else
                        this.DC_k_DY_k_i_[RealClusterIndex][VectorIndex1] = FindDC_k_DY_k_[VectorIndex1].TotalVectorSum[LocalActiveClusterIndex];
                    if (Program.CalculateCorrelationMatrix)
                    {
                        for (int VectorIndex2 = VectorIndex1; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                        {
                            if (this.DistributedExecutionMode)
                                result = Find2Components.TotalVectorSum[AccumulationPosition][BeginCorrelations + countdoubleindex];
                            else
                                result = FindCorrelations[countdoubleindex].TotalVectorSum[LocalActiveClusterIndex];
                            this.Correlation_k_i_j[RealClusterIndex][VectorIndex1, VectorIndex2] = result / (this.Temperature * Math.Sqrt(this.Sigma_k_i_[RealClusterIndex][VectorIndex1] * this.Sigma_k_i_[RealClusterIndex][VectorIndex2]));
                            if (VectorIndex1 != VectorIndex2)
                                this.Correlation_k_i_j[LocalActiveClusterIndex][VectorIndex2, VectorIndex1]
                                = this.Correlation_k_i_j[LocalActiveClusterIndex][VectorIndex1, VectorIndex2];
                            ++countdoubleindex;
                        }
                    }
                }
            }
            this.CorrelationsSet = true;
            DAVectorUtility.StopSubTimer(4);

        }   // End SetClusterCorrelations()

        public void FindOccupationCounts()
        {
            DAVectorUtility.StartSubTimer(4); 
            GlobalReductions.FindDoubleArraySum GetOccupationCounts = null;
            DistributedReductions.FindIndirectMultiVectorDoubleSum FindOccupationCounts = null;
            int BeginOccupationCounts = -1;

            if (this.DistributedExecutionMode)
            {
                FindOccupationCounts = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginOccupationCounts = FindOccupationCounts.AddComponents(1);
                FindOccupationCounts.NodeInitialize();
            }
            else
                GetOccupationCounts = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, ClusteringSolution.NumberLocalActiveClusters);

            //  Parallel Section setting cluster occupation counts
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                if (this.DistributedExecutionMode)
                    FindOccupationCounts.ThreadInitialize((int)ThreadNo);
                else
                    GetOccupationCounts.startthread((int)ThreadNo);
                int ThreadStorePosition = -1;

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    double MaxMvalue = -1.0;
                    int NearestClusterIndirectIndex = -1;
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        LegalCluster(alpha, IndirectClusterIndex);
                        if (this.M_alpha_kpointer_[alpha][IndirectClusterIndex] > MaxMvalue)
                        {
                            MaxMvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            NearestClusterIndirectIndex = IndirectClusterIndex;
                        }
                    }

                    if( NearestClusterIndirectIndex < 0)
                    {
                        string msg = "";
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                            msg +=  " " + this.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4");
                        Exception e = DAVectorUtility.SALSAError("Error due to no nearby Cluster for Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " " + IndirectSize.ToString() 
                            + " Process " + DAVectorUtility.MPI_Rank.ToString() + msg);
                        throw e;
                    }
                    int RealClusterIndex = -1;
                    int RemoteIndex = -1;
                    int ActiveClusterIndex = -1;
                    VectorAnnealIterate.ClusterPointersforaPoint(alpha, NearestClusterIndirectIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                    if (RemoteIndex < 0)
                    {
                        if (this.DistributedExecutionMode)
                            ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
                    }
                    else
                    {   // Only possible for distributed execution
                        ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo];
                    }
                    if (this.DistributedExecutionMode)
                        FindOccupationCounts.addapoint(ThreadNo, ThreadStorePosition, 1.0);
                    else
                        GetOccupationCounts.addapoint((int)ThreadNo, ActiveClusterIndex);
                }
            });  // End Parallel Section setting cluster labels

            double Result = 0.0;
            if (this.DistributedExecutionMode)
                FindOccupationCounts.sumoverthreadsandmpi();
            else
                GetOccupationCounts.sumoverthreadsandmpi();

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                if (this.DistributedExecutionMode)
                    Result = FindOccupationCounts.TotalVectorSum[ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex]][0];
                else
                    Result = GetOccupationCounts.TotalSum[LocalActiveClusterIndex];
                this.OccupationCounts_k_[ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex]] = (int)Math.Floor(Result + 0.5);
            }
            DAVectorUtility.StopSubTimer(4);

        }   // End FindOccupationCounts()

        public static void SetGlobalClusterNumbers()
        {   // Find Cluster Numbers starting at 0 in Host 0 and incrementing through hosts systematically
            // Set TotalClusterSummary

            if (TotalClusterSummary.IterationSetAt == ParallelClustering.RunningSolution.IterationSetAt)
                return;
            TotalClusterSummary.IterationSetAt = ParallelClustering.RunningSolution.IterationSetAt;
            ParallelClustering.RunningSolution.FindOccupationCounts();

            int[] LocalClusterRealIndixes = new int[ClusteringSolution.NumberLocalActiveClusters];
            int CountLocal = 0;
            int CountGlobal = 0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 0)
                    TotalClusterSummary.SpongeCluster = CountGlobal;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                {
                    LocalClusterRealIndixes[CountLocal] = RealClusterIndex;
                    ++CountLocal;
                }
                else
                {
                    UniversalMapping[ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]].GlobalClusterNumber = CountGlobal;
                    TotalClusterSummary.CreatedIndex[CountGlobal] = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        TotalClusterSummary.CenterPosition[CountGlobal][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                    TotalClusterSummary.CurrentNode[CountGlobal] = 0;
                    TotalClusterSummary.OccupationCount[CountGlobal] = ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex];
                    ++CountGlobal;
                }
            }
            int BaseCountGlobal = CountGlobal;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            DAVectorUtility.MPI_communicator.Broadcast<int>(ref CountGlobal, 0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                TotalClusterSummary.NumberofCenters = CountGlobal;
                return;
            }

            int Maxsize = CountLocal;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming2);
            Maxsize = DAVectorUtility.MPI_communicator.Allreduce<int>(Maxsize, Operation<int>.Max);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming2);

            MPIPacket<int> FromAfarInt = new MPIPacket<int>(Maxsize);
            MPIPacket<double> FromAfarDouble = new MPIPacket<double>(Maxsize);

            for (int SourceNode = 0; SourceNode < DAVectorUtility.MPI_Size; SourceNode++)
            {
                int SaveCountGlobal = CountGlobal;

                // Transport Created Index
                if (SourceNode == DAVectorUtility.MPI_Rank)
                {
                    FromAfarInt.NumberofPoints = CountLocal;
                    FromAfarInt.FirstPoint = 0;
                    for (int FromAfarIndex = 0; FromAfarIndex < FromAfarInt.NumberofPoints; FromAfarIndex++)
                        FromAfarInt.Marray[FromAfarIndex] = ParallelClustering.RunningSolution.LocalCreatedIndex[LocalClusterRealIndixes[FromAfarIndex]];
                }
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<MPIPacket<int>>(ref FromAfarInt, SourceNode);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                for (int FromAfarIndex = 0; FromAfarIndex < FromAfarInt.NumberofPoints; FromAfarIndex++)
                {
                    int CreatedIndex = FromAfarInt.Marray[FromAfarIndex];
                    if (ClusteringSolution.UniversalMapping[CreatedIndex] == null)
                        ClusteringSolution.UniversalMapping[CreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, 0);

                    UniversalMapping[CreatedIndex].GlobalClusterNumber = CountGlobal;
                    TotalClusterSummary.CreatedIndex[CountGlobal] = CreatedIndex;
                    TotalClusterSummary.CurrentNode[CountGlobal] = SourceNode;
                    ++CountGlobal;
                }

                //  Transport Occupation Counts
                if (SourceNode == DAVectorUtility.MPI_Rank)
                {
                    FromAfarInt.NumberofPoints = CountLocal;
                    FromAfarInt.FirstPoint = 0;
                    for (int FromAfarIndex = 0; FromAfarIndex < CountLocal; FromAfarIndex++)
                        FromAfarInt.Marray[FromAfarIndex] = ParallelClustering.RunningSolution.OccupationCounts_k_[LocalClusterRealIndixes[FromAfarIndex]];
                }
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<MPIPacket<int>>(ref FromAfarInt, SourceNode);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                for (int FromAfarIndex = 0; FromAfarIndex < FromAfarInt.NumberofPoints; FromAfarIndex++)
                    TotalClusterSummary.OccupationCount[SaveCountGlobal + FromAfarIndex] = FromAfarInt.Marray[FromAfarIndex];

                //  Transport Center positions
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {  
                    if (SourceNode == DAVectorUtility.MPI_Rank)
                    {
                        FromAfarDouble.NumberofPoints = CountLocal;
                        FromAfarDouble.FirstPoint = 0;
                        for (int FromAfarIndex = 0; FromAfarIndex < CountLocal; FromAfarIndex++)
                            FromAfarDouble.Marray[FromAfarIndex] = ParallelClustering.RunningSolution.Y_k_i_[LocalClusterRealIndixes[FromAfarIndex]][VectorIndex];
                    }
                    DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                    DAVectorUtility.MPI_communicator.Broadcast<MPIPacket<double>>(ref FromAfarDouble, SourceNode);
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                    for (int FromAfarIndex = 0; FromAfarIndex < FromAfarDouble.NumberofPoints; FromAfarIndex++)
                        TotalClusterSummary.CenterPosition[SaveCountGlobal + FromAfarIndex][VectorIndex] = FromAfarDouble.Marray[FromAfarIndex];
                }
            }
            TotalClusterSummary.NumberofCenters = CountGlobal;

        }   // End SetGlobalClusterNumbers()

        //    Return Pointer associate with given Active Cluster Index associated with given point
        //    Return -1 if Point not associated with this cluster
        public int MapClusterToIndirect(int LocalPointPosition, int ActiveClusterIndex)
        {
            int CreatedIndex = ClusteringSolution.MapActivetoCreatedIndex(ActiveClusterIndex, this);

            int IndirectSize = this.NumClusters_alpha_[LocalPointPosition];
            for(int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
            {
                if(CreatedIndex == this.Map_alpha_PointertoCreatedIndex[LocalPointPosition][IndirectClusterIndex])
                    return IndirectClusterIndex;
            }
            return -1;

        }   // End  MapClusterToIndirect(int LocalPointPosition)

        //   Set Clusters for a given point. Tries to preserve existing Malpha -- not always good idea 
        public void SetClustersforaPoint(int LocalPointPosition)
        {

            int IndirectClusterSize = this.NumClusters_alpha_[LocalPointPosition];

//          Following test removed as cases where clusters removed need Malpha reset even if Clusters perfect
            if ((!DistributedExecutionMode) && (IndirectClusterSize == this.Ncent_ThisNode) && (IndirectClusterSize <= TargetClustersperPoint))
               return;

            int TargetNumber1 = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, TargetClustersperPoint);
            int TargetNumber = TargetNumber1;   // Number of Target non-sponge clusters
            if(this.SpongeCluster >= 0 )
                TargetNumber--;

            // DAVectorUtility.StartSubTimer(13); // Wrong as per thread
            int[] ProposedClusters = new int[TargetNumber];
            double[] ClusterDistances = new double[TargetNumber];

            int PositionofWorstCluster = -1;
            for(int listindex = 0; listindex < TargetNumber; listindex++)
            {
                ProposedClusters[listindex] = -1;
                ClusterDistances[listindex] = 0.0;
            }
            int NearestClusterIndex = -1;    // Cluster Index of Nearest Cluster
            double NearestClusterDistance = 0.0;
            int ActiveSpongeIndex = -1;
            if (this.SpongeCluster >= 0)
                ActiveSpongeIndex = ClusteringSolution.ActiveClusterIndices[this.SpongeCluster];

            //  Note ActiveClusterIndex runs over local and remote clusters
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberAvailableActiveClusters; ActiveClusterIndex++)
            {
                if (ActiveClusterIndex == ActiveSpongeIndex)
                    continue;   // We will add back Sponge Cluster later

                double CurrentClusterScaledDistance = DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(LocalPointPosition, ActiveClusterIndex, this);
                if( (NearestClusterIndex < 0) || (CurrentClusterScaledDistance < NearestClusterDistance))
                {
                    NearestClusterIndex = ActiveClusterIndex;
                    NearestClusterDistance = CurrentClusterScaledDistance;
                }
                GlobalReductions.FindMinimumSet(CurrentClusterScaledDistance, ActiveClusterIndex, ref PositionofWorstCluster, ClusterDistances, ProposedClusters, TargetNumber);
            }

            //  Make of clusters in order -- starting with sponge if exists
            int[] ClusterIndicesTobeUsed = new int[TargetNumber1];
            int NumberofClustersTobeUsed = 0;
            if (ActiveSpongeIndex >= 0)
            {
                ClusterIndicesTobeUsed[0] = ActiveSpongeIndex;
                NumberofClustersTobeUsed = 1;
            }

            //  Add in Non sponge clusters
            int Initial_ActualNumberofClusters = NumberofClustersTobeUsed;
            for(int LoopOverMinimumSet = 0; LoopOverMinimumSet < TargetNumber; LoopOverMinimumSet++)
            {
                int ActiveClusterIndex = ProposedClusters[LoopOverMinimumSet];
                if(ActiveClusterIndex < 0 )
                    continue;
                if ((ClusteringSolution.ExponentArgumentCut > 0.0) && ( (ClusterDistances[LoopOverMinimumSet] - NearestClusterDistance) > (2.0 * ExponentArgumentCut * this.Temperature) ))
                    continue;
                ClusterIndicesTobeUsed[NumberofClustersTobeUsed] = ActiveClusterIndex;
                NumberofClustersTobeUsed++;
            }

            // Add Nearest Cluster if none (other than sponge) selected
            if(NumberofClustersTobeUsed == Initial_ActualNumberofClusters)
            {
                if ((NumberofClustersTobeUsed < 0) || (NumberofClustersTobeUsed >= TargetNumber1) || (NumberofClustersTobeUsed >= ClusterIndicesTobeUsed.Length))
                {
                    Exception e = DAVectorUtility.SALSAError("Error " + NumberofClustersTobeUsed.ToString() + " " + TargetNumber1.ToString()
                         + " Old Number " + IndirectClusterSize.ToString() + " " + Initial_ActualNumberofClusters.ToString() + " " + TargetClustersperPoint.ToString()
                         + " " + ClusterIndicesTobeUsed.Length + " " + ClusteringSolution.NumberAvailableActiveClusters.ToString() 
                         + DistributedClusteringSolution.StorageforTransportedClusters.SizeOfTransportedArray.ToString());
                    throw e;
                }
                ClusterIndicesTobeUsed[NumberofClustersTobeUsed] = NearestClusterIndex;
                ++NumberofClustersTobeUsed;
            }

            //  Now make list for storing back
            double[] Mvalues = new double[NumberofClustersTobeUsed];
            double Msum = 0.0;
            int NearestClusterIndirectIndex = -1;
            for(int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
            {
                Mvalues[IndirectClusterIndex] = 0.0;
                int ActiveClusterIndex = ClusterIndicesTobeUsed[IndirectClusterIndex];
                if( ActiveClusterIndex == NearestClusterIndex )
                    NearestClusterIndirectIndex = IndirectClusterIndex;
                int OldIndirectClusterIndex = this.MapClusterToIndirect(LocalPointPosition, ActiveClusterIndex);
                if(OldIndirectClusterIndex < 0)
                    continue;
                Mvalues[IndirectClusterIndex] = this.M_alpha_kpointer_[LocalPointPosition][OldIndirectClusterIndex];
                Msum += Mvalues[IndirectClusterIndex];
            }
            if(Msum < 0.2)
            {
                for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                    this.M_alpha_kpointer_[LocalPointPosition][IndirectClusterIndex] = 1.0 / (double)NumberofClustersTobeUsed;
                Msum = 1.0;
            }

            this.NumClusters_alpha_[LocalPointPosition] = NumberofClustersTobeUsed;
            if (NumberofClustersTobeUsed <= 0)
            {
                Exception e = DAVectorUtility.SALSAError("Error due to zero New Number of Clusters Point " + 
                    (LocalPointPosition + DAVectorUtility.PointStart_Process).ToString() + " Old Number " + IndirectClusterSize.ToString());
                throw e;
            }

            for(int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
            {
                this.M_alpha_kpointer_[LocalPointPosition][IndirectClusterIndex] = Mvalues[IndirectClusterIndex] / Msum;
                this.Map_alpha_PointertoCreatedIndex[LocalPointPosition][IndirectClusterIndex] = ClusteringSolution.MapActivetoCreatedIndex(ClusterIndicesTobeUsed[IndirectClusterIndex], this);
            }

            this.DiffMalpha_k_Set = -1;
            // DAVectorUtility.StopSubTimer(13);  // Wrong as per thread
            return;

        }   // End SetClustersforaPoint(int LocalPointPosition)

        public static void CopySolution(ClusteringSolution From, ClusteringSolution To)
        {
            To.DistributedExecutionMode = From.DistributedExecutionMode;
            To.SpongeCluster = From.SpongeCluster;
            To.PairwiseHammy = From.PairwiseHammy; // Value of Hamiltonian
            To.OldHammy = From.OldHammy;  //  Previous value of Hamiltonian
            To.Ncent_ThisNode = From.Ncent_ThisNode;    //the current number of clusters in node
            To.Ncent_Global = From.Ncent_Global;    //the current total number of clusters
            To.ActualCoolingFactor = From.ActualCoolingFactor;  // Actual Cooling Factor
            To.Temperature = From.Temperature; //the current temperature
            To.Eigenvectorset = From.Eigenvectorset;
            To.SolutionSet = From.SolutionSet;
            To.IterationSetAt = From.IterationSetAt;
            To.CorrelationsSet = From.CorrelationsSet;
            To.DiffMalpha_k_Set = From.DiffMalpha_k_Set;
            To.YPreviousSet = From.YPreviousSet;
            To.ClustertoSplit = From.ClustertoSplit;
            To.TotaloverVectorIndicesAverageWidth = From.TotaloverVectorIndicesAverageWidth;
            for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                To.AverageWidth[VectorIndex1] = From.AverageWidth[VectorIndex1];

            int NumberClusters = From.Ncent_ThisNode;
            for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
            {
                To.C_k_[ClusterIndex] = From.C_k_[ClusterIndex];
                To.ClusterScaledSquaredWidth_k_[ClusterIndex] = From.ClusterScaledSquaredWidth_k_[ClusterIndex];
                To.P_k_[ClusterIndex] = From.P_k_[ClusterIndex];
                To.FreezingMeasure_k_[ClusterIndex] = From.FreezingMeasure_k_[ClusterIndex];
                To.Splittable_k_[ClusterIndex] = From.Splittable_k_[ClusterIndex];
                To.SplitPriority_k_[ClusterIndex] = From.SplitPriority_k_[ClusterIndex];
                To.Eigenvalue_k[ClusterIndex] = From.Eigenvalue_k[ClusterIndex];
                To.DiffMsummed_k_[ClusterIndex] = From.DiffMsummed_k_[ClusterIndex];
                To.OccupationCounts_k_[ClusterIndex] = From.OccupationCounts_k_[ClusterIndex];
                To.LocalCreatedIndex[ClusterIndex] = From.LocalCreatedIndex[ClusterIndex];  
                To.LocalStatus[ClusterIndex] = From.LocalStatus[ClusterIndex];   
                To.LocalSplitCreatedIndex[ClusterIndex] = From.LocalSplitCreatedIndex[ClusterIndex];
                To.YPreviousActuallySet[ClusterIndex] = From.YPreviousActuallySet[ClusterIndex];

                for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                {
                    To.Y_k_i_[ClusterIndex][VectorIndex1] = From.Y_k_i_[ClusterIndex][VectorIndex1];
                    To.YPrevious_k_i_[ClusterIndex][VectorIndex1] = From.YPrevious_k_i_[ClusterIndex][VectorIndex1];
                    To.Eigenvector_k_i[ClusterIndex][VectorIndex1] = From.Eigenvector_k_i[ClusterIndex][VectorIndex1];
                    To.DC_k_DY_k_i_[ClusterIndex][VectorIndex1] = From.DC_k_DY_k_i_[ClusterIndex][VectorIndex1];
                    To.Sigma_k_i_[ClusterIndex][VectorIndex1] = From.Sigma_k_i_[ClusterIndex][VectorIndex1];

                    if (From.CorrelationsSet)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                        {
                            From.Correlation_k_i_j[ClusterIndex][VectorIndex1, VectorIndex2]
                                = To.Correlation_k_i_j[ClusterIndex][VectorIndex1, VectorIndex2];
                        }
                    }
                }
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = From.NumClusters_alpha_[alpha];
                    To.NumClusters_alpha_[alpha] = IndirectSize;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        To.M_alpha_kpointer_[alpha][IndirectClusterIndex] = From.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        To.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = From.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                    }
                }
            }); // End loop copying Point dependent quantities
        }   // End CopySolution

        public void CleanupClusters()
        {   // Must follow by SetActiveClusters()

            int NewClusterIndex = 0;
            this.ClustertoSplit = -1;
            int[] ClusterMapping = new int[this.Ncent_ThisNode];

            for (int OldClusterIndex = 0; OldClusterIndex < this.Ncent_ThisNode; OldClusterIndex++)
            {
                int status = this.LocalStatus[OldClusterIndex];
                if ((status < 0) || (status > 2))
                {
                    if (this.SpongeCluster == OldClusterIndex)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Attempt to Remove Sponge Cluster " + OldClusterIndex.ToString());
                        throw (e);
                    }
                    ClusterMapping[OldClusterIndex] = -1;
                    continue;
                }

                ClusterMapping[OldClusterIndex] = NewClusterIndex;
                if (this.SpongeCluster == OldClusterIndex)
                    this.SpongeCluster = NewClusterIndex;
                this.DiffMsummed_k_[NewClusterIndex] = this.DiffMsummed_k_[OldClusterIndex];
                this.C_k_[NewClusterIndex] = this.C_k_[OldClusterIndex];
                this.ClusterScaledSquaredWidth_k_[NewClusterIndex] = this.ClusterScaledSquaredWidth_k_[OldClusterIndex];
                this.P_k_[NewClusterIndex] = this.P_k_[OldClusterIndex];
                this.FreezingMeasure_k_[NewClusterIndex] = this.FreezingMeasure_k_[OldClusterIndex];
                this.Splittable_k_[NewClusterIndex] = this.Splittable_k_[OldClusterIndex];
                this.SplitPriority_k_[NewClusterIndex] = this.SplitPriority_k_[OldClusterIndex];
                this.Eigenvalue_k[NewClusterIndex] = this.Eigenvalue_k[OldClusterIndex];
                this.OccupationCounts_k_[NewClusterIndex] = this.OccupationCounts_k_[OldClusterIndex];
                this.LocalCreatedIndex[NewClusterIndex] = this.LocalCreatedIndex[OldClusterIndex];
                this.LocalStatus[NewClusterIndex] = this.LocalStatus[OldClusterIndex];
                this.LocalSplitCreatedIndex[NewClusterIndex] = this.LocalSplitCreatedIndex[OldClusterIndex];
                this.YPreviousActuallySet[NewClusterIndex] = this.YPreviousActuallySet[OldClusterIndex];

                for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                {
                    this.Y_k_i_[NewClusterIndex][VectorIndex1] = this.Y_k_i_[OldClusterIndex][VectorIndex1];
                    this.YPrevious_k_i_[NewClusterIndex][VectorIndex1] = this.YPrevious_k_i_[OldClusterIndex][VectorIndex1];
                    this.Sigma_k_i_[NewClusterIndex][VectorIndex1] = this.Sigma_k_i_[OldClusterIndex][VectorIndex1];
                    this.Eigenvector_k_i[NewClusterIndex][VectorIndex1] = this.Eigenvector_k_i[OldClusterIndex][VectorIndex1];
                    this.DC_k_DY_k_i_[NewClusterIndex][VectorIndex1] = this.DC_k_DY_k_i_[OldClusterIndex][VectorIndex1];

                    if (this.CorrelationsSet)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                        {
                            this.Correlation_k_i_j[NewClusterIndex][VectorIndex1, VectorIndex2]
                                = this.Correlation_k_i_j[OldClusterIndex][VectorIndex1, VectorIndex2];
                        }
                    }
                }
                int CreatedIndex = this.LocalCreatedIndex[OldClusterIndex];
                UniversalMapping[CreatedIndex].Availability = 1 + NewClusterIndex;
                ++NewClusterIndex;
            }
            this.Ncent_ThisNode = NewClusterIndex;

        }   // End CleanupClusters()

        public void RemoveCluster( int RemovedIndex)
        {
            //  Process Child and Parent Clusters before we change numbers
            if (this.DistributedExecutionMode)
            {
                this.LocalStatus[RemovedIndex] = -1;
                return;
            }
            this.DiffMalpha_k_Set = -1;

            // Reduce Number of Clusters
            //  Global Ncent set later in call to SetActiveClusters
            --this.Ncent_ThisNode;
            if (this.SpongeCluster >= 0)
            {
                if (this.SpongeCluster == RemovedIndex)
                {
                    Exception e = DAVectorUtility.SALSAError(" Attempt to Remove Sponge Cluster " + RemovedIndex.ToString());
                    throw (e);
                }
                if (this.SpongeCluster > RemovedIndex)
                    --this.SpongeCluster;
            }
            int DeletedCreatedIndex = this.LocalCreatedIndex[RemovedIndex];
            UniversalMapping[DeletedCreatedIndex].Availability = 0;
            this.ClustertoSplit = -1;
            

            for (int ClusterIndex = RemovedIndex; ClusterIndex < this.Ncent_ThisNode; ClusterIndex++)
            {
                this.DiffMsummed_k_[ClusterIndex] = this.DiffMsummed_k_[ClusterIndex+1];
                this.C_k_[ClusterIndex] = this.C_k_[ClusterIndex+1];
                this.ClusterScaledSquaredWidth_k_[ClusterIndex] = this.ClusterScaledSquaredWidth_k_[ClusterIndex+1];
                this.P_k_[ClusterIndex] = this.P_k_[ClusterIndex+1];
                this.FreezingMeasure_k_[ClusterIndex] = this.FreezingMeasure_k_[ClusterIndex+1];
                this.Splittable_k_[ClusterIndex] = this.Splittable_k_[ClusterIndex+1];
                this.SplitPriority_k_[ClusterIndex] = this.SplitPriority_k_[ClusterIndex + 1];
                this.Eigenvalue_k[ClusterIndex] = this.Eigenvalue_k[ClusterIndex+1];
                this.OccupationCounts_k_[ClusterIndex] = this.OccupationCounts_k_[ClusterIndex + 1];
                int CurrentCreatedIndex = this.LocalCreatedIndex[ClusterIndex + 1];
                this.LocalCreatedIndex[ClusterIndex] = CurrentCreatedIndex;
                UniversalMapping[CurrentCreatedIndex].Availability = 1 + ClusterIndex;
                this.LocalStatus[ClusterIndex] = this.LocalStatus[ClusterIndex + 1];
                this.LocalSplitCreatedIndex[ClusterIndex] = this.LocalSplitCreatedIndex[ClusterIndex + 1];
                this.YPreviousActuallySet[ClusterIndex] = this.YPreviousActuallySet[ClusterIndex + 1];

                for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                {
                    this.Y_k_i_[ClusterIndex][VectorIndex1] = this.Y_k_i_[ClusterIndex+1][VectorIndex1];
                    this.YPrevious_k_i_[ClusterIndex][VectorIndex1] = this.YPrevious_k_i_[ClusterIndex + 1][VectorIndex1];
                    this.Sigma_k_i_[ClusterIndex][VectorIndex1] = this.Sigma_k_i_[ClusterIndex+1][VectorIndex1];
                    this.Eigenvector_k_i[ClusterIndex][VectorIndex1] = this.Eigenvector_k_i[ClusterIndex + 1][VectorIndex1];
                    this.DC_k_DY_k_i_[ClusterIndex][VectorIndex1] = this.DC_k_DY_k_i_[ClusterIndex + 1][VectorIndex1];

                    if (this.CorrelationsSet)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                        {
                            this.Correlation_k_i_j[ClusterIndex][VectorIndex1, VectorIndex2]
                                = this.Correlation_k_i_j[ClusterIndex + 1][VectorIndex1, VectorIndex2];
                        }
                    }
                }
            }

            //  Record changes of count in ActiveCluster and Ncent_Global
            this.YPreviousActuallySet[this.Ncent_ThisNode] = false;
            ParallelClustering.RunningSolution.SetActiveClusters();

            if (Program.UseTriangleInequality_DA > 0)
                DATriangleInequality.DeleteCenter(RemovedIndex, DeletedCreatedIndex);

            GlobalReductions.FindDoubleArraySum OldClusterNumberHistogram = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, this.Ncent_ThisNode + 2);
            GlobalReductions.FindDoubleArraySum NewClusterNumberHistogram = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, this.Ncent_ThisNode + 1);
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                if (Program.RemovalDiagnosticPrint)
                {
                    OldClusterNumberHistogram.startthread(ThreadNo);
                    NewClusterNumberHistogram.startthread(ThreadNo);
                }
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = this.NumClusters_alpha_[alpha];
                    int decrement = 0;
                    Double Msum = 0.0;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if (CreatedIndex == DeletedCreatedIndex)
                        {
                            if (decrement > 0)
                            {
                                Exception e = DAVectorUtility.SALSAError(" Double Cluster Entry " + RemovedIndex.ToString()
                                    + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                                throw (e);
                            }
                            decrement = 1;
                            continue;
                        }
                        if( CreatedIndex <= 0)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Illegal Created Index in Point List " + CreatedIndex.ToString() + " " + IndirectClusterIndex.ToString() 
                                + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                            throw (e);
                        }

                        int listmember = UniversalMapping[CreatedIndex].Availability - 1;
                        if( UniversalMapping[CreatedIndex].IterationSet != CurrentIteration)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Out of Date Cluster in Point List " + listmember.ToString() + " Iterations " + UniversalMapping[CreatedIndex].IterationSet.ToString() + " " + CurrentIteration.ToString()
                                + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                            throw (e);
                        }

                        if (listmember >= this.Ncent_ThisNode)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Bad List Member " + listmember.ToString() + " Position " + IndirectClusterIndex.ToString() + " Decrement " + decrement.ToString()
                                 + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                            throw (e);
                        }

                        if (decrement == 1)
                        {
                            this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex - decrement] = CreatedIndex; 
                            this.M_alpha_kpointer_[alpha][IndirectClusterIndex - decrement] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        }
                        Msum += this.M_alpha_kpointer_[alpha][IndirectClusterIndex - decrement];

                    }
                    this.NumClusters_alpha_[alpha] = IndirectSize - decrement;
                    int NewIndirectSize = this.NumClusters_alpha_[alpha];
                    if ( (NewIndirectSize <= 0) && (Program.UseTriangleInequality_DA == 0))
                    {
                        Exception e = DAVectorUtility.SALSAError(" Zero Cluster Count after removing cluster for Point " + (alpha + DAVectorUtility.PointStart_Process).ToString()
                             + " Originally " + IndirectSize.ToString());
                        throw (e);
                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                    {
                        int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if (CreatedIndex <= 0)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Illegal Created Index in Corrected Point List " + CreatedIndex.ToString() + " " + IndirectClusterIndex.ToString()
                                + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                            throw (e);
                        }
                        int listmember = UniversalMapping[CreatedIndex].Availability - 1;
                        if ( (listmember >= this.Ncent_ThisNode) || (listmember < 0) )
                        {
                            Exception e = DAVectorUtility.SALSAError(" Bad Corrected List Member " + listmember.ToString() + " Position " + IndirectClusterIndex.ToString() 
                                 + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                            throw (e);
                        }
                    }
                    //  If Program.UseTriangleInequality_DA > 0, zero cluster count will be fixed in DistributedClusteringSolution.ManageMajorSynchronization(true) called later
                    if (this.NumClusters_alpha_[alpha] < 1)
                    {
                        if( Program.UseTriangleInequality_DA == 0)
                            this.SetClustersforaPoint(alpha);
                    }
                    else
                    {
                        if (Msum < 0.2)
                        {
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                                this.M_alpha_kpointer_[alpha][IndirectClusterIndex] = 1.0 / (double)NewIndirectSize;
                        }
                        else
                        {
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                                this.M_alpha_kpointer_[alpha][IndirectClusterIndex] = this.M_alpha_kpointer_[alpha][IndirectClusterIndex] / Msum;
                        }

                    }

                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                        ParallelClustering.RunningSolution.LegalCluster(alpha, IndirectClusterIndex);

                    if (Program.RemovalDiagnosticPrint)
                    {
                        OldClusterNumberHistogram.addapoint(ThreadNo, IndirectSize);
                        NewClusterNumberHistogram.addapoint(ThreadNo, NewIndirectSize);
                    }
                }
            }); // End loop copying Point dependent quantities

            if (Program.RemovalDiagnosticPrint)
            {
                OldClusterNumberHistogram.sumoverthreadsandmpi();
                NewClusterNumberHistogram.sumoverthreadsandmpi();

                string message = "\nOld Counts ";
                int ArraySize = 2 + this.Ncent_ThisNode;
                for (int ClusterIndex = 0; ClusterIndex < ArraySize; ClusterIndex++)
                    message += OldClusterNumberHistogram.TotalSum[ClusterIndex].ToString() + " ";
                ArraySize--;
                message += "\nNew Counts ";
                for (int ClusterIndex = 0; ClusterIndex < ArraySize; ClusterIndex++)
                    message += NewClusterNumberHistogram.TotalSum[ClusterIndex].ToString() + " ";
                DAVectorUtility.SALSAPrint(1, "Cluster " + RemovedIndex.ToString() + " Removed" + message);
            }
            return;

        }   // End RemoveCluster

        public bool LegalCluster(int alpha, int IndirectClusterIndex)
        {
            int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha]; 
            int CreatedIndex = this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
            int MappedClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
            int ClusterIterationNo = ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet;
            if ((ClusterIterationNo < ClusteringSolution.CurrentIteration) || (MappedClusterIndex == 0) )
            {
                string errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString() 
                    + " Created Index " + CreatedIndex.ToString() + " Actual Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Cluster Iteration " + ClusterIterationNo.ToString()
                    + " Mapped Index " + MappedClusterIndex.ToString() + " Full Set Created Indices ";
                for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
                    errormessage += ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop].ToString() + " ";
                Exception e = DAVectorUtility.SALSAError(errormessage);
                throw (e);
            }
            double Mvalue = this.M_alpha_kpointer_[alpha][IndirectClusterIndex];
            if (Double.IsNaN(Mvalue) || Double.IsInfinity(Mvalue))
            {
                string message = "";
                for (int loop = 0; loop < IndirectSize; loop++)
                    message += this.Map_alpha_PointertoCreatedIndex[alpha][loop].ToString() + " Malpha " + ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][loop].ToString("F5") + " * ";
                Exception e = DAVectorUtility.SALSAError("Arithmetic Error Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Number of Clusters " + IndirectSize.ToString()
                    + " Indirect Index " + IndirectClusterIndex.ToString() + " Created Index " + this.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex].ToString()
                    + "\n" + message);
                throw (e);
            }
            return true;

        } // End LegalCluster

    }   // End ClusteringSolution

    //  Data structure to record quality of clusters
    //  SpongePointsinCut Number of Sponge Points within Distance Cut
    //  NumberClustersinCut Number of Clusters within Distance Cut
    //  NearbyClusterIndices Labels of Nearby Clusters sorted by Distance
    //  NearbyClusterDistances Distance of Nearby Clusters
    public class ClusterQuality
    {
        int NumberClusters = 0;
        int SpongeOption = -1;
        int NumberNearbyClusters = 0;
        double NearbyCut = 0.0;
        double TotalConfusedClusterPoints = 0.0;
        double TotalConfusedSpongePoints1 = 0.0;
        double TotalConfusedSpongePoints2 = 0.0;
        double AverageClustersinCut = 0.0;
        int TotalIllegalPoints = 0;
        int HistogramOccupationMax = 202;
        int HistogramDistancesMax = 200;

        int[] SpongePointsinCut;
        int[] NonSpongePointsinCut;
        int[] NumberClustersinCut;
        double[] ClusterMaxWidths;
        int[] ClusterIllegalPoints;
        int[][] NearbyClusterIndices;
        double[][] NearbyClusterDistances;
        int[] HistogramOccupationValues;
        int[] HistogramDistancesValues;

        //  InputNumberClusters -- Number of Clusters
        //  InputSpongeOption -- Index of Sponge or if < 0 No sponge
        //  InputNumberNearbyClusters -- List this number of nearby clusters
        //  InputNearbyCut -- Count Sponge Points and Clusters within this cut times sigma
        public ClusterQuality(int InputNumberClusters, int InputSpongeOption, int InputNumberNearbyClusters, double InputNearbyCut)
        {
            NumberClusters = InputNumberClusters;
            SpongeOption = InputSpongeOption;
            NumberNearbyClusters = Math.Min(InputNumberNearbyClusters, NumberClusters-1);
            NearbyCut = InputNearbyCut * InputNearbyCut;

            HistogramOccupationValues = new int[HistogramOccupationMax];
            HistogramDistancesValues = new int[HistogramDistancesMax];

            SpongePointsinCut = new int[NumberClusters];
            NonSpongePointsinCut = new int[NumberClusters];
            NumberClustersinCut = new int[NumberClusters];
            NearbyClusterIndices = new int[NumberClusters][];
            ClusterIllegalPoints = new int[NumberClusters];
            ClusterMaxWidths = new double[NumberClusters];
            NearbyClusterDistances = new double[NumberClusters][];

            if (InputNumberNearbyClusters > 0)
            {
                for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                {
                    NearbyClusterIndices[ClusterIndex] = new int[NumberNearbyClusters];
                    NearbyClusterDistances[ClusterIndex] = new double[NumberNearbyClusters];
                }
                if (SpongeOption >= 0)
                {
                    for (int spongespecial = 0; spongespecial < NumberNearbyClusters; spongespecial++)
                    {
                        NearbyClusterIndices[SpongeOption][spongespecial] = 0;
                        NearbyClusterDistances[SpongeOption][spongespecial] = 0.0;
                    }
                    NumberClustersinCut[SpongeOption] = 0;
                }
            }
        }   // End Initializer ClusterQuality(int InputNumberClusters, int InputSpongeOption, int InputNumberNearbyClusters, double InputNearbyCut)

        public void SetClusterStatistics()
        {
            GlobalReductions.FindDoubleArraySum PointsperCluster = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, HistogramOccupationMax);
            GlobalReductions.FindDoubleArraySum DistancesfromCluster = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, HistogramDistancesMax);
            GlobalReductions.FindIntSum TotalClustersinCut = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);

            double DistanceCut = (double) HistogramDistancesMax;

            Range[] ParallelClusterRange = RangePartitioner.Partition(NumberClusters, DAVectorUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ClusterThreadNo) =>
            {
                PointsperCluster.startthread(ClusterThreadNo);
                DistancesfromCluster.startthread(ClusterThreadNo);
                if (DAVectorUtility.MPI_Rank == 0)
                {
                    int beginindex = ParallelClusterRange[ClusterThreadNo].StartIndex;
                    int indexlength = ParallelClusterRange[ClusterThreadNo].Length;

                    double[] ClusterSigma = new double[Program.ParameterVectorDimension];
                    for (int ClusterIndex1 = beginindex; ClusterIndex1 < beginindex + indexlength; ClusterIndex1++)
                    {
                        int count = 0;
                        if (ClusterIndex1 == SpongeOption)
                            continue;

                        int histcount = ClusteringSolution.TotalClusterSummary.OccupationCount[ClusterIndex1];
                        if (histcount > (HistogramOccupationMax - 1))
                            histcount = HistogramOccupationMax - 1;
                        PointsperCluster.addapoint(ClusterThreadNo, histcount);

                        for (int ClusterIndex2 = 0; ClusterIndex2 < NumberClusters; ClusterIndex2++)
                        {
                            if (ClusterIndex2 == SpongeOption)
                                continue;
                            if (ClusterIndex1 == ClusterIndex2)
                                continue;
                            Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1], ref ClusterSigma);
                            double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                                ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex2], ClusterSigma);

                            if (tmp <= DistanceCut)
                            {
                                int itmp = (int)Math.Floor(tmp + 0.5);
                                itmp = Math.Min(itmp, HistogramDistancesMax - 1);
                                DistancesfromCluster.addapoint(ClusterThreadNo, itmp);
                            }
                            if (tmp < NearbyCut)
                                ++count;
                        }
                        NumberClustersinCut[ClusterIndex1] = count;
                        TotalClustersinCut.addapoint(ClusterThreadNo, count);
                    }
                }
            });

            TotalClustersinCut.sumoverthreadsandmpi();
            double ClusterSum = NumberClusters;
            if (SpongeOption >= 0)
                ClusterSum -= 1.0;
            AverageClustersinCut = (double) TotalClustersinCut.TotalInt / (2.0 * ClusterSum);

            PointsperCluster.sumoverthreadsandmpi();
            DistancesfromCluster.sumoverthreadsandmpi();
            for (int histloop = 0; histloop < HistogramOccupationMax; histloop++)
                HistogramOccupationValues[histloop] = (int) (PointsperCluster.TotalSum[histloop] + 0.001);
            for (int histloop = 0; histloop < HistogramDistancesMax; histloop++)
                HistogramDistancesValues[histloop] = (int)(DistancesfromCluster.TotalSum[histloop] + 0.001);

        }   // End SetClusterStatistics()

        public void SetPointStatistics()
        {
            TotalIllegalPoints = 0;
            
            GlobalReductions.FindDoubleArraySum AccumulateSpongePoints = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, NumberClusters);
            GlobalReductions.FindDoubleArraySum AccumulateNonSpongePoints = new GlobalReductions.FindDoubleArraySum(DAVectorUtility.ThreadCount, NumberClusters);
            GlobalReductions.FindDoubleSum ConfusedSpongePoints1 = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleSum ConfusedSpongePoints2 = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleSum ConfusedClusterPoints = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleMax[] FindMaxExtension = new GlobalReductions.FindDoubleMax[NumberClusters];
            GlobalReductions.FindIntSum[] FindIllegalPoints = new GlobalReductions.FindIntSum[NumberClusters];
            for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
            {
                FindMaxExtension[ClusterIndex] = new GlobalReductions.FindDoubleMax(DAVectorUtility.ThreadCount);
                FindIllegalPoints[ClusterIndex] = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);
            }
            
            //  Parallel Section over points
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] ClusterSigma = new double[Program.ParameterVectorDimension]; 

                AccumulateSpongePoints.startthread((int)ThreadNo);
                AccumulateNonSpongePoints.startthread((int)ThreadNo);

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;

                for (int index = beginpoint; index < indexlen + beginpoint; index++)
                {
                    int AssignedCluster = Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process];
                    bool isSponge = false;
                    if (AssignedCluster == SpongeOption)
                        isSponge = true;

                    //  Remove Sponge Confused Points
                    if (isSponge)
                    {
                        int nearbyrealcluster = -1;
                        double NearesetDistce = -1.0;
                        for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                        {
                            if (ClusterIndex == SpongeOption)
                                continue;
                            if (ClusterIndex == AssignedCluster)
                                continue;
                            Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ref ClusterSigma);
                            double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(Program.PointPosition[index], ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ClusterSigma);
                            if (tmp < NearbyCut)
                            {
                                if ( (nearbyrealcluster < 0) || (tmp < NearesetDistce))
                                {
                                    nearbyrealcluster = ClusterIndex;
                                    NearesetDistce = tmp;
                                }
                            }
                        }
                        if (nearbyrealcluster >= 0)
                        {
                            Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process] = nearbyrealcluster;
                            AssignedCluster = nearbyrealcluster;
                            isSponge = false;
                            ConfusedSpongePoints1.addapoint((int)ThreadNo, 1.0);
                        }
                    }

                    //  Not in Sponge
                    if (!isSponge)
                    {
                        if ((AssignedCluster < 0) || (AssignedCluster >= NumberClusters))
                        {
                            Exception e = DAVectorUtility.SALSAError("Point " + (index + DAVectorUtility.PointStart_Process).ToString()
                                + " Assigned " + AssignedCluster + " Max " + NumberClusters.ToString());
                            throw (e);
                        }
                        Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[AssignedCluster], ref ClusterSigma);
                        double distce = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(Program.PointPosition[index], ClusteringSolution.TotalClusterSummary.CenterPosition[AssignedCluster], ClusterSigma);
                        FindMaxExtension[AssignedCluster].addapoint(ThreadNo, distce);
                        if (distce > NearbyCut)
                        {
                            FindIllegalPoints[AssignedCluster].addapoint(ThreadNo, 1);
                            if (SpongeOption >= 0)
                            {
                                AssignedCluster = SpongeOption;
                                Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process] = AssignedCluster;
                            }
                        }
                    }
                    bool SpongeConfused = false;
                    bool ClusterConfused = false;
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                    {
                        if( ClusterIndex == SpongeOption)
                            continue;
                        if (ClusterIndex == AssignedCluster)
                            continue;
                        Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ref ClusterSigma);
                        double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(Program.PointPosition[index], ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex], ClusterSigma);
                        if (tmp < NearbyCut)
                        {
                            if (isSponge)
                            {   // Impossible to reach
                                AccumulateSpongePoints.addapoint((int)ThreadNo, ClusterIndex);
                                SpongeConfused = true;
                            }
                            else
                            {
                                AccumulateNonSpongePoints.addapoint((int)ThreadNo, ClusterIndex);
                                ClusterConfused = true;
                            }
                        }
                    }
                    if(SpongeConfused)
                        ConfusedSpongePoints2.addapoint((int)ThreadNo, 1.0);
                    if(ClusterConfused)
                        ConfusedClusterPoints.addapoint((int)ThreadNo, 1.0);
                }
            });  // End Parallel Section

            AccumulateSpongePoints.sumoverthreadsandmpi();
            AccumulateNonSpongePoints.sumoverthreadsandmpi();
            ConfusedClusterPoints.sumoverthreadsandmpi();
            ConfusedSpongePoints1.sumoverthreadsandmpi();
            ConfusedSpongePoints2.sumoverthreadsandmpi();

            for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
            {
                FindIllegalPoints[ClusterIndex].sumoverthreadsandmpi();
                FindMaxExtension[ClusterIndex].sumoverthreadsandmpi();
                SpongePointsinCut[ClusterIndex] = (int)Math.Floor(AccumulateSpongePoints.TotalSum[ClusterIndex] + 0.5);
                NonSpongePointsinCut[ClusterIndex] = (int)Math.Floor(AccumulateNonSpongePoints.TotalSum[ClusterIndex] + 0.5);
                ClusterIllegalPoints[ClusterIndex] = FindIllegalPoints[ClusterIndex].TotalInt;
                ClusterMaxWidths[ClusterIndex] = FindMaxExtension[ClusterIndex].TotalMax;
                TotalIllegalPoints += ClusterIllegalPoints[ClusterIndex];
            }
            TotalConfusedClusterPoints = ConfusedClusterPoints.Total;
            TotalConfusedSpongePoints1 = ConfusedSpongePoints1.Total;
            TotalConfusedSpongePoints2 = ConfusedSpongePoints2.Total;

        }   // End SetPointStatistics()

        public void SetNearbyClusters()
        {   
            Range[] ParallelClusterRange = RangePartitioner.Partition(NumberClusters, DAVectorUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ClusterThreadNo) =>
            {
                int beginindex = ParallelClusterRange[ClusterThreadNo].StartIndex;
                int indexlength = ParallelClusterRange[ClusterThreadNo].Length;

                double[] ClusterSigma = new double[Program.ParameterVectorDimension];
                int[] localindices = new int[NumberNearbyClusters];
                double[] localdistces = new double[NumberNearbyClusters];

                for (int ClusterIndex1 = beginindex; ClusterIndex1 < beginindex + indexlength; ClusterIndex1++)
                {
                    if (ClusterIndex1 == SpongeOption)
                        continue;
                    for (int nearbyloop = 0; nearbyloop < NumberNearbyClusters; nearbyloop++)
                    {
                        localindices[nearbyloop] = -1;
                        localdistces[nearbyloop] = 0.0;
                    }
                    int worstposition = -1;
                    for (int ClusterIndex2 = 0; ClusterIndex2 < NumberClusters; ClusterIndex2++)
                    {
                        if (ClusterIndex2 == SpongeOption)
                            continue;
                        if (ClusterIndex1 == ClusterIndex2)
                            continue;
                        Program.CalculateSigma(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1], ref ClusterSigma);
                        double tmp = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex1],
                            ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex2], ClusterSigma);
                        GlobalReductions.FindMinimumSet(tmp, ClusterIndex2, ref worstposition, localdistces, localindices, NumberNearbyClusters);
                    }

                    for (int countpos = 0; countpos < NumberNearbyClusters; countpos++)
                    {
                        int nearbyloop = -1;
                        for (int nearbyloop1 = 0; nearbyloop1 < NumberNearbyClusters; nearbyloop1++)
                        {
                            if (localindices[nearbyloop1] < 0)
                                continue;
                            if ((nearbyloop < 0) || (localdistces[nearbyloop1] <= localdistces[nearbyloop]))
                                nearbyloop = nearbyloop1;
                        }
                        NearbyClusterIndices[ClusterIndex1][countpos] = localindices[nearbyloop];
                        NearbyClusterDistances[ClusterIndex1][countpos] = localdistces[nearbyloop];
                        localindices[nearbyloop] = -1;
                        localdistces[nearbyloop] = 0.0;
                    }
                }
            });

        }   // End SetNearbyClusters()

        public void OutputStatus()
        {
            if (DAVectorUtility.MPI_Rank != 0)
                return;
            string nextline = "\n---------------- Cluster Status with Cut " + NearbyCut.ToString("F4") + " Nearby Cluster Limit " + NumberNearbyClusters.ToString() + " WRITTEN to status file " + " Sponge " + SpongeOption.ToString()
                + "\nIllegal Points outside cut of center MOVED to Sponge " + TotalIllegalPoints.ToString() + " Confused Sponge Points MOVED to Real Cluster " + TotalConfusedSpongePoints1.ToString("F0") + " After Correction " + TotalConfusedSpongePoints2.ToString("F0")
                + " Confused Cluster Points near >1 cluster " + TotalConfusedClusterPoints.ToString("F0") + " Confused Cluster Fraction (Clusters within Cut)" + AverageClustersinCut.ToString("F2");
            DAVectorUtility.SALSAPrint(0, nextline);

            nextline = "Cluster Occupation Count Histogram\n";
            for (int histloop = 0; histloop < HistogramOccupationMax; histloop++)
                nextline += HistogramOccupationValues[histloop].ToString() + ", ";
            nextline += "\n\nInter Cluster Distance Histogram\n";
            for (int histloop = 0; histloop < HistogramDistancesMax; histloop++)
                nextline += HistogramDistancesValues[histloop].ToString() + ", ";
            DAVectorUtility.SALSAPrint(0, nextline);

            //  Output Global Counts versus Temperature
            ClusterQuality.CaculateTemperatureClusterCountPlot();

            string directory = Path.GetDirectoryName(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string file = Path.GetFileNameWithoutExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile) + "Status1" + "-M" + Program.maxNcentperNode.ToString()
                + "-C" + ParallelClustering.RunningSolution.Ncent_Global.ToString() + Path.GetExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string StatusFileName = Path.Combine(directory, file);

            using (FileStream stream = new FileStream(StatusFileName, FileMode.Create, FileAccess.Write, FileShare.None))
            {
                using (StreamWriter writer = new StreamWriter(stream))
                {
                    nextline = "\nLabel\tY-0\tY-1\tCount\tNearby Clusters\tNearby Sponge Points\tIllegal Points\tMax Extension\tNearby Cluster Points";
                    writer.WriteLine(nextline);
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                    {
                        if (ClusterIndex == SpongeOption)
                            continue;
                        nextline = DAVectorUtility.PrintFixedInteger(ClusterIndex, 5);
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        {
                            double tmp = ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex][VectorIndex];
                            if ((Program.SigmaMethod >= 2) && (VectorIndex == 0))
                                tmp = Math.Log(tmp);
                            tmp = tmp / Program.SigmaVectorParameters_i_[VectorIndex];
                            nextline += "\t" + DAVectorUtility.PadString(tmp.ToString("F4"), 10);
                        }
                        nextline +=
                             "\t" + DAVectorUtility.PrintFixedInteger(ClusteringSolution.TotalClusterSummary.OccupationCount[ClusterIndex], 4)
                            + "\t" + DAVectorUtility.PrintFixedInteger(NumberClustersinCut[ClusterIndex], 4)
                            + "\t" + DAVectorUtility.PrintFixedInteger(SpongePointsinCut[ClusterIndex], 4)
                            + "\t" + DAVectorUtility.PrintFixedInteger(ClusterIllegalPoints[ClusterIndex], 4)
                            + "\t" + DAVectorUtility.PadString(ClusterMaxWidths[ClusterIndex].ToString("F3"), 8)
                            + "\t" + DAVectorUtility.PrintFixedInteger(NonSpongePointsinCut[ClusterIndex], 4);

                        string inmessage = "";
                        for (int nearbyloop = 0; nearbyloop < NumberNearbyClusters; nearbyloop++)
                        {
                            int ClusterNeighbor = NearbyClusterIndices[ClusterIndex][nearbyloop];
                            if ( ClusterIndex != ClusterNeighbor )
                                inmessage += "\t" + ClusterNeighbor.ToString() + "(" + NearbyClusterDistances[ClusterIndex][nearbyloop].ToString("F2") + ") ";
                        }
                        if (inmessage.Length > 0)
                            nextline += inmessage;
                        writer.WriteLine(nextline);
                    }
                }
            }

        }   // End OutputStatus()

        public static void CaculateTemperatureClusterCountPlot()
        {
            //  Output Global Counts versus Temperature
            int TemperatureCount = DAVectorUtility.TemperatureValues.Count;
            int step = TemperatureCount / 500;
            step = Math.Max(1, step);
            int loopnumber = 1 + (TemperatureCount - 1) / step;
            int Tindex = 0;
            string OutputLine = "\nCluster Count v Temperature\n";
            for (int loop = 0; loop < loopnumber; loop++)
            {
                if (Tindex >= TemperatureCount)
                    Tindex = TemperatureCount - 1;
                double value = (double)DAVectorUtility.TemperatureValues[Tindex];
                OutputLine += value.ToString("E4") + ", ";
                Tindex += step;
            }
            OutputLine += "\n";
            Tindex = 0;
            for (int loop = 0; loop < loopnumber; loop++)
            {
                if (Tindex >= TemperatureCount)
                    Tindex = TemperatureCount - 1;
                int ivalue = (int)DAVectorUtility.ClusterCountValues[Tindex];
                OutputLine += ivalue.ToString() + ", ";
                Tindex += step;
            }
            DAVectorUtility.SALSAPrint(0, OutputLine);
        }

    }   // End Class ClusterQuality

}