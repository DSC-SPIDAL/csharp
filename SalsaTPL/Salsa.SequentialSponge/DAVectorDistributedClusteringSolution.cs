using System;
using System.Threading;
using System.Threading.Tasks;
using MPI;
using Salsa.Core;


using Salsa.DAVectorSponge;
using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    // Hybrid Point and Cluster Parallelism
    // We exploit a natural one dimensional decomposition in target problem but most ideas in method are general
    //  One dimensional decomposition simplifies MPI communication and determination of clusters of interest
    //
    //  On input assume data (in Y0) is sorted with Processor 0 holding lowest values in sorted decomposed quantity
    //  Generate cuts in Y0 based on equal number of points per hread
    //  Clusters are assigned based on these point generated cuts. One processor makes these decisions and so ambiguity at cuts is not important
    //  A processor can control a cluster stored just outside it quite correctly
    //
    //  Clusters are labelled by Created Index which is 1 + CreatingProcessor + Creation Number in Processor * MULTIPLIER
    //  There is a look up table in every node that maps CreatedIndex to 
    //  Global Cluster Numbers
    //  Clusters controlled by THIS node and 
    //  Clusters controlled by other nodes but transported to this node
    //  Nothing for rest of clusters
    //  Created Indices are Unique and information set labelled by iteration count to identify stale data
    //
    //  Set up all key initial parameters and arrays
    //  Start in Global mode
    //  When clusters reach some limit, switch to distributed mode
    //  Do Major Synchronization involving
    //      Transmission of Cluster values and Changes
    //      Iteration of Point--Cluster assignments
    //      Initialize distributed reductions
    //
    //  Iterate Little loops updating Y T using distributed reductions with "minor synchronization"
    //
    //  Examine quality of clusters and splitting (involves distributed computation of correlations)
    //  Do Major Synchronization if clusters deleted or split or just have changed a lot
 
    //  Data structures used to store ALL (local and remote)Clusters in Node in Reductions
    public class NodeAccumulationMetaData
    {
        public int NumberofPointsperNode = 0;   // Count of Positions in Node Accumulation
        public int[] NumberofPointsperThread;   // Count of Positions in Thread Accumulations
        public int[] NodeAccumulationCreatedIndices;    // Current list of CreatedIndices for Node Accumulations
        public int[] NodeAccumulationClusterStatus;    //  -2 -1 Deleted 0 Sponge 1 Global 2 Local Distribututed 3 From Another Node
        public int[] NodeAccumulationClusterHosts;       // Original Hosts in form H + PACKINGMULTIPLIER (H1 + PACKINGMULTIPLIER * H2 )
        public int[][] AccumulationNodetoThreadClusterAssociations;    // Associated Thread Locations for Node Accumulations

        public NodeAccumulationMetaData(int MaxAccumulations)
        {
            NumberofPointsperThread = new int[DAVectorUtility.ThreadCount];
            NodeAccumulationCreatedIndices = new int[MaxAccumulations];
            NodeAccumulationClusterStatus = new int[MaxAccumulations];
            NodeAccumulationClusterHosts = new int[MaxAccumulations];
            AccumulationNodetoThreadClusterAssociations = new int[MaxAccumulations][];
            for (int AccumulationIndex = 0; AccumulationIndex < MaxAccumulations; AccumulationIndex++)
                AccumulationNodetoThreadClusterAssociations[AccumulationIndex] = new int[DAVectorUtility.ThreadCount];
        }

    }   // End NodeAccumulationMetaData

    //  Used to temporarily store data on Distributed Clusters hosted on this node for transportation to other nodes
    public class TemporaryLocalClusterInfoforTransport
    {
        public int SizeOfTransportedArray;   // Current Size of Transported Array
        public double[][] TotalTransportedY_t_i;   // Center Y_t_i of Transported Cluster with last index P[k]so dimension 
        public int[] TotalTransportedCreatedIndex;   // Created Index of Transported Cluster
        public int[] TotalTransportedOriginalHost;   // Host (i.e. node or MPI process) where cluster home in form H + MULTIPLIER (H1 + MULTIPLIER * H2 )
        public int[][] TotalTransportedStatus;   // Status of Transported Cluster [0] -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here (Not possible) 3 Distributed Stored Elsewhere (default)
                                                //  Status[1] 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED

        public TemporaryLocalClusterInfoforTransport(int MaxStorage)
        {
            SizeOfTransportedArray = 0;
            TotalTransportedOriginalHost = new int[MaxStorage];
            TotalTransportedCreatedIndex = new int[MaxStorage];
            TotalTransportedStatus = new int[MaxStorage][];
            TotalTransportedY_t_i = new double[MaxStorage][];
            for (int TransportIndex = 0; TransportIndex < MaxStorage; TransportIndex++)
            {
                TotalTransportedY_t_i[TransportIndex] = new double[Program.ParameterVectorDimension + 1];
                TotalTransportedStatus[TransportIndex] = new int[2];
            }
        }

    }   // End LocalDistributedClusterInformation

    //  Used to store data on Clusters transported from other nodes
    public class TransportedClusterStoredInformation
    {
        public int SizeOfTransportedArray;   // Current Size of Transported Array
        public double[][] TotalTransportedY_t_i;   // Center Y_t_i of Transported Cluster (with one extra index to load P_k)
        public double[] TotalTransportedP_t;   // Multiplier P_t of Transported Cluster
        public double[][] TotalTransportedSigma_t_i; // Sigma of Transported Cluster (calculated NOT transported)
        public int[] TotalTransportedCreatedIndex;   // Created Index of Transported Cluster
        public int[] TotalTransportedOriginalHost;   // Host (i.e. node or MPI process) where cluster home in form H + MULTIPLIER (H1 + MULTIPLIER * H2 )
        public int[][] TotalTransportedStatus;   // Status of Transported Cluster [0] -2 Moved -1 deleted 0 Sponge 1 Global 2 Distributed Stored here (Not possible) 3 Distributed Stored Elsewhere (default)
                                                 //  Status[1] 1 + CreatedIndex of Parent or -1 -CreatedIndex of child if Split UNTIL SPLIT PROPAGATED
        public int[] TransportedNodeAccPosition;    // Position in Node Accumulation of this remote cluster; -1 NOT used
        public int[][] TransportedThreadAccPosition;    // Position in Thread Accumulation of this remote cluster; -1 NOT used

        public TransportedClusterStoredInformation(int MaxStorage)
        {
            SizeOfTransportedArray = 0;
            TotalTransportedOriginalHost = new int[MaxStorage];
            TotalTransportedCreatedIndex = new int[MaxStorage];
            TotalTransportedStatus = new int[MaxStorage][];
            TotalTransportedY_t_i = new double[MaxStorage][];
            TotalTransportedP_t = new double[MaxStorage];
            TotalTransportedSigma_t_i = new double[MaxStorage][];
            TransportedNodeAccPosition = new int[MaxStorage];
            TransportedThreadAccPosition = new int[MaxStorage][];

            for (int TransportIndex = 0; TransportIndex < MaxStorage; TransportIndex++)
            {
                TotalTransportedY_t_i[TransportIndex] = new double[Program.ParameterVectorDimension + 1];
                TotalTransportedSigma_t_i[TransportIndex] = new double[Program.ParameterVectorDimension];
                TotalTransportedStatus[TransportIndex] = new int[2];
                TransportedThreadAccPosition[TransportIndex] = new int[DAVectorUtility.ThreadCount];
            }

        }
    }   // End TransportedClusterStoredInformation

    public class DistributedClusteringSolution
    {

        public static int MaxMPITransportBuffer;  // Maximum size of buffer for input and output transport
        public static int MaxTransportedClusterStorage;  // Maximum size of Storage for input and output transport
        public static int MaxNumberAccumulationsperNode;  // Maximum size of buffer for Node Accumulations
        public static int MaxDoubleComponents;    // Maximum Number of Double Components Allowed
        public static int MaxIntegerComponents;    // Maximum Number of Integer Components Allowed

        public static double[] NodeCutUpperLimit;   // Upper limit of one dimemsional cut for each node
        public static double[] ThreadCutUpperLimit;    // Upper limit of one dimensional cut for each thread in this node; thread cuts in other nodes are NOT known

        public static TransportedClusterStoredInformation StorageforTransportedClusters;    // Structure Storing Distributed Clusters hosted elsewhere
        public static TemporaryLocalClusterInfoforTransport TemporaryLocalClustersforSynchronization;  // Structure for storing Local Cluster Data for Synchronization

        public static Range[] ParallelNodeAccumulationRanges;   // Ranges for parallel loops over Node Accumulations
        public static NodeAccumulationMetaData NodeAccMetaData; // Hold all the metadata for Node Accumulatio

        public static MPI.CompletedStatus MPISecStatus; // MPI Send Receive Status

        public static string HostLinkageMessage;    // Holds Information on Linkage

        public DistributedClusteringSolution(int MaxBuf, int MaxStorage, int MaxAcc, int _MaxDoubleComponents, int _MaxIntegerComponents)
        {   // Initialize Distributed Clusters

            MaxTransportedClusterStorage = MaxStorage;
            MaxNumberAccumulationsperNode = MaxAcc;
            MaxDoubleComponents = _MaxDoubleComponents;
            MaxIntegerComponents = _MaxIntegerComponents;
            MaxMPITransportBuffer = MaxBuf;

            //  This is place on node where information on clusters external to this node are stored
            StorageforTransportedClusters = new TransportedClusterStoredInformation(MaxStorage);

            //  This is temporary storage for local cluster data so it can be pipelined to other nodes
            TemporaryLocalClustersforSynchronization = new TemporaryLocalClusterInfoforTransport(MaxStorage);

            //  Accumulation Arrays to support distributed reductions
            NodeAccMetaData = new NodeAccumulationMetaData(MaxNumberAccumulationsperNode);

            //  Define Responsibity of Nodes and Threads
            SetOneDimensionalUpperLimits();

            //  Array for RECEIVING MPI Information
            DistributedSynchronization.TransportComponent = new MPITransportComponentPacket(MaxMPITransportBuffer, MaxDoubleComponents, MaxIntegerComponents);

        }   // End DistributedClusteringSolution

        //  Note that in ambiguous cases (Cluster value = Upper limit, the current host node makes the decision and broadcasts that
        //  It won't matter which decision it makes 
        public static void SetOneDimensionalUpperLimits()
        {
            double[] NodeCutLowerLimit = new double[DAVectorUtility.MPI_Size];
            NodeCutUpperLimit = new double[DAVectorUtility.MPI_Size];
            ThreadCutUpperLimit = new double[DAVectorUtility.ThreadCount];

            double Lower = Program.PointPosition[0][0];
            double Upper = Program.PointPosition[DAVectorUtility.PointCount_Process - 1][0];

            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIGATHERTiming);
            NodeCutLowerLimit = DAVectorUtility.MPI_communicator.Allgather<double>(Lower);
            NodeCutUpperLimit = DAVectorUtility.MPI_communicator.Allgather<double>(Upper);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIGATHERTiming);

            string message = "\nPosition Cuts";
            for (int rank = 0; rank < DAVectorUtility.MPI_Size - 1; rank++)
            {
                NodeCutUpperLimit[rank] = 0.5 * (NodeCutUpperLimit[rank] + NodeCutLowerLimit[rank + 1]);
                message += " " + NodeCutUpperLimit[rank].ToString("E4");
            }
            DAVectorUtility.SALSAPrint(1, message);

            for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
            {
                int indextop = DAVectorUtility.PointsperThread[ThreadNo] + DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process - 1;
                ThreadCutUpperLimit[ThreadNo] = Program.PointPosition[indextop][0];
                if (ThreadNo != (DAVectorUtility.ThreadCount - 1))
                    ThreadCutUpperLimit[ThreadNo] = 0.5 * (ThreadCutUpperLimit[ThreadNo] + Program.PointPosition[indextop + 1][0]);
            }

        }   // End SetOneDimensionalUpperLimits()

        //  Set Cluster Host associations needed before any Major Synchronization
        public static void ClusterHostlinkage(bool useClustertype1)
        {
            int TestType = 2;
            if (useClustertype1)
                TestType = 1;

            string MoveMessage = "";
            string UpMessage = "";
            string DownMessage = "";
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                int H = DAVectorUtility.MPI_Rank;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 0)
                {
                    ClusteringSolution.LocalHost[ActiveClusterIndex] = H | (H | (H << ClusteringSolution.PACKINGSHIFT)) << ClusteringSolution.PACKINGSHIFT;
                    ClusteringSolution.UniversalMapping[ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
                    continue;
                }
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] != TestType)
                    continue;

                // New Children MUST be identical to Parent
                int Parent = ParallelClustering.RunningSolution.LocalSplitCreatedIndex[RealClusterIndex];
                if (Parent > 0)
                {
                    ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[Parent - 1].PackedHost;
                    ClusteringSolution.UniversalMapping[ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost = ClusteringSolution.UniversalMapping[Parent - 1].PackedHost;
                    continue;
                }
                double TestPosition = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][0];

                H = 0;
                while (H < DAVectorUtility.MPI_Size - 1)
                {
                    if (TestPosition < NodeCutUpperLimit[H])
                        break;
                    H++;
                }

                int H1 = H;
                while (H1 >= 1)
                {
                    double tmp = TestPosition - NodeCutUpperLimit[H1 - 1];
                    if (tmp > 0.0)
                    {
                        double scaled1Ddistance = tmp * tmp / ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex][0];
                        if (scaled1Ddistance >= (2.0 * Program.ExpArgumentCut3 * ParallelClustering.RunningSolution.Temperature))
                            break;
                    }
                    H1--;
                }

                int H2 = H;
                while (H2 < DAVectorUtility.MPI_Size - 1)
                {
                    double tmp = TestPosition - NodeCutUpperLimit[H2];
                    if (tmp < 0.0)
                    {
                        double scaled1Ddistance = tmp * tmp / ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex][0];
                        if (scaled1Ddistance >= (2.0 * Program.ExpArgumentCut3 * ParallelClustering.RunningSolution.Temperature))
                            break;
                    }
                    H2++;
                }
                int CreatedIndex = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                if (useClustertype1)
                {
                    if (H == DAVectorUtility.MPI_Rank)
                    {
                        if (H1 != DAVectorUtility.MPI_Rank)
                            DownMessage += CreatedIndex.ToString() + " ";
                        if (H2 != DAVectorUtility.MPI_Rank)
                            UpMessage += CreatedIndex.ToString() + " ";
                    }
                }
                else
                {
                    if (H != DAVectorUtility.MPI_Rank)
                        MoveMessage += CreatedIndex.ToString() + " ";
                    if (H1 != DAVectorUtility.MPI_Rank)
                        DownMessage += CreatedIndex.ToString() + " ";
                    if (H2 != DAVectorUtility.MPI_Rank)
                        UpMessage += CreatedIndex.ToString() + " ";
                }

                H1 = H1 | (H2 << ClusteringSolution.PACKINGSHIFT);
                ClusteringSolution.LocalHost[ActiveClusterIndex] = H | (H1 << ClusteringSolution.PACKINGSHIFT);
                ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
            }

            HostLinkageMessage = "";
            if(MoveMessage.Length !=0)
                HostLinkageMessage += "Move " + MoveMessage;
            if(DownMessage.Length != 0)
                HostLinkageMessage += "Down " + DownMessage;
            if(UpMessage.Length != 0)
                HostLinkageMessage += "Up " + UpMessage;
            if (!useClustertype1)
                return;

            //  If Processing Status 1 clusters, then take results from Processor 0
            //  If Processing Status 2 clusters, then each Processor decides on utility and placement of its clusters
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            DAVectorUtility.MPI_communicator.Broadcast<int>(ref ClusteringSolution.LocalHost, 0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

        }   // End ClusterHostlinkage()

        //  Return Array pointers associated with any cluster based on CreatedIndex hashing
        public static int IndicesperCluster(int CreatedIndex, int ThreadNo, ref int LocalClusterIndex, ref int TransportedClusterIndex_instorage, ref int NodeAccumulationPosition, ref int ThreadAccumulationPosition)
        {
            if (ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet < ClusteringSolution.CurrentIteration)
                return -1;
            int position = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
            if (position == 0)
                return -2;
            if (position > 0)
            {
                TransportedClusterIndex_instorage = -1;
                LocalClusterIndex = position - 1;
                int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[LocalClusterIndex];
                NodeAccumulationPosition = ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex];
                if (NodeAccumulationPosition == -1)
                    return -3;
                ThreadAccumulationPosition = -1;
                if (ThreadNo >= 0)
                    ThreadAccumulationPosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo];
                return 0;
            }
            LocalClusterIndex = -1;
            TransportedClusterIndex_instorage = -position - 1;
            NodeAccumulationPosition = StorageforTransportedClusters.TransportedNodeAccPosition[TransportedClusterIndex_instorage];
            if (NodeAccumulationPosition == -1)
                return -4;
            ThreadAccumulationPosition = -1;
            if (ThreadNo >= 0)
                ThreadAccumulationPosition = StorageforTransportedClusters.TransportedThreadAccPosition[TransportedClusterIndex_instorage][ThreadNo];
            return 1;

        }   // End IndicesperCluster

        //  Only called in distributed mode
        //  Perform Minor Synchronization when cluster positions Updated
        public static void MinorSynchronizationTransportDistributedClusterCenters()
        {
            DAVectorUtility.StartSubTimer(7);
            int NumberofLocalDistributedClusters = 0;
            ++Program.NumberMinorSynchs;
            if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
            {
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    DAVectorUtility.SALSAPrint(1, "Minor Synchronization " + Program.NumberMinorSynchs.ToString() + " Iteration " + ClusteringSolution.CurrentIteration.ToString()
                        + " Temperature " + ParallelClustering.RunningSolution.Temperature.ToString("F4") + " Cluster Count " + ParallelClustering.RunningSolution.Ncent_Global.ToString());
            }
            int H = DAVectorUtility.MPI_Rank;

            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                    continue;
                int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
                if ((H1 == H) && (H2 == H))
                {
                    if ((PackedHost & ClusteringSolution.PACKINGMASK) != H)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Inconsistent Host Range " + (PackedHost & ClusteringSolution.PACKINGMASK).ToString() + " in Node " + H.ToString());
                        throw (e);
                    }
                    continue;
                }
                TemporaryLocalClustersforSynchronization.TotalTransportedOriginalHost[NumberofLocalDistributedClusters] = PackedHost;
                TemporaryLocalClustersforSynchronization.TotalTransportedCreatedIndex[NumberofLocalDistributedClusters] = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i[NumberofLocalDistributedClusters][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i[NumberofLocalDistributedClusters][Program.ParameterVectorDimension] = ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                ++NumberofLocalDistributedClusters;
            }   // End Loop over LocalClusterIndex

            TemporaryLocalClustersforSynchronization.SizeOfTransportedArray = NumberofLocalDistributedClusters;
            Program.ActualMaxTransportedClusterStorage = Math.Max(Program.ActualMaxTransportedClusterStorage, NumberofLocalDistributedClusters);
            DistributedSynchronization.TransportviaPipeline DoSynchronization = new DistributedSynchronization.TransportviaPipeline(0, true, 3, 0, NumberofLocalDistributedClusters, 2);

            int FinalClusterCount = 0;  // Dummy
            DoSynchronization.PipelineDistributedBroadcast(TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i, StorageforTransportedClusters.TotalTransportedY_t_i, null, null,
                TemporaryLocalClustersforSynchronization.TotalTransportedCreatedIndex,
                StorageforTransportedClusters.TotalTransportedCreatedIndex,
                TemporaryLocalClustersforSynchronization.TotalTransportedOriginalHost, StorageforTransportedClusters.TotalTransportedOriginalHost, ref FinalClusterCount);

            // Load P_t and Reset Sigmas if necessary
            for (int Clustersfromafar = 0; Clustersfromafar < StorageforTransportedClusters.SizeOfTransportedArray; Clustersfromafar++)
            {
                StorageforTransportedClusters.TotalTransportedP_t[Clustersfromafar] = StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar][Program.ParameterVectorDimension];
                if ( (Program.SigmaMethod > 1) && (StorageforTransportedClusters.TotalTransportedStatus[Clustersfromafar][0] == 3) )
                    Program.CalculateSigma(StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar],
                    ref StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar]);
            }
            DAVectorUtility.StopSubTimer(7);

        }   // End MinorSynchronizationTransportDistributedClusterCenters()

        //  Perform Major Synchronization Data Transport needed when Cluster structure changes
        //  The CurrentIteration value is incremented at such events and information about splits and deletes are propagated
        public static void MajorSynchronizationTransportDistributedClusterCenters()
        {

            int NumberofLocalDistributedClusters = 0;
            int H = DAVectorUtility.MPI_Rank;
            // string message = "";

            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                    continue;
                int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
                if ((H1 == H) && (H2 == H))
                {
                    if ((PackedHost & ClusteringSolution.PACKINGMASK) != H)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Inconsistent Host Range " + (PackedHost & ClusteringSolution.PACKINGMASK).ToString() + " in Node " + H.ToString());
                        throw (e);
                    }
                    continue;
                }
                // message += RealClusterIndex.ToString() + "(" + ParallelClustering.RunningDAVectorSolt.LocalCreatedIndex[RealClusterIndex].ToString() + ") ";
                TemporaryLocalClustersforSynchronization.TotalTransportedOriginalHost[NumberofLocalDistributedClusters] = PackedHost;
                TemporaryLocalClustersforSynchronization.TotalTransportedCreatedIndex[NumberofLocalDistributedClusters] = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                TemporaryLocalClustersforSynchronization.TotalTransportedStatus[NumberofLocalDistributedClusters][0] = 3;
                TemporaryLocalClustersforSynchronization.TotalTransportedStatus[NumberofLocalDistributedClusters][1] = ParallelClustering.RunningSolution.LocalSplitCreatedIndex[RealClusterIndex];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i[NumberofLocalDistributedClusters][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i[NumberofLocalDistributedClusters][Program.ParameterVectorDimension] = ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                ++NumberofLocalDistributedClusters;
            }   // End Loop over LocalClusterIndex

            TemporaryLocalClustersforSynchronization.SizeOfTransportedArray = NumberofLocalDistributedClusters;
            Program.ActualMaxTransportedClusterStorage = Math.Max(Program.ActualMaxTransportedClusterStorage, NumberofLocalDistributedClusters);
            DistributedSynchronization.TransportviaPipeline DoSynchronization = new DistributedSynchronization.TransportviaPipeline(2, true, 3, 1, NumberofLocalDistributedClusters, 0);

            int FinalClusterCount = 0;
            DoSynchronization.PipelineDistributedBroadcast(TemporaryLocalClustersforSynchronization.TotalTransportedY_t_i, StorageforTransportedClusters.TotalTransportedY_t_i,
                TemporaryLocalClustersforSynchronization.TotalTransportedStatus, StorageforTransportedClusters.TotalTransportedStatus,
                TemporaryLocalClustersforSynchronization.TotalTransportedCreatedIndex,
                StorageforTransportedClusters.TotalTransportedCreatedIndex,
                TemporaryLocalClustersforSynchronization.TotalTransportedOriginalHost, StorageforTransportedClusters.TotalTransportedOriginalHost, ref FinalClusterCount);
            StorageforTransportedClusters.SizeOfTransportedArray = FinalClusterCount;
            Program.ActualMaxTransportedClusterStorage = Math.Max(FinalClusterCount + 1, Program.ActualMaxTransportedClusterStorage);

            // Set P_t and sigmas which were not transported
            // Add Created Index if new
            for (int Clustersfromafar = 0; Clustersfromafar < StorageforTransportedClusters.SizeOfTransportedArray; Clustersfromafar++)
            {
                int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[Clustersfromafar];
                if (ClusteringSolution.UniversalMapping[CreatedIndex] == null)
                    ClusteringSolution.UniversalMapping[CreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, -1 - Clustersfromafar);
                ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[Clustersfromafar];
                StorageforTransportedClusters.TotalTransportedP_t[Clustersfromafar] = StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar][Program.ParameterVectorDimension];
                if (StorageforTransportedClusters.TotalTransportedStatus[Clustersfromafar][0] == 3)
                    Program.CalculateSigma(StorageforTransportedClusters.TotalTransportedY_t_i[Clustersfromafar],
                        ref StorageforTransportedClusters.TotalTransportedSigma_t_i[Clustersfromafar]);
            }

        }   // End MajorSynchronizationTransportDistributedClusterCenters()

        // public static int PrintOuts = 0;
        public static void SetClustersforaPoint(ClusteringSolution Solution)
        {
            DAVectorUtility.StartSubTimer(13);
            GlobalReductions.FindVectorDoubleSum FindDiagnosticSums_Points = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, 8);


            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int[] ProposedClusters = new int[Program.maxNcentperPoint];
                double[] ProposedDistances = new double[Program.maxNcentperPoint];

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {   // Loop over points

                    int NumberofClustersTobeUsed = 0;
                    // double SmallDistce = -1.0;
                    FindDiagnosticSums_Points.addapoint(ThreadNo, 1.0, 0);
                    // int AssignmentStatus = DATriangleInequality.SetAssociatedCenters(alpha, out NumberofClustersTobeUsed, ProposedClusters, out SmallDistce, ProposedDistances );
                    int AssignmentStatus = DATriangleInequality.SetAssociatedCenters(alpha, out NumberofClustersTobeUsed, ProposedClusters);
                    FindDiagnosticSums_Points.addapoint(ThreadNo, 1.0, 2 + AssignmentStatus);
                    if (AssignmentStatus < 0)
                    {
                            Exception e = DAVectorUtility.SALSAError("Error in Triangle Inequality at Point " +
                                (alpha + DAVectorUtility.PointStart_Process).ToString() + " Old Number " + Solution.NumClusters_alpha_[alpha].ToString());
                            throw e;
                    }

                    if (NumberofClustersTobeUsed <= 0)
                    {
                        Exception e = DAVectorUtility.SALSAError("Error due to zero New Number of Clusters Point " +
                            (alpha + DAVectorUtility.PointStart_Process).ToString() + " Old Number " + Solution.NumClusters_alpha_[alpha].ToString());
                        throw e;
                    }

                    //  Now make list for storing back
                    double[] Mvalues = new double[NumberofClustersTobeUsed];
                    double Msum = 0.0;

                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                    {
                        Mvalues[IndirectClusterIndex] = 0.0;
                      
                        if( Solution.NumClusters_alpha_[alpha] == 0 )
                            continue;

                        int ActiveClusterIndex = ProposedClusters[IndirectClusterIndex];
                        int OldIndirectClusterIndex = Solution.MapClusterToIndirect(alpha, ActiveClusterIndex);
                        if (OldIndirectClusterIndex < 0)
                            continue;
                        Mvalues[IndirectClusterIndex] = Solution.M_alpha_kpointer_[alpha][OldIndirectClusterIndex];
                        Msum += Mvalues[IndirectClusterIndex];
                    }
                    if (Msum < 0.2)
                    {
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                            Mvalues[IndirectClusterIndex] = 1.0 / (double)NumberofClustersTobeUsed;
                        Msum = 1.0;
                        FindDiagnosticSums_Points.addapoint(ThreadNo, 1.0, 6);
                    }

                    //  Reset List
                    int NumberChange = Solution.NumClusters_alpha_[alpha] - NumberofClustersTobeUsed;
                    /*
                    if ( (Solution.Ncent_Global == 140) && (Solution.Temperature < 0.1) && (PrintOuts < 100) && ( NumberChange != 0) && (DAVectorUtility.MPI_Rank == 0) )
                    {
                        int bestindex = ClusteringSolution.RealClusterIndices[ProposedClusters[0]];
                        double deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[bestindex],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[bestindex]);

                        string message = " alpha " + alpha.ToString() + " size " + NumberofClustersTobeUsed.ToString() + " # Chg " + NumberChange.ToString()
                            + " Small " + SmallDistce.ToString("E4") + " Best " +  bestindex.ToString() + " " + deltaY.ToString("E4") + " new ";
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                        {
                            int ActiveClusterIndex = ProposedClusters[IndirectClusterIndex];
                            int OldIndirectClusterIndex = Solution.MapClusterToIndirect(alpha, ActiveClusterIndex);
                            if (OldIndirectClusterIndex >= 0)
                                continue;
                            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                            deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealClusterIndex],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[RealClusterIndex]);
                            message += RealClusterIndex.ToString() + " " + deltaY.ToString("E4") + " " + ProposedDistances[IndirectClusterIndex].ToString("E4") + " ";
                        }
                        message += " old ";

                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < Solution.NumClusters_alpha_[alpha]; IndirectClusterIndex++)
                        {
                            int RealClusterIndex1 = -1;
                            int ActiveClusterIndex1 = -1;
                            int RemoteIndex = 0;
                            VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex1, ref ActiveClusterIndex1, ref RemoteIndex);
                            for (int IndirectClusterIndex1 = 0; IndirectClusterIndex1 < NumberofClustersTobeUsed; IndirectClusterIndex1++)
                            {
                                if (ActiveClusterIndex1 == ProposedClusters[IndirectClusterIndex1])
                                {
                                    ActiveClusterIndex1 = -1;
                                    break;
                                }
                            }
                            if (ActiveClusterIndex1 < 0)
                                continue;
                            deltaY = DAVectorParallelism.getNOTSquaredScaledDistancebetweenVectors(DATriangleInequality.CenterY_k_i_Current[RealClusterIndex1],
                                DATriangleInequality.PointPosition[alpha], DATriangleInequality.CenterSigma_k_i_Current[RealClusterIndex1]);
                            message += RealClusterIndex1.ToString() + " " + deltaY.ToString("E4") + " " ;
                        }
                        DAVectorUtility.SALSAPrint(0, message);
                        ++PrintOuts;
                    }
                    */

                    FindDiagnosticSums_Points.addapoint(ThreadNo, (double) NumberChange, 7);
                    Solution.NumClusters_alpha_[alpha] = NumberofClustersTobeUsed;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                    {
                        Solution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = Mvalues[IndirectClusterIndex] / Msum;
                        Solution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = ClusteringSolution.MapActivetoCreatedIndex(ProposedClusters[IndirectClusterIndex], Solution);
                    }

                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NumberofClustersTobeUsed; IndirectClusterIndex++)
                        ParallelClustering.RunningSolution.LegalCluster(alpha, IndirectClusterIndex);
                }
            }); // End loop over Point dependent quantities

            Solution.DiffMalpha_k_Set = -1;

            FindDiagnosticSums_Points.sumoverthreadsandmpi();
            for (int DiagnosticLoop = 0; DiagnosticLoop < 8; DiagnosticLoop++)
                Program.ClustersperCenterDiagnostics[DiagnosticLoop] += FindDiagnosticSums_Points.TotalVectorSum[DiagnosticLoop];
            VectorAnnealIterate.NumberCountChanges = Program.ClustersperCenterDiagnostics[7] / Program.ClustersperCenterDiagnostics[0];
            DAVectorUtility.StopSubTimer(13);
            return;

        }   // End SetClustersforaPoint(ClusteringSolution Solution)

        // Control Major Synchronization for distributed and global case
        //  ongoing = true is normal case and works whether distributed or global mode
        //  ongoing = false switches INTO distributed mode
        public static void ManageMajorSynchronization(bool ongoing)
        {
            //  Set up Active Cluster Lists
            DAVectorUtility.StartSubTimer(6);
            ParallelClustering.RunningSolution.SetActiveClusters();
            VectorAnnealIterate.NumberCountChanges = 0.0;

            //  Case of non distributed operation
            if (ongoing && (!ParallelClustering.RunningSolution.DistributedExecutionMode))
            {
                ++ClusteringSolution.CurrentIteration; // Increment Iteration Number
                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    int CreatedIndex = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                    ClusteringSolution.UniversalMapping[CreatedIndex].Availability = 1 + RealClusterIndex;
                    ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
                    ClusteringSolution.LocalHost[LocalActiveClusterIndex] = 0;
                }
                if (Program.UseTriangleInequality_DA > 0)
                {
                    ParallelClustering.RunningSolution.SetClusterWidths();
                    DATriangleInequality.NextIteration();
                    DistributedClusteringSolution.SetClustersforaPoint(ParallelClustering.RunningSolution);
                }
                ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;

                ParallelClustering.RunningSolution.SetClusterSizes();
                ParallelClustering.RunningSolution.SetClusterWidths();
                ++Program.NumberMajorSynchs1;
                DAVectorUtility.StopSubTimer(6);
                return;
            }

            ++Program.NumberMajorSynchs2;
            Program.CountClusters2 += ParallelClustering.RunningSolution.Ncent_ThisNode;

            // Reset Host Information
            bool UseCluster1 = true;
            if (ongoing)
                UseCluster1 = false;
            ClusterHostlinkage(UseCluster1);

            // If update, Transport Cluster Information and populate local storage in StorageforTransportedClusters
            if (ongoing)
                MajorSynchronizationTransportDistributedClusterCenters();

            //  Set up Distributed Mode for first time setting StorageforTransportedClusters from local data
            else
                SetupDistributedMode();

            //  Process Clusters moved between nodes
            ProcessChangedClusters(ongoing);
            DAVectorUtility.StopSubTimer(6);

            //  Set up a new iteration
            SetupNewIteration();

            //  Initialize storage pointers
            DAVectorUtility.StartSubTimer(11);
            int LocalNumberClusters = ClusteringSolution.NumberLocalActiveClusters;
            int TotalNumberClusters = StorageforTransportedClusters.SizeOfTransportedArray + LocalNumberClusters;
            ClusteringSolution.NumberAvailableActiveClusters = TotalNumberClusters;
            if( TotalNumberClusters <= 0 )
            {
                Exception e = DAVectorUtility.SALSAError("No Clusters in this Process " + DAVectorUtility.MPI_Rank.ToString() + " Remote " + StorageforTransportedClusters.SizeOfTransportedArray.ToString()
                    + " Local " + LocalNumberClusters.ToString() );
                throw e;
            }


            Range[] TotalClusterRanges = RangePartitioner.Partition(TotalNumberClusters, DAVectorUtility.ThreadCount);

            //  Initialize Mapping of Local and Remote Clusters to Node and Thread Accumulations
            //  Initialize ClusteringSolution.UniversalMapping for this Node at this NEW Iteration number
            //  This assumes that any child cluster is created fully but values of Malpha are not set any where
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (AccumulationThreadNo) =>
            {
                int BeginClusterIndex = TotalClusterRanges[AccumulationThreadNo].StartIndex;
                int ClusterIndexLength = TotalClusterRanges[AccumulationThreadNo].Length;

                for (int LocalActiveClusterIndex = BeginClusterIndex; LocalActiveClusterIndex < BeginClusterIndex + ClusterIndexLength; LocalActiveClusterIndex++)
                {
                    int RemoteIndex = -1;
                    int RealClusterIndex = -1;
                    if (LocalActiveClusterIndex >= LocalNumberClusters)
                        RemoteIndex = LocalActiveClusterIndex - LocalNumberClusters;
                    else
                        RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    if (RemoteIndex == -1)
                    {
                        ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] = 0;   // INSIST all clusters hosted on this node have a node index
                        //  It may not have any thread indices if points for this cluster are outside this node
                        int CreatedIndex = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                        ClusteringSolution.UniversalMapping[CreatedIndex].Availability = 1 + RealClusterIndex;
                        ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
                    }
                    else
                    {
                        StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = -1;
                        int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
                        ClusteringSolution.UniversalMapping[CreatedIndex].Availability = -1 - RemoteIndex;
                        ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet = ClusteringSolution.CurrentIteration;
                    }
                    for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                    {
                        if (RemoteIndex == -1)
                            ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] = -1;
                        else
                            StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = -1;
                    }
                }

            }); // End Parallel Section over Total Number of Clusters

            //  Loop over points in this node
            //  Now find which clusters are actually used and remove clusters that have disappeared
            //  Process Pending Cluster Splits
            //  Set up Accumulation arrays by setting any used to 0 (initialized to -1 above)
         
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] Malpha_perpoint = new double[ClusteringSolution.NumberAvailableActiveClusters];
                int[] CreatedIndex_perpoint = new int[ClusteringSolution.NumberAvailableActiveClusters];

                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    int NewIndirectSize = 0;
                    double LostSize = 0.0;
                    bool changed = false;
                    for(int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int CreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        if ((CreatedIndex < 0 ) || (ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet != ClusteringSolution.CurrentIteration))
                        {   // A deleted Cluster or one no longer transported here -- perhaps because changed position
                            LostSize += ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            changed = true;
                            /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                            {
                                DAVectorUtility.SALSAFullPrint(0, " Lost " + IndirectClusterIndex.ToString() + " " + ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                    + " Mapped " + ClusteringSolution.UniversalMapping[CreatedIndex].Availability);
                            } */
                            continue;
                        }
                        /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                        {
                            DAVectorUtility.SALSAFullPrint(0, " OK " + IndirectClusterIndex.ToString() + " " + ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                + " Mapped " + ClusteringSolution.UniversalMapping[CreatedIndex].Availability);
                        } */
                        Malpha_perpoint[NewIndirectSize] = ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                        CreatedIndex_perpoint[NewIndirectSize] = CreatedIndex;
                        int AvailabilityIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                        int RemoteIndex = -1;
                        int SplitInfo = 0;
                        if (AvailabilityIndex > 0)
                        {
                            --AvailabilityIndex;
                            if (ParallelClustering.RunningSolution.LocalStatus[AvailabilityIndex] < 0)
                            {   // Actually deleted cluster -- remove it from list. Not certain if this could happen!
                                LostSize += ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                                changed = true;
                                continue;
                            }
                            int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[AvailabilityIndex];
                            SplitInfo = ParallelClustering.RunningSolution.LocalSplitCreatedIndex[AvailabilityIndex];
                            ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex] = 0;
                            ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo] = 0;
                        }
                        else
                        {
                            RemoteIndex = -AvailabilityIndex - 1;
                            SplitInfo = StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1];
                            StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = 0;
                            StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = 0;
                        }
                        /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                        {
                            DAVectorUtility.SALSAFullPrint(0, " OK Split " + IndirectClusterIndex.ToString() + " " + RemoteIndex.ToString() + " " + SplitInfo.ToString());
                        } */
                        if(SplitInfo > 0 )
                        {   // This indicates a child which is impossible
                            Exception e = DAVectorUtility.SALSAError(" Created Index " + CreatedIndex.ToString() + " has illegal split flag " + (SplitInfo-1).ToString());
                            throw (e);
                        }
                        if (SplitInfo < 0)
                        {   // This indicates a parent of a child
                            Malpha_perpoint[NewIndirectSize] *= 0.5;
                            ++NewIndirectSize;
                            Malpha_perpoint[NewIndirectSize] = Malpha_perpoint[NewIndirectSize - 1];
                            int ChildCreatedIndex = -SplitInfo - 1;
                            CreatedIndex_perpoint[NewIndirectSize] = ChildCreatedIndex;
                            changed = true;

                            if (ClusteringSolution.UniversalMapping[ChildCreatedIndex].IterationSet != ClusteringSolution.CurrentIteration)
                            {
                                Exception e = DAVectorUtility.SALSAError(" Child Created Index " + ChildCreatedIndex.ToString() + " is set up wrong -- parent " + CreatedIndex.ToString());
                                throw (e);
                            }
                            int ChildClusterIndex = ClusteringSolution.UniversalMapping[ChildCreatedIndex].Availability;
                            if (ChildClusterIndex > 0)
                            {
                                int ActiveChildClusterIndex = ClusteringSolution.ActiveClusterIndices[ChildClusterIndex - 1];
                                ClusteringSolution.LocalNodeAccPosition[ActiveChildClusterIndex] = 0;
                                ClusteringSolution.LocalThreadAccPosition[ActiveChildClusterIndex][ThreadNo] = 0;
                            }
                            else
                            {
                                StorageforTransportedClusters.TransportedNodeAccPosition[-ChildClusterIndex - 1] = 0;
                                StorageforTransportedClusters.TransportedThreadAccPosition[-ChildClusterIndex - 1][ThreadNo] = 0;
                            }
                        }
                        ++NewIndirectSize;
                    }

                    //  Reset Cluster Information for this point
                    if (changed)
                    {
                        if (NewIndirectSize > ClusteringSolution.TargetMinimumClustersperPoint)
                        {
                            ParallelClustering.RunningSolution.NumClusters_alpha_[alpha] = NewIndirectSize;
                            bool multiplication = true;
                            double fudge;
                            if (LostSize >= 0.99)
                            {
                                multiplication = false;
                                fudge = 0.0;
                            }
                            else
                                fudge = 1.0 / (1.0 - LostSize);
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                            {
                                ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex] = CreatedIndex_perpoint[IndirectClusterIndex];
                                if (multiplication)
                                    ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = Malpha_perpoint[IndirectClusterIndex] * fudge;
                                else
                                    ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = 1.0 / NewIndirectSize;
                            }
                        }   // End case where enough points

                        else
                        {   // Not enough Clusters for this Point
                            ParallelClustering.RunningSolution.SetClustersforaPoint(alpha);
                            NewIndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                            for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                            {
                                int CreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                                int AvailabilityIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                                /* if ((alpha + DAVectorUtility.PointStart_Process == 22625) && (ClusteringSolution.CurrentIteration > 6500))
                                {
                                    DAVectorUtility.SALSAFullPrint(0, " Redone " + IndirectClusterIndex.ToString() + " " + ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex].ToString("F4")
                                        + " Created Index " + CreatedIndex.ToString() + " Mapped " + AvailabilityIndex.ToString() );
                                } */
                                if (AvailabilityIndex > 0)
                                {
                                    --AvailabilityIndex;
                                    int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[AvailabilityIndex];
                                    ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex] = 0;
                                    ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadNo] = 0;
                                }
                                else
                                {
                                    int RemoteIndex = -AvailabilityIndex - 1;
                                    StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = 0;
                                    StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = 0;
                                }
                            }
                        }   // End case where not enough Clusters for this Point

                    }   // End Changed case

                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < NewIndirectSize; IndirectClusterIndex++)
                        ParallelClustering.RunningSolution.LegalCluster(alpha, IndirectClusterIndex);
                }
            }); // End loop Processing Point dependent quantities

            //  Unset any split information and Set Accumulation Orders and start setting NodeAccMetaData
            //  Done Sequentially in each node
            int NumberNodeClusterAccIndices = 0;
            int[] NumberThreadClusterAccIndices = new int[DAVectorUtility.ThreadCount];
            for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                NumberThreadClusterAccIndices[ThreadNo] = 0;

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                ParallelClustering.RunningSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
                if( ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] < 0 )
                    continue;
                ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex] = NumberNodeClusterAccIndices;
                NodeAccMetaData.NodeAccumulationCreatedIndices[NumberNodeClusterAccIndices] = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                ++NumberNodeClusterAccIndices;

                for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                {
                    if( ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] < 0 )
                        continue;
                    ClusteringSolution.LocalThreadAccPosition[LocalActiveClusterIndex][ThreadNo] = NumberThreadClusterAccIndices[ThreadNo];
                    ++NumberThreadClusterAccIndices[ThreadNo];
                }
            }

            for (int RemoteIndex = 0; RemoteIndex < StorageforTransportedClusters.SizeOfTransportedArray; RemoteIndex++)
            {
                StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1] = 0;
                if (StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] < 0)
                    continue;
                StorageforTransportedClusters.TransportedNodeAccPosition[RemoteIndex] = NumberNodeClusterAccIndices;
                NodeAccMetaData.NodeAccumulationCreatedIndices[NumberNodeClusterAccIndices] = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
                ++NumberNodeClusterAccIndices;

                for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                {
                    if (StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] < 0)
                        continue;
                    StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadNo] = NumberThreadClusterAccIndices[ThreadNo];
                    ++NumberThreadClusterAccIndices[ThreadNo];
                }
            }

            //  Node Accumulation Book Keeping
            NodeAccMetaData.NumberofPointsperNode = NumberNodeClusterAccIndices;
            Program.ActualMaxNumberAccumulationsperNode = Math.Max(Program.ActualMaxNumberAccumulationsperNode, NumberNodeClusterAccIndices);
            ParallelNodeAccumulationRanges = RangePartitioner.Partition(NumberNodeClusterAccIndices, DAVectorUtility.ThreadCount);
            for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                NodeAccMetaData.NumberofPointsperThread[ThreadNo] = NumberThreadClusterAccIndices[ThreadNo];

            //  Associate Thread accumulation to Node Accumulation arrays
            AssociateThreadtoNodeAccumulation();
            DAVectorUtility.StopSubTimer(11);

            DAVectorUtility.StartSubTimer(12);
            ParallelClustering.RunningSolution.SetClusterSizes();
            ParallelClustering.RunningSolution.SetClusterWidths();
            DAVectorUtility.StopSubTimer(12);

        }   // End ManageMajorSynchronization()

        //  Associate Thread Accumnulation Positions to Node Accumulation Positions
        //  Also set Status and Hosts
        //  Parallel over Node Accumulation Positions
        public static void AssociateThreadtoNodeAccumulation()
        {

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (AccumulationThreadNo) =>
            {
                int beginindex = ParallelNodeAccumulationRanges[AccumulationThreadNo].StartIndex;
                int indexlength = ParallelNodeAccumulationRanges[AccumulationThreadNo].Length;
                int NodeStoragePosition = -1;
                int TransportedStoragePosition = -1;
                int NodeAccumulationPosition = -1;
                int ThreadAccumulationPosition = -1;

                for (int NodeAccumulationIndex = beginindex; NodeAccumulationIndex < beginindex + indexlength; NodeAccumulationIndex++)
                {
                    int createdindex = NodeAccMetaData.NodeAccumulationCreatedIndices[NodeAccumulationIndex];
                    for (int ThreadNo = 0; ThreadNo < DAVectorUtility.ThreadCount; ThreadNo++)
                    {
                        int isitOK = IndicesperCluster(createdindex, ThreadNo, ref NodeStoragePosition, ref TransportedStoragePosition, ref NodeAccumulationPosition, ref ThreadAccumulationPosition);
                        if (isitOK < 0)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Inconsistent Created Index " + createdindex.ToString() + " in Node  Accumulation " + NodeAccumulationIndex.ToString());
                            throw (e);
                        }
                        if (NodeAccumulationPosition != NodeAccumulationIndex)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Inconsistent Created Index " + createdindex.ToString() + " " + DAVectorUtility.MPI_Rank.ToString() + " " + ClusteringSolution.NumberLocalActiveClusters.ToString()
                                + " in Node  Accumulation 1 " + NodeAccumulationIndex.ToString() + " and Node  Accumulation 2 " + NodeAccumulationPosition.ToString());
                            throw (e);
                        }
                        NodeAccMetaData.AccumulationNodetoThreadClusterAssociations[NodeAccumulationIndex][ThreadNo] = ThreadAccumulationPosition;  // Set Mapping between Node and Thread Storage
                    }
                    //  Set role and message structure of Node Accumulation Points
                    if (NodeStoragePosition >= 0)
                    {   // Cluster controlled locally or is a global cluster
                        int ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[NodeStoragePosition];
                        int LocalClusterStatus = ParallelClustering.RunningSolution.LocalStatus[NodeStoragePosition];
                        if ((LocalClusterStatus <= -1) || (LocalClusterStatus == 3))
                        {
                            Exception e = DAVectorUtility.SALSAError(" Created Index " + createdindex.ToString() + " in Node  Accumulation " + NodeAccumulationIndex.ToString()
                                + " and Node  Storage " + NodeStoragePosition.ToString() + " Deleted");
                            throw (e);
                        }
                        NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex] = LocalClusterStatus;    
                        NodeAccMetaData.NodeAccumulationClusterHosts[NodeAccumulationIndex] = ClusteringSolution.LocalHost[ActiveClusterIndex];
                        continue;
                    }
                    else
                    {
                        if (TransportedStoragePosition < 0)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Created Index " + createdindex.ToString() + " in Node Accumulation " + NodeAccumulationIndex.ToString()
                                + " in neither category");
                            throw (e);
                        }
                        NodeAccMetaData.NodeAccumulationClusterStatus[NodeAccumulationIndex] = 3;
                        NodeAccMetaData.NodeAccumulationClusterHosts[NodeAccumulationIndex] =
                            StorageforTransportedClusters.TotalTransportedOriginalHost[TransportedStoragePosition];  // Packed Node which controls this cluster
                    }
                }

            }); //  End Parallel Section over NodeAccumulationIndex

        }   // End AssociateThreadtoNodeAccumulation()

        public static void SetupNewIteration()
        {
            ++ClusteringSolution.CurrentIteration; // Increment Iteration Number
            int Testiteration = ClusteringSolution.CurrentIteration;
            DAVectorUtility.MPI_communicator.Broadcast<int>(ref Testiteration, 0);
            if (Testiteration != ClusteringSolution.CurrentIteration)
            {
                Exception e = DAVectorUtility.SALSAError(" Inconsistent Iteration " + Testiteration.ToString() + " " + ClusteringSolution.CurrentIteration.ToString() + " with rank " + DAVectorUtility.MPI_Rank.ToString());
                throw (e);
            }
            int NumberTransported = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                    continue;
                int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
                if ((H1 == DAVectorUtility.MPI_Rank) && (H2 == DAVectorUtility.MPI_Rank))
                    continue;
                ++NumberTransported;
            }
            if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
            {
                string endinfo = "";
                if (!Program.CalculateIndividualWidths)
                    endinfo = " Average Width " + ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth.ToString("E3");
                else
                {
                    endinfo = " Average Widths ";
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        endinfo += ParallelClustering.RunningSolution.AverageWidth[VectorIndex].ToString("E3") + " ";
                }
                string LocalMessage = "Clstrs " + ParallelClustering.RunningSolution.Ncent_ThisNode.ToString() + " Dltd " + ClusteringSolution.ClustersDeleted.ToString()
                    + " Splt " + ClusteringSolution.ClustersSplit.ToString() + " Mvd " + ClusteringSolution.ClustersMoved.ToString()
                    + " FromRmt " + StorageforTransportedClusters.SizeOfTransportedArray.ToString() + " ToRmt " + NumberTransported.ToString();
                string OverallMessage = "Major Synch Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Total Clstrs " + ParallelClustering.RunningSolution.Ncent_Global.ToString()
                    + " Temperature " + ParallelClustering.RunningSolution.Temperature.ToString("F4") + " " + endinfo;
                if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                    OverallMessage = "Sponge Count " + ParallelClustering.RunningSolution.C_k_[ParallelClustering.RunningSolution.SpongeCluster].ToString("F2") + " Factor "
                        + Program.SpongeFactor.ToString("F2") + " " + OverallMessage;
                DAVectorUtility.SALSASyncPrint(1, OverallMessage, LocalMessage);
                DAVectorUtility.SALSASyncPrint(2, "Host", HostLinkageMessage);
            }

            ClusteringSolution.ClustersDeleted = 0;
            ClusteringSolution.ClustersMoved = 0;
            ClusteringSolution.ClustersSplit = 0;
            return;

        }   // End SetupNewIteration()

        public static void ProcessChangedClusters(bool ongoing)
        {
            // Deal with deleted clusters removed as zero size
            // Process Any Clusters Moved from afar to here
            // This involves moving clusters from transported to local storagr and deleting moved clusters from local
            bool changedclusters = false;
            ClusteringSolution.ClustersDeleted = 0;
            ClusteringSolution.ClustersMoved = 0;

            //  Delete Clusters transported elsewhere or those deleted
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                    continue;
                int CreatedIndex = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0)
                {
                    changedclusters = true;
                    continue;
                }
                int NewHost = ClusteringSolution.LocalHost[RealClusterIndex] & ClusteringSolution.PACKINGMASK;
                if ( NewHost == DAVectorUtility.MPI_Rank)
                {
                    ClusteringSolution.UniversalMapping[CreatedIndex].Ymapping = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][0];
                    continue;
                }

                ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] = -1;
                ++ClusteringSolution.ClustersMoved;
                changedclusters = true;
            }
            bool GlobalChangedClusters = changedclusters;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming5);
            GlobalChangedClusters = DAVectorUtility.MPI_communicator.Allreduce<bool>(GlobalChangedClusters, Operation<bool>.LogicalOr);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming5);

            if (changedclusters)
                ParallelClustering.RunningSolution.CleanupClusters();
            if(GlobalChangedClusters)
            {
                ParallelClustering.RunningSolution.SetActiveClusters();
                for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                    ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost;
                }
            }

            //  Move Clusters from afar into local storage with status 2
            if (DAVectorUtility.MPI_Size == 1)
                return;

            int NewTransportedCount = 0;
            for (int RemoteIndex = 0; RemoteIndex < StorageforTransportedClusters.SizeOfTransportedArray; RemoteIndex++)
            {
                int RemotePackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];
                if ((RemotePackedHost & ClusteringSolution.PACKINGMASK) == DAVectorUtility.MPI_Rank)
                {   // Move this cluster to be a local distributed cluster
                    //  C_k_ must be set later
                    int NewCenterIndex = ParallelClustering.RunningSolution.Ncent_ThisNode;
                    int CreatedIndex = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
                    ParallelClustering.RunningSolution.LocalCreatedIndex[NewCenterIndex] = CreatedIndex;
                    ParallelClustering.RunningSolution.P_k_[NewCenterIndex]  = StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
                    ParallelClustering.RunningSolution.LocalStatus[NewCenterIndex] = 2;
                    ParallelClustering.RunningSolution.LocalSplitCreatedIndex[NewCenterIndex] = StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1];
                    ParallelClustering.RunningSolution.Splittable_k_[NewCenterIndex] = -1;
                    ParallelClustering.RunningSolution.SplitPriority_k_[NewCenterIndex] = 2;
                    ClusteringSolution.UniversalMapping[CreatedIndex].PackedHost = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];

                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        ParallelClustering.RunningSolution.Y_k_i_[NewCenterIndex][VectorIndex] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex];
                        ParallelClustering.RunningSolution.Sigma_k_i_[NewCenterIndex][VectorIndex] = StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex];
                    }
                    ++ParallelClustering.RunningSolution.Ncent_ThisNode;
                    continue;
                }
                if (NewTransportedCount < RemoteIndex)
                {
                    StorageforTransportedClusters.TotalTransportedCreatedIndex[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex];
                    StorageforTransportedClusters.TotalTransportedOriginalHost[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex];
                    StorageforTransportedClusters.TotalTransportedP_t[NewTransportedCount] = StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
                    for (int VectorIndex = 0; VectorIndex <= Program.ParameterVectorDimension; VectorIndex++)
                    {
                        StorageforTransportedClusters.TotalTransportedY_t_i[NewTransportedCount][VectorIndex] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex];
                        if(VectorIndex < Program.ParameterVectorDimension) 
                            StorageforTransportedClusters.TotalTransportedSigma_t_i[NewTransportedCount][VectorIndex] = StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex];
                    }
                    StorageforTransportedClusters.TotalTransportedY_t_i[NewTransportedCount][Program.ParameterVectorDimension] = StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][Program.ParameterVectorDimension];
                    for (int StatusIndex = 0; StatusIndex < 2; StatusIndex++)
                        StorageforTransportedClusters.TotalTransportedStatus[NewTransportedCount][StatusIndex] = StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][StatusIndex];
                }
                ++NewTransportedCount;
            }

            // Set new counts
            if (NewTransportedCount < StorageforTransportedClusters.SizeOfTransportedArray)
            {
                changedclusters = true;
                StorageforTransportedClusters.SizeOfTransportedArray = NewTransportedCount;
                Program.ActualMaxTransportedClusterStorage = Math.Max(NewTransportedCount + 1, Program.ActualMaxTransportedClusterStorage);
            }
            GlobalChangedClusters = changedclusters;
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming5);
            GlobalChangedClusters = DAVectorUtility.MPI_communicator.Allreduce<bool>(GlobalChangedClusters, Operation<bool>.LogicalOr);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming5);

            if (GlobalChangedClusters)
            {
                ParallelClustering.RunningSolution.SetActiveClusters();
                for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                    ClusteringSolution.LocalHost[ActiveClusterIndex] = ClusteringSolution.UniversalMapping[ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]].PackedHost;
                }
            }

        }   // ProcessMovedClusters()

        //  Switch from global to distributed mode
        //  Note later ProcessMovedClusters() will delete clusters not hosted on each node
        //  This routine sets in StorageforTransportedClusters array all clusters NOT hosted here but relevant
        public static void SetupDistributedMode()
        {
            int RemoteIndex = 0;
            int LocalType2 = 0;
            int SendUp = 0;
            int SendDown = 0;
            int Type01 = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                {
                    Type01++;
                    continue;
                }
                int PackedHost = ClusteringSolution.LocalHost[ActiveClusterIndex];
                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
                if (H == DAVectorUtility.MPI_Rank)
                {
                    ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] = 2;
                    ++LocalType2;
                    if( H2 > H)
                        ++SendUp;
                    if (H1 < H)
                        ++SendDown;
                    continue;
                }

                if(DAVectorUtility.MPI_Rank < H1)
                    continue;
                if(DAVectorUtility.MPI_Rank > H2)
                    continue;

                // A Cluster that affects this node but not hosted here
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0)
                    continue;
                if (RemoteIndex >= Program.MaxMPITransportBuffer)
                {
                    Exception e = DAVectorUtility.SALSAError("MPI Buffer over limit " + RemoteIndex.ToString() + " Limit " + Program.MaxMPITransportBuffer
                        + " Clusters " + ParallelClustering.RunningSolution.Ncent_Global.ToString() + " Local " + ParallelClustering.RunningSolution.Ncent_ThisNode.ToString());
                    throw (e);
                }
                StorageforTransportedClusters.TotalTransportedCreatedIndex[RemoteIndex] = ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex];
                StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex] = ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                StorageforTransportedClusters.TotalTransportedOriginalHost[RemoteIndex] = PackedHost;
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                    StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex][VectorIndex] = ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex][VectorIndex];
                }
                StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex][Program.ParameterVectorDimension] = ParallelClustering.RunningSolution.P_k_[RealClusterIndex]; ;

                StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][0] = 3;
                StorageforTransportedClusters.TotalTransportedStatus[RemoteIndex][1] = ParallelClustering.RunningSolution.LocalSplitCreatedIndex[RealClusterIndex];
                ++RemoteIndex;
            }   // End Loop over Clusters being dispersed

            StorageforTransportedClusters.SizeOfTransportedArray = RemoteIndex;
            Program.ActualMaxTransportedClusterStorage = Math.Max(RemoteIndex + 1, Program.ActualMaxTransportedClusterStorage);

            int[] ClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
            int[] UpClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
            int[] DownClusterCountsperNode = new int[DAVectorUtility.MPI_Size];
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIGATHERTiming);
            ClusterCountsperNode = DAVectorUtility.MPI_communicator.Allgather<int>(LocalType2);
            UpClusterCountsperNode = DAVectorUtility.MPI_communicator.Allgather<int>(SendUp);
            DownClusterCountsperNode = DAVectorUtility.MPI_communicator.Allgather<int>(SendDown);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIGATHERTiming);
            string message = "Global " + Type01.ToString();
            int TotalClusterCount = Type01;
            for (int nodeindex = 0; nodeindex < DAVectorUtility.MPI_Size; nodeindex++)
            {
                TotalClusterCount += ClusterCountsperNode[nodeindex];
                message += " Node " + nodeindex.ToString() + ":" + ClusterCountsperNode[nodeindex].ToString() + "(Up:" + UpClusterCountsperNode[nodeindex].ToString()
                    + " Down:" + DownClusterCountsperNode[nodeindex].ToString() +")";
            }
            DAVectorUtility.SALSAPrint(1, "Clusters Distributed " + message);
            if (TotalClusterCount != ParallelClustering.RunningSolution.Ncent_Global)
            {
                Exception e = DAVectorUtility.SALSAError("Inconsistent Cluster Counts " + TotalClusterCount.ToString() + " " + ParallelClustering.RunningSolution.Ncent_Global.ToString());
                throw(e);
            }

        }   //  End SetupDistributedMode()
        
    }   // End DistributedClusteringSolution

}   // End Namespace Salsa.DAVectorSponge