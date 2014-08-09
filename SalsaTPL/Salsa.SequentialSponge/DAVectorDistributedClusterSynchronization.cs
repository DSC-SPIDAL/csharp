using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MPI;
using Salsa.DAVectorSponge;
using SALSALibrary;

namespace SALSALibrary
{

    [Serializable()]
    public class MPITransportComponentPacket
    {
        public int NumberofClusters;  // Number of Clusters in Package
        public int NumberofDoubleComponents;  // Number of Double Components in Package
        public int NumberofIntegerComponents;  // Number of Integer Components in Package
        public double[] ClusterDoubleComponents;  // Transported Double Components
        public int[] ClusterIntegerComponents;  // Transported Integer Components
        public int[] AssociatedCreatedIndex;    // Created Index of Cluster
        public int[] ClusterHostRange;   // Hosting node Intended for Cluster in form H + MULTIPLIER (H1 + MULTIPLIER * H2 )
        //   Note H1 <= H <= H2
        //  Near End either H1 = H2 = H or H1 = H - 1 and/or H2 = H + 1

        public MPITransportComponentPacket(int _NumberofClusters, int _NumberofDoubleComponents, int _NumberofIntegerComponents)
        {
            int Maxlength = _NumberofClusters * _NumberofDoubleComponents;
            Maxlength = Math.Max(Maxlength, 1);
            ClusterDoubleComponents = new double[Maxlength];

            Maxlength = _NumberofClusters * _NumberofIntegerComponents;
            Maxlength = Math.Max(Maxlength, 1);
            ClusterIntegerComponents = new int[Maxlength];

            Maxlength = Math.Max(_NumberofClusters, 1);
            AssociatedCreatedIndex = new int[Maxlength];
            ClusterHostRange = new int[Maxlength];
            NumberofClusters = _NumberofClusters;
            NumberofDoubleComponents = _NumberofDoubleComponents;
            NumberofIntegerComponents = _NumberofIntegerComponents;
        }
    }

    public class DistributedSynchronization
    {   // Do the pipeline distributed broadcast with several different options
 
        public static MPITransportComponentPacket TransportComponent;   // Place to receive transported components

        public class TransportviaPipeline
        {
            private int InitialArraySize;
            private int NumberofDoubleComponents;
            private int NumberofIntegerComponents;
            private int UpdateMode;   // = 0 Replace in place; used to send round new values of cluster parameters;
            // UpdateMode = 1 Increment in place; used in distributed reduction;
            // UpdateMode = 2 Gather targeted clusters when setting up remote clusters and propagating deleted, moved and split(new) clusters
            private int StorageMode;    // = 0 use with UpdateMode=2; 
            // StorageMode = 1 Local Cluster Storage (NodeStoragePosition) NOT USED
            // StorageMode = 2 Remote Clusters Stored Locally (TransportedStoragePosition)
            // StorageMode = 3 Node Accumulation Array (NodeAccumulationPosition)
            private bool HostRangeProcessing;   // True use Host Range for destination; False Send only to host

            public TransportviaPipeline(int _UpdateMode, bool _HostRangeProcessing, int _NumberofDoubleComponents, 
                int _NumberofIntegerComponents, int _InitialArraySize, int _StorageMode )
            {
                NumberofDoubleComponents = _NumberofDoubleComponents;
                NumberofIntegerComponents = _NumberofIntegerComponents;
                UpdateMode = _UpdateMode;
                InitialArraySize = _InitialArraySize;
                HostRangeProcessing = _HostRangeProcessing;
                StorageMode = _StorageMode;
            }

            //  Note FinalClusterCount, FinalCreatedIndex and FinalHostSpecification ONLY used in UpdateMode 2 and only set in this case
            public void PipelineDistributedBroadcast(double[][] InitialDoubleComponents, double[][] FinalDoubleComponents, int[][] InitialIntegerComponents, int[][] FinalIntegerComponents, 
                int[] InitialCreatedIndex, int[] FinalCreatedIndex, int[] InitialHostSpecification, int[] FinalHostSpecification, ref int FinalClusterCount )
            {
                FinalClusterCount = 0;
                if (DAVectorUtility.MPI_Size <= 1)
                    return;
                    
                // Now process distributed clusters
                // Variables for processing createdindex
                int NodeStoragePosition = -1;
                int TransportedStoragePosition = -1;
                int NodeAccumulationPosition = -1;
                int ThreadAccumulationPosition = -1;

                //  Place where received data stored
                int FinalDataLocationIndex = -1;

                ++Program.NumberPipelineGroups; // Increment calls of this routine

                int[] DownbySteps = new int[DAVectorUtility.MPI_Size];
                int[] UpbySteps = new int[DAVectorUtility.MPI_Size];
                int[] DownbyStepsTotal = new int[DAVectorUtility.MPI_Size];
                int[] UpbyStepsTotal = new int[2*DAVectorUtility.MPI_Size];
                for (int PipelineSteps = 0; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++)
                {
                    DownbySteps[PipelineSteps] = 0;
                    UpbySteps[PipelineSteps] = 0;
                }

                //  Set NumberUp and NumberDown
                for (int ClusterIndirectIndex = 0; ClusterIndirectIndex < InitialArraySize; ClusterIndirectIndex++)
                {
                    int PackedHost = InitialHostSpecification[ClusterIndirectIndex];
                    for (int PipelineSteps = 1; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++)
                    {
                        if (HostRangeProcessing)
                        {
                            int H1 = PackedHost >> ClusteringSolution.PACKINGSHIFT;
                            int H2 = H1 >> ClusteringSolution.PACKINGSHIFT;
                            H1 = H1 & ClusteringSolution.PACKINGMASK;
                            if (H2 > (DAVectorUtility.MPI_Rank + PipelineSteps - 1) )
                                ++UpbySteps[PipelineSteps];
                            if (H1 < (DAVectorUtility.MPI_Rank -PipelineSteps + 1) )
                                ++DownbySteps[PipelineSteps];
                        }
                        else
                        {
                            int H = PackedHost & ClusteringSolution.PACKINGMASK;
                            if (H > (DAVectorUtility.MPI_Rank + PipelineSteps - 1) )
                                ++UpbySteps[PipelineSteps];
                            if (H < (DAVectorUtility.MPI_Rank - PipelineSteps + 1) )
                                ++DownbySteps[PipelineSteps];
                        }
                    }
                }
                for (int PipelineSteps = 0; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++)
                {
                    UpbyStepsTotal[PipelineSteps] = UpbySteps[PipelineSteps];
                    UpbyStepsTotal[PipelineSteps + DAVectorUtility.MPI_Size] = DownbySteps[PipelineSteps];
                }
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming4);
                UpbyStepsTotal = DAVectorUtility.MPI_communicator.Allreduce<int>(UpbyStepsTotal, Operation<int>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming4);
                for (int PipelineSteps = 0; PipelineSteps < DAVectorUtility.MPI_Size; PipelineSteps++)
                    DownbyStepsTotal[PipelineSteps] = UpbyStepsTotal[PipelineSteps + DAVectorUtility.MPI_Size];


                // Variables Used for Up and Down Sections
                bool Initialstep;
                int CurrentNode = DAVectorUtility.MPI_Rank;
                int ReceivedTotal = 0;
                int NumberClustertoSend = 0;
                int NumberDoubletoSend = 0;
                int NumberIntegertoSend = 0;
                int IsItOK;

                // Process Clusters going Up the Chain
                Initialstep = true;
                int LocalTotal = UpbySteps[1]; ;
                int StepsUp = 0;

                while(true)
                {   
                    // Decide if ANY node needs to communicate Up
                    ++StepsUp;
                    if (StepsUp >= DAVectorUtility.MPI_Size)
                        break;
                    int JobTotal = UpbyStepsTotal[StepsUp];
                    if(JobTotal == 0)
                        break;
                    

                    // Some Nodes want to go up the line
                    int SourceProc = MPI.Intercommunicator.Null;
                    int DestProc = MPI.Intercommunicator.Null;
                    int SourceTag = 0; // Random Number
                    int DestTag = 0;
                    SourceProc = DAVectorUtility.MPI_Size - 1;
                    DestProc = 0;
                    if( CurrentNode != 0)
                        SourceProc = CurrentNode - 1;
                    if (CurrentNode != (DAVectorUtility.MPI_Size - 1))
                        DestProc = CurrentNode + 1;
                    else
                        LocalTotal = 0;
                    MPITransportComponentPacket SendBuffer = new MPITransportComponentPacket(LocalTotal, NumberofDoubleComponents, NumberofIntegerComponents); // Sent Buffer is EXACT size
                    NumberClustertoSend = 0;
                    NumberDoubletoSend = 0;
                    NumberIntegertoSend = 0;

                    if( LocalTotal > 0 )
                    {   // If no data here, just send dummy packet
                        if(Initialstep)
                        {   // Construct message to send from Initial Arrays
                            for (int ClusterSendPointer = 0; ClusterSendPointer < InitialArraySize; ClusterSendPointer++ )
                            {
                                int PackedHost = InitialHostSpecification[ClusterSendPointer];
                                if (HostRangeProcessing)
                                {
                                    int H2 = PackedHost >> (2 * ClusteringSolution.PACKINGSHIFT);
                                    if (H2 <= DAVectorUtility.MPI_Rank)
                                        continue;
                                }
                                else
                                {
                                    int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                    if (H <= DAVectorUtility.MPI_Rank)
                                        continue;
                                }
                                SendBuffer.AssociatedCreatedIndex[NumberClustertoSend] = InitialCreatedIndex[ClusterSendPointer];
                                SendBuffer.ClusterHostRange[NumberClustertoSend] = InitialHostSpecification[ClusterSendPointer];
                                if (NumberofDoubleComponents > 0)
                                {
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterDoubleComponents[NumberDoubletoSend] = InitialDoubleComponents[ClusterSendPointer][ComponentIndex];
                                        ++NumberDoubletoSend;
                                    }
                                }
                                if (NumberofIntegerComponents > 0)
                                {
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterIntegerComponents[NumberIntegertoSend] = InitialIntegerComponents[ClusterSendPointer][ComponentIndex];
                                        ++NumberIntegertoSend;
                                    }
                                }
                                ++NumberClustertoSend;
                                if (NumberClustertoSend >= LocalTotal)
                                    break;
                            }
                        }
                        else
                        {   // Construct message to send from en passant data
                            for(int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++)
                            {
                                int PackedHost = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                                if(HostRangeProcessing)
                                {
                                    int H2 = PackedHost >> (2*ClusteringSolution.PACKINGSHIFT);
                                    if( H2 <= DAVectorUtility.MPI_Rank)
                                        continue;
                                }
                                else
                                {
                                    int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                    if( H <= DAVectorUtility.MPI_Rank)
                                        continue;
                                }  
                                SendBuffer.AssociatedCreatedIndex[NumberClustertoSend] = TransportComponent.AssociatedCreatedIndex[ReceivedClusterIndex];
                                SendBuffer.ClusterHostRange[NumberClustertoSend] = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                                if(NumberofDoubleComponents > 0 )
                                {
                                    int OverallIndex = NumberofDoubleComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterDoubleComponents[NumberDoubletoSend] = TransportComponent.ClusterDoubleComponents[OverallIndex];
                                        ++NumberDoubletoSend;
                                        ++OverallIndex;
                                    }
                                }
                                if(NumberofIntegerComponents > 0 )
                                {
                                    int OverallIndex = NumberofIntegerComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterIntegerComponents[NumberIntegertoSend] = TransportComponent.ClusterIntegerComponents[OverallIndex];
                                        ++NumberIntegertoSend;
                                        ++OverallIndex;
                                    }
                                }
                                ++NumberClustertoSend;
                                if(NumberClustertoSend >= LocalTotal)
                                    break;
                            }

                        }
                    }   // End Case where there is Local Data to Send

                    // Send data in a pipeline forward
                    DAVectorUtility.StartSubTimer(DAVectorUtility.MPISENDRECEIVETiming);
                    DAVectorUtility.MPI_communicator.SendReceive<MPITransportComponentPacket>(SendBuffer, DestProc, DestTag, SourceProc, SourceTag, 
                        out TransportComponent, out DistributedClusteringSolution.MPISecStatus);
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPISENDRECEIVETiming);

                    ++Program.NumberPipelineSteps;
                    Program.NumberofPipelineClusters += SendBuffer.NumberofClusters;

                    //  Examine Data passed from lower ranked processor
                    //  Set new LocalTotal and Store Data
                    ReceivedTotal = TransportComponent.NumberofClusters;
                    Program.ActualMaxMPITransportBuffer = Math.Max(Program.ActualMaxMPITransportBuffer, ReceivedTotal);
                    LocalTotal = 0; // Count Number of Clusters on next step

                    if(NumberofDoubleComponents != TransportComponent.NumberofDoubleComponents)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Double Components Inconsistent " + NumberofDoubleComponents.ToString() + " " + TransportComponent.NumberofDoubleComponents.ToString()
                            + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad");
                        throw (e);
                    }
                    if(NumberofIntegerComponents != TransportComponent.NumberofIntegerComponents)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Integer Components Inconsistent " + NumberofIntegerComponents.ToString() + " " + TransportComponent.NumberofIntegerComponents.ToString()
                            + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad");
                        throw (e);
                    }

                    if(ReceivedTotal > 0)
                    {
                        for(int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++)
                        {
                            int PackedHost = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                            if(HostRangeProcessing)
                            {
                                int H2 = PackedHost >> (2*ClusteringSolution.PACKINGSHIFT);
                                if(H2 < DAVectorUtility.MPI_Rank)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported host " + PackedHost.ToString() + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad Up Range");
                                    throw (e);
                                }
                                if( H2 > DAVectorUtility.MPI_Rank)
                                    ++LocalTotal;
                            }
                            else
                            {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if(H < DAVectorUtility.MPI_Rank)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported host " + PackedHost.ToString() + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad Up Host");
                                    throw (e);
                                }
                                if( H > DAVectorUtility.MPI_Rank)
                                {
                                    ++LocalTotal;
                                    continue;
                                }
                            }
                            int host = PackedHost & ClusteringSolution.PACKINGMASK;
                            int CreatedIndex = TransportComponent.AssociatedCreatedIndex[ReceivedClusterIndex];
                            if(UpdateMode < 2)
                            {
                                FinalDataLocationIndex = -1;
                                IsItOK = DistributedClusteringSolution.IndicesperCluster(CreatedIndex, -1, ref NodeStoragePosition, ref TransportedStoragePosition, 
                                    ref NodeAccumulationPosition, ref ThreadAccumulationPosition);
                                if( StorageMode == 1)
                                    FinalDataLocationIndex = NodeStoragePosition;
                                if( StorageMode == 2)
                                    FinalDataLocationIndex = TransportedStoragePosition;
                                if( StorageMode == 3)
                                    FinalDataLocationIndex = NodeAccumulationPosition;
                                if( (host == DAVectorUtility.MPI_Rank) && (IsItOK != 0))
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported Created Index " + CreatedIndex.ToString() + " in Rank "
                                        + DAVectorUtility.MPI_Rank.ToString() + " Bad with code " + IsItOK.ToString() + " host " + host.ToString() + " Update mode " + UpdateMode.ToString());
                                    throw (e);
                                }
                            }
                            else
                            {   // UpdateMode 2
                                IsItOK = 0;
                                FinalDataLocationIndex = FinalClusterCount;
                                ++FinalClusterCount;
                                FinalCreatedIndex[FinalDataLocationIndex] = CreatedIndex;
                                FinalHostSpecification[FinalDataLocationIndex] = PackedHost;
                            }
                            if(IsItOK >= 0 )
                            {
                                if( FinalDataLocationIndex == -1)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported Created Index " + CreatedIndex.ToString() + " in Rank "
                                        + DAVectorUtility.MPI_Rank.ToString() + " Bad with Storage Mode " + StorageMode.ToString() + " host " + host.ToString() + " Update mode " + UpdateMode.ToString());
                                    throw (e);
                                }
                                
                                if(NumberofDoubleComponents > 0)
                                {
                                    string message = "";
                                    int OverallDoubleIndex = NumberofDoubleComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        if( (UpdateMode == 0) || (UpdateMode == 2))
                                            FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.ClusterDoubleComponents[OverallDoubleIndex];
                                        if(UpdateMode == 1)
                                            FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.ClusterDoubleComponents[OverallDoubleIndex];
                                        message += " * " + TransportComponent.ClusterDoubleComponents[OverallDoubleIndex].ToString("E3") + " " + FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex].ToString("F3");
                                        ++OverallDoubleIndex;
                                    }
                                    /* if (CreatedIndex == 901)
                                        DAVectorUtility.SALSAFullPrint(1, "Up901 Transport " + UpdateMode.ToString() + " " + FinalDataLocationIndex.ToString()  + message); */
                                }
                                
                                if(NumberofIntegerComponents > 0)
                                {
                                    int OverallIntegerIndex = NumberofIntegerComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        if( (UpdateMode == 0) || (UpdateMode == 2))
                                            FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.ClusterIntegerComponents[OverallIntegerIndex];
                                        if(UpdateMode == 1)
                                            FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.ClusterIntegerComponents[OverallIntegerIndex];
                                        ++OverallIntegerIndex;
                                    }
                                }
                            }   // End case where location found IsItOK >= 0
                        }   // end Loop over ReceivedClusterIndex
                    }   // End case where ReceivedTotal > 0
                    Initialstep = false;

                }   // End While over MPI pipeline steps for pipeline going UP the chain

                // Process Clusters going Down the Chain
                Initialstep = true;
                int StepsDown = 0;
                LocalTotal = DownbySteps[1];

                while(true)
                {
                    StepsDown++;
                    if (StepsDown >= DAVectorUtility.MPI_Size)
                        break;
                    int JobTotal = DownbyStepsTotal[StepsDown];
                    if(JobTotal == 0)
                        break;

                    // Some Nodes want to go down the line
                    int SourceProc = MPI.Intercommunicator.Null;
                    int DestProc = MPI.Intercommunicator.Null;
                    DestProc = DAVectorUtility.MPI_Size - 1;
                    SourceProc = 0;
                    int SourceTag = 22; // Random Number
                    int DestTag = 22;
                    if (CurrentNode != 0)
                        DestProc = CurrentNode - 1;
                    else
                        LocalTotal = 0;
                    if( CurrentNode != (DAVectorUtility.MPI_Size-1))
                        SourceProc = CurrentNode + 1;
                    MPITransportComponentPacket SendBuffer = new MPITransportComponentPacket(LocalTotal, NumberofDoubleComponents, NumberofIntegerComponents); // Sent Buffer is EXACT size
                    NumberClustertoSend = 0;
                    NumberDoubletoSend = 0;
                    NumberIntegertoSend = 0;

                    if( LocalTotal > 0 )
                    {   // If no data here, just send dummy packet

                        if(Initialstep)
                        {   // Construct message to send from local accumulation arrays
                            for (int ClusterIndirectIndex = 0; ClusterIndirectIndex < InitialArraySize; ClusterIndirectIndex++)
                            {
                                int PackedHost = InitialHostSpecification[ClusterIndirectIndex];
                                if(HostRangeProcessing)
                                {
                                    int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                                    if( H1 >= DAVectorUtility.MPI_Rank)
                                        continue;
                                }
                                else
                                {
                                    int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                    if( H >= DAVectorUtility.MPI_Rank)
                                        continue;
                                }
                                SendBuffer.AssociatedCreatedIndex[NumberClustertoSend] = InitialCreatedIndex[ClusterIndirectIndex];
                                SendBuffer.ClusterHostRange[NumberClustertoSend] = InitialHostSpecification[ClusterIndirectIndex];
                                if(NumberofDoubleComponents > 0 )
                                {
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterDoubleComponents[NumberDoubletoSend] = InitialDoubleComponents[ClusterIndirectIndex][ComponentIndex];
                                        ++NumberDoubletoSend;
                                    }
                                }
                                if(NumberofIntegerComponents > 0 )
                                {
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterIntegerComponents[NumberIntegertoSend] = InitialIntegerComponents[ClusterIndirectIndex][ComponentIndex];
                                        ++NumberIntegertoSend;
                                    }
                                }
                                ++NumberClustertoSend;
                                if(NumberClustertoSend >= LocalTotal)
                                    break;
                            }
                        }
                        else
                        {   // Construct message to send from en passant data
                            for(int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++)
                            {
                                int PackedHost = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                                if(HostRangeProcessing)
                                {
                                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                                if( H1 >= DAVectorUtility.MPI_Rank)
                                    continue;
                                }
                                else
                                {
                                    int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                    if( H >= DAVectorUtility.MPI_Rank)
                                        continue;
                                }  
                                SendBuffer.AssociatedCreatedIndex[NumberClustertoSend] = TransportComponent.AssociatedCreatedIndex[ReceivedClusterIndex];
                                SendBuffer.ClusterHostRange[NumberClustertoSend] = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                                if(NumberofDoubleComponents > 0 )
                                {
                                    int OverallIndex = NumberofDoubleComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterDoubleComponents[NumberDoubletoSend] = TransportComponent.ClusterDoubleComponents[OverallIndex];
                                        ++NumberDoubletoSend;
                                        ++OverallIndex;
                                    }
                                }
                                if(NumberofIntegerComponents > 0 )
                                {
                                    int OverallIndex = NumberofIntegerComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        SendBuffer.ClusterIntegerComponents[NumberIntegertoSend] = TransportComponent.ClusterIntegerComponents[OverallIndex];
                                        ++NumberIntegertoSend;
                                        ++OverallIndex;
                                    }
                                }
                                ++NumberClustertoSend;
                                if(NumberClustertoSend >= LocalTotal)
                                    break;
                            }
                        }   // end en passant data is source of information

                    }   // End Case where there is Local Data to Send

                    // Send data in a pipeline backwards
                    DAVectorUtility.StartSubTimer(DAVectorUtility.MPISENDRECEIVETiming);
                    DAVectorUtility.MPI_communicator.SendReceive<MPITransportComponentPacket>(SendBuffer, DestProc, DestTag, SourceProc, SourceTag, out TransportComponent, out DistributedClusteringSolution.MPISecStatus);;
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPISENDRECEIVETiming);
                    ++Program.NumberPipelineSteps;
                    Program.NumberofPipelineClusters += SendBuffer.NumberofClusters;

                    //  Examine Data passed from higher ranked processor
                    ReceivedTotal = TransportComponent.NumberofClusters;
                    Program.ActualMaxMPITransportBuffer = Math.Max(Program.ActualMaxMPITransportBuffer, ReceivedTotal);
                    LocalTotal = 0;

                    if(NumberofDoubleComponents != TransportComponent.NumberofDoubleComponents)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Double Components Inconsistent " + NumberofDoubleComponents.ToString() + " " + TransportComponent.NumberofDoubleComponents.ToString()
                            + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad");
                        throw (e);
                    }
                    if(NumberofIntegerComponents != TransportComponent.NumberofIntegerComponents)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Integer Components Inconsistent " + NumberofIntegerComponents.ToString() + " " + TransportComponent.NumberofIntegerComponents.ToString()
                            + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad");
                        throw (e);
                    }

                    if(ReceivedTotal > 0)
                    {
                        for(int ReceivedClusterIndex = 0; ReceivedClusterIndex < ReceivedTotal; ReceivedClusterIndex++)
                        {
                            int PackedHost = TransportComponent.ClusterHostRange[ReceivedClusterIndex];
                            if(HostRangeProcessing)
                            {
                                int H1 = (PackedHost >> ClusteringSolution.PACKINGSHIFT) & ClusteringSolution.PACKINGMASK;
                                if(H1 > DAVectorUtility.MPI_Rank)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported host " + PackedHost.ToString() + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad Down Range");
                                    throw (e);
                                }
                                if( H1 < DAVectorUtility.MPI_Rank)
                                    ++LocalTotal;
                            }
                            else
                            {
                                int H = PackedHost & ClusteringSolution.PACKINGMASK;
                                if(H > DAVectorUtility.MPI_Rank)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported host " + PackedHost.ToString() + " in Rank " + DAVectorUtility.MPI_Rank.ToString() + " Bad Down Not Range");
                                    throw (e);
                                }
                                if( H < DAVectorUtility.MPI_Rank)
                                {
                                    ++LocalTotal;
                                    continue;
                                }
                            }
                            int host = PackedHost & ClusteringSolution.PACKINGMASK;
                            int CreatedIndex = TransportComponent.AssociatedCreatedIndex[ReceivedClusterIndex];
                            if(UpdateMode < 2)
                            {
                                FinalDataLocationIndex = -1;
                                IsItOK = DistributedClusteringSolution.IndicesperCluster(CreatedIndex, -1, ref NodeStoragePosition, ref TransportedStoragePosition, 
                                    ref NodeAccumulationPosition, ref ThreadAccumulationPosition);
                                if( StorageMode == 1)
                                    FinalDataLocationIndex = NodeStoragePosition;
                                if( StorageMode == 2)
                                    FinalDataLocationIndex = TransportedStoragePosition;
                                if( StorageMode == 3)
                                    FinalDataLocationIndex = NodeAccumulationPosition;
                                if( (host == DAVectorUtility.MPI_Rank) && (IsItOK != 0))
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported Created Index " + CreatedIndex.ToString() + " in Rank "
                                        + DAVectorUtility.MPI_Rank.ToString() + " Bad with code " + IsItOK.ToString() + " host " + host.ToString() + " Update mode " + UpdateMode.ToString());
                                    throw (e);
                                }
                            }
                            else
                            {
                                IsItOK = 0;
                                FinalDataLocationIndex = FinalClusterCount;
                                ++FinalClusterCount;
                                FinalCreatedIndex[FinalDataLocationIndex] = CreatedIndex;
                                FinalHostSpecification[FinalDataLocationIndex] = PackedHost;
                            }
                            if(IsItOK >= 0 )
                            {
                                if( FinalDataLocationIndex == -1)
                                {
                                    Exception e = DAVectorUtility.SALSAError(" Transported Created Index " + CreatedIndex.ToString() + " in Rank "
                                        + DAVectorUtility.MPI_Rank.ToString() + " Bad with Storage Mode " + StorageMode.ToString() + " host " + host.ToString() + " Update mode " + UpdateMode.ToString());
                                    throw (e);
                                }
                                
                                if(NumberofDoubleComponents > 0)
                                {
                                    string message = "";
                                    int OverallDoubleIndex = NumberofDoubleComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofDoubleComponents; ComponentIndex++)
                                    {
                                        if( (UpdateMode == 0) || (UpdateMode == 2))
                                            FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.ClusterDoubleComponents[OverallDoubleIndex];
                                        if( UpdateMode == 1 )
                                            FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.ClusterDoubleComponents[OverallDoubleIndex];
                                        message += " * " + TransportComponent.ClusterDoubleComponents[OverallDoubleIndex].ToString("E3") + " " + FinalDoubleComponents[FinalDataLocationIndex][ComponentIndex].ToString("F3");
                                        ++OverallDoubleIndex;
                                    }
                                    /* if (CreatedIndex == 901)
                                        DAVectorUtility.SALSAFullPrint(1, "Dn901 Transport " + UpdateMode.ToString() + " " + FinalDataLocationIndex.ToString() + message); */
                                }
                                
                                if(NumberofIntegerComponents > 0)
                                {
                                    int OverallIntegerIndex = NumberofIntegerComponents*ReceivedClusterIndex;
                                    for (int ComponentIndex = 0; ComponentIndex < NumberofIntegerComponents; ComponentIndex++)
                                    {
                                        if( (UpdateMode == 0) || (UpdateMode == 2))
                                            FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] = TransportComponent.ClusterIntegerComponents[OverallIntegerIndex];
                                        if(UpdateMode == 1)
                                            FinalIntegerComponents[FinalDataLocationIndex][ComponentIndex] += TransportComponent.ClusterIntegerComponents[OverallIntegerIndex];
                                        ++OverallIntegerIndex;
                                    }
                                }
                            }
                        }   // End Loop over Received Clusters
                    }   // End case when data received ReceivedTotal > 0

                    Initialstep = false;

                }   // End While over MPI pipeline steps going DOWN the chain
                
            }   // End PipelineDistributedBroadcast

        }   // End TransportviaPipeline

    }   // End class DistributedSynchronization

}   // End Namespace SALSALIbrary
