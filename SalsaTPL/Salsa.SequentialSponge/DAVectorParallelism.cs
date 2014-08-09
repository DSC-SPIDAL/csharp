using System;
using System.IO;
using MPI;
using Salsa.Core.Blas;
using Salsa.Core;

using Salsa.DAVectorSponge;

namespace SALSALibrary
{
    public class DAVectorParallelism
    {
        public static void SetupParallelism(ref string[] args)
        {
            //  Set up MPI
            DAVectorUtility.MPI_Environment = new MPI.Environment(ref args);
            DAVectorUtility.MPI_communicator = Communicator.world;         //initializing MPI world communicator
            DAVectorUtility.MPI_Rank = DAVectorUtility.MPI_communicator.Rank; // Rank of this process
            DAVectorUtility.MPI_Size = DAVectorUtility.MPI_communicator.Size; // Number of MPI Processes

            // Set up MPI
            DAVectorUtility.MPIperNodeCount = DAVectorUtility.MPI_Size / DAVectorUtility.NodeCount;

            if ((DAVectorUtility.MPIperNodeCount * DAVectorUtility.NodeCount) != DAVectorUtility.MPI_Size)
            {
                Exception e = DAVectorUtility.SALSAError("Inconsistent MPI counts Nodes "
                    + DAVectorUtility.NodeCount.ToString() + " Size " + DAVectorUtility.MPI_Size.ToString());

                throw (e);
            }

            DAVectorUtility.ParallelPattern = "---------------------------------------------------------\nMachine:" 
                + MPI.Environment.ProcessorName.ToString() + " " + DAVectorUtility.ThreadCount.ToString() 
                + "x" + DAVectorUtility.MPIperNodeCount.ToString() + "x" + DAVectorUtility.NodeCount.ToString();
            if (DAVectorUtility.MPI_Rank == 0)
                DAVectorUtility.SALSAPrint(0, DAVectorUtility.ParallelPattern);

        }   // End SetupParallelism

        public static void TearDownParallelism()
        {
            if (DAVectorUtility.MPI_Environment != null)
            {
                DAVectorUtility.MPI_Environment.Dispose();
                DAVectorUtility.MPI_Environment = null;
            }
        }   // End TearDownParallelism

        public static void SetParallelDecomposition()
        {
            //	First divide points among processes
            Range[] processRanges = RangePartitioner.Partition(DAVectorUtility.PointCount_Global, DAVectorUtility.MPI_Size);
            Range processRange = processRanges[DAVectorUtility.MPI_Rank];  // The answer for this process

            DAVectorUtility.PointStart_Process = processRange.StartIndex;
            DAVectorUtility.PointCount_Process = processRange.Length;
            DAVectorUtility.PointCount_Largest = int.MinValue;

            foreach (Range r in processRanges)
            {
                DAVectorUtility.PointCount_Largest = Math.Max(r.Length, DAVectorUtility.PointCount_Largest);
            }

            // We need points per process for all processes as used by load balancing algorithm and then to reorder distance matrix consistently
            DAVectorUtility.PointsperProcess = new int[DAVectorUtility.MPI_Size];
            DAVectorUtility.PointsperThreadperProcess = new int[DAVectorUtility.MPI_Size][];

            for (int i = 0; i < DAVectorUtility.MPI_Size; i++)
            {
                DAVectorUtility.PointsperProcess[i] = processRanges[i].Length;
                Range[] threadprocessRanges = RangePartitioner.Partition(processRanges[i], DAVectorUtility.ThreadCount);
                DAVectorUtility.PointsperThreadperProcess[i] = new int[DAVectorUtility.ThreadCount];
                for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
                {
                    DAVectorUtility.PointsperThreadperProcess[i][j] = threadprocessRanges[j].Length;
                }
            }

            //	Now divide points among threads for this process
            Range[] threadRanges = RangePartitioner.Partition(processRanges[DAVectorUtility.MPI_Rank], DAVectorUtility.ThreadCount);
            DAVectorUtility.PointsperThread = new int[DAVectorUtility.ThreadCount];
            DAVectorUtility.StartPointperThread = new int[DAVectorUtility.ThreadCount];

            for (int j = 0; j < DAVectorUtility.ThreadCount; j++)
            {
                DAVectorUtility.PointsperThread[j] = threadRanges[j].Length;
                DAVectorUtility.StartPointperThread[j] = threadRanges[j].StartIndex;
            }
        }   // End SetParallelDecomposition()

        public static double getDistancefromPoint(int LocalToProcessIndex, double[] ClusterPosition)
        {
            double tmp0 = Program.PointPosition[LocalToProcessIndex][0] - ClusterPosition[0];
            double tmp1 = Program.PointPosition[LocalToProcessIndex][1] - ClusterPosition[1];
            double tmp = tmp0 * tmp0 + tmp1 * tmp1;
            if ( Program.ParameterVectorDimension == 2 )
                return Math.Sqrt(tmp);
            if (Program.ParameterVectorDimension == 3)
            {
                double tmp2 = Program.PointPosition[LocalToProcessIndex][2] - ClusterPosition[2];
                return Math.Sqrt(tmp + tmp2 * tmp2);
            }
            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                double tmp3 = Program.PointPosition[LocalToProcessIndex][VectorIndex] - ClusterPosition[VectorIndex];
                tmp += tmp3 * tmp3;
            }

            return Math.Sqrt(tmp);

        }   // End getClusterDistancefromPoint

        public static double getSquaredScaledDistancePointActiveCluster(int LocalToProcessIndex, int ActiveClusterIndex, ClusteringSolution Solution)
        {
            if (Solution.DistributedExecutionMode && (ActiveClusterIndex >= ClusteringSolution.NumberLocalActiveClusters) )
            {
                int RemoteIndex = ActiveClusterIndex - ClusteringSolution.NumberLocalActiveClusters;
                return getSquaredScaledDistancebetweenVectors(Program.PointPosition[LocalToProcessIndex],
                    DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex], 
                    DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex]);
            }
            int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
            return getSquaredScaledDistancebetweenVectors(Program.PointPosition[LocalToProcessIndex], Solution.Y_k_i_[RealClusterIndex], Solution.Sigma_k_i_[RealClusterIndex]);

        }   // End getSquaredScaledDistancefromPoint

        public static double getSquaredScaledDistanceTweenActiveClusters(int ActiveClusterIndex1, int ActiveClusterIndex2, ClusteringSolution Solution)
        {
            double[] firstarg;
            double[] secondarg;
            double[] sigma;
            if (Solution.DistributedExecutionMode && (ActiveClusterIndex1 >= ClusteringSolution.NumberLocalActiveClusters))
            {
                int RemoteIndex1 = ActiveClusterIndex1 - ClusteringSolution.NumberLocalActiveClusters;
                firstarg = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex1];
                sigma = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedSigma_t_i[RemoteIndex1];
            }
            else
            {
                int RealClusterIndex1 = ClusteringSolution.RealClusterIndices[ActiveClusterIndex1];
                firstarg = Solution.Y_k_i_[RealClusterIndex1];
                sigma = Solution.Sigma_k_i_[RealClusterIndex1];
            }
            if (Solution.DistributedExecutionMode && (ActiveClusterIndex1 >= ClusteringSolution.NumberLocalActiveClusters))
            {
                int RemoteIndex2 = ActiveClusterIndex2 - ClusteringSolution.NumberLocalActiveClusters;
                secondarg = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedY_t_i[RemoteIndex2];
            }
            else
            {
                int RealClusterIndex2 = ClusteringSolution.RealClusterIndices[ActiveClusterIndex2];
                secondarg = Solution.Y_k_i_[RealClusterIndex2];
            }
            return getSquaredScaledDistancebetweenVectors(firstarg, secondarg, sigma);

        }   // End getSquaredScaledDistancefromCluster

        public static double getSquaredScaledDistancebetweenVectors(double[] vector1, double[] vector2, double[] sigma)
        {
            if (Program.SigmaMethod == 0)
            {
                double tmp4 = 0.0;
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    double tmp3 = (vector1[VectorIndex] - vector2[VectorIndex]);
                    tmp4 += tmp3 * tmp3;
                }
                return tmp4;
            }
            
            double tmp0 = (vector1[0] - vector2[0]);
            double tmp1 = (vector1[1] - vector2[1]);
            double tmp = (tmp0 * tmp0 / sigma[0]) + (tmp1 * tmp1 / sigma[1]);
            if (Program.ParameterVectorDimension == 2)
                return tmp;
            if (Program.ParameterVectorDimension == 3)
            {
                double tmp2 = (vector1[2] - vector2[2]);
                return tmp + tmp2 * tmp2 / sigma[2];
            }
            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                double tmp3 = (vector1[VectorIndex] - vector2[VectorIndex]);
                tmp += tmp3 * tmp3 / sigma[VectorIndex];
            }

            return tmp;

        }   // End getSquaredScaledDistancefromVector

        public static double getNOTSquaredScaledDistancebetweenVectors(double[] vector1, double[] vector2, double[] sigma)
        {
            if (Program.SigmaMethod == 0)
            {
                double tmp4 = 0.0;
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    double tmp3 = (vector1[VectorIndex] - vector2[VectorIndex]);
                    tmp4 += tmp3 * tmp3;
                }
                return Math.Sqrt(tmp4);
            }

            double tmp0 = (vector1[0] - vector2[0]);
            double tmp1 = (vector1[1] - vector2[1]);
            double tmp = (tmp0 * tmp0 / sigma[0]) + (tmp1 * tmp1 / sigma[1]);
            if (Program.ParameterVectorDimension == 2)
                return Math.Sqrt(tmp);
            if (Program.ParameterVectorDimension == 3)
            {
                double tmp2 = (vector1[2] - vector2[2]);
                return Math.Sqrt(tmp + tmp2 * tmp2 / sigma[2]);
            }
            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                double tmp3 = (vector1[VectorIndex] - vector2[VectorIndex]);
                tmp += tmp3 * tmp3 / sigma[VectorIndex];
            }

            return Math.Sqrt(tmp);

        }   // End getNOTSquaredScaledDistancefromVector

        public static double getNOTSquaredUNScaledDistancebetweenVectors(double[] vector1, double[] vector2)
        {
            double tmp0 = (vector1[0] - vector2[0]);
            double tmp1 = (vector1[1] - vector2[1]);
            double tmp = (tmp0 * tmp0) + (tmp1 * tmp1);
            if (Program.ParameterVectorDimension == 2)
                return Math.Sqrt(tmp);
            if (Program.ParameterVectorDimension == 3)
            {
                double tmp2 = (vector1[2] - vector2[2]);
                return Math.Sqrt(tmp + tmp2 * tmp2);
            }
            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                double tmp3 = (vector1[VectorIndex] - vector2[VectorIndex]);
                tmp += tmp3 * tmp3;
            }

            return Math.Sqrt(tmp);

        }   // End getSquaredScaledDistancefromVector

    }   // End class DAVectorParallelism
}
