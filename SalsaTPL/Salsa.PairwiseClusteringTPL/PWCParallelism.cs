using System;
using System.IO;
using MPI;
using Salsa.Core.Blas;
using Salsa.Core;

#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using TDistance = System.Int16;
#else
using TDistance = System.Double;
#endif

namespace SALSALibrary
{
    public class PWCParallelism
    {
        public static void SetupParallelism(ref string[] args)
        {
            //  Set up MPI
            PWCUtility.MPI_Environment = new MPI.Environment(ref args);
            PWCUtility.MPI_communicator = Communicator.world;         //initializing MPI world communicator
            PWCUtility.MPI_Rank = PWCUtility.MPI_communicator.Rank; // Rank of this process
            PWCUtility.MPI_Size = PWCUtility.MPI_communicator.Size; // Number of MPI Processes

            // Set up MPI
            PWCUtility.MPIperNodeCount = PWCUtility.MPI_Size / PWCUtility.NodeCount;

            if ((PWCUtility.MPIperNodeCount * PWCUtility.NodeCount) != PWCUtility.MPI_Size)
            {
                Exception e = PWCUtility.SALSAError("Inconsistent MPI counts Nodes "
                    + PWCUtility.NodeCount.ToString() + " Size " + PWCUtility.MPI_Size.ToString());

                throw (e);
            }

            PWCUtility.ParallelPattern = "Machine:" + MPI.Environment.ProcessorName.ToString() + " " + PWCUtility.ThreadCount.ToString() + "x" + PWCUtility.MPIperNodeCount.ToString() + "x" + PWCUtility.NodeCount.ToString();
            if (PWCUtility.MPI_Rank == 0)
            {
                PWCUtility.SALSAPrint(0, " Distance Data Type: " + typeof(TDistance));
                PWCUtility.SALSAPrint(0, PWCUtility.ParallelPattern);
            }

        }   // End SetupParallelism

        public static void TearDownParallelism()
        {
            if (PWCUtility.MPI_Environment != null)
            {
                PWCUtility.MPI_Environment.Dispose();
                PWCUtility.MPI_Environment = null;
            }
        }   // End TearDownParallelism

        public static void SetParallelDecomposition()
        {
            //	First divide points among processes
            Range[] processRanges = RangePartitioner.Partition(PWCUtility.PointCount_Global, PWCUtility.MPI_Size);
            Range processRange = processRanges[PWCUtility.MPI_Rank];  // The answer for this process

            PWCUtility.PointStart_Process = processRange.StartIndex;
            PWCUtility.PointCount_Process = processRange.Length;
            PWCUtility.PointCount_Largest = int.MinValue;

            foreach (Range r in processRanges)
            {
                PWCUtility.PointCount_Largest = Math.Max(r.Length, PWCUtility.PointCount_Largest);
            }

            // We need points per process for all processes as used by load balancing algorithm and then to reorder distance matrix consistently
            PWCUtility.PointsperProcess = new int[PWCUtility.MPI_Size];
            PWCUtility.PointsperThreadperProcess = new int[PWCUtility.MPI_Size][];

            for (int i = 0; i < PWCUtility.MPI_Size; i++)
            {
                PWCUtility.PointsperProcess[i] = processRanges[i].Length;
                Range[] threadprocessRanges = RangePartitioner.Partition(processRanges[i], PWCUtility.ThreadCount);
                PWCUtility.PointsperThreadperProcess[i] = new int[PWCUtility.ThreadCount];
                for (int j = 0; j < PWCUtility.ThreadCount; j++)
                {
                    PWCUtility.PointsperThreadperProcess[i][j] = threadprocessRanges[j].Length;
                }
            }

            //	Now divide points among threads for this process
            Range[] threadRanges = RangePartitioner.Partition(processRanges[PWCUtility.MPI_Rank], PWCUtility.ThreadCount);
            PWCUtility.PointsperThread = new int[PWCUtility.ThreadCount];
            PWCUtility.StartPointperThread = new int[PWCUtility.ThreadCount];

            for (int j = 0; j < PWCUtility.ThreadCount; j++)
            {
                PWCUtility.PointsperThread[j] = threadRanges[j].Length;
                PWCUtility.StartPointperThread[j] = threadRanges[j].StartIndex;
            }
        }

        /// <summary>
        /// Read data from partial distance matrices. Here the assumption is that per node parallelism
        /// is same as for the program that generated these matrices. Also the total number of nodes
        /// are assumed to be both same in number and same in physically. 
        /// </summary>
        /// <param name="path">Local path where partial matrices are stored</param>
        /// <param name="prefix">Name prefix of partial matrix files</param>
        /// <param name="ext">File extension of partial matrix files</param>
        /// <param name="totalThreads">Number of threads per process times number of processes</param>
        public static void ReadDataFromFiles(string path, string prefix, string ext)
        {
            int totalThreads = PWCUtility.MPI_Size * PWCUtility.ThreadCount;

            // Divide the total number of points among total thread count
            Range[] threadRanges = RangePartitioner.Partition(PWCUtility.PointCount_Global, totalThreads);
            PWCUtility.StartPointperThread = new int[PWCUtility.ThreadCount];
            PWCUtility.PointsperThread = new int[PWCUtility.ThreadCount];

            for (int i = 0; i < PWCUtility.ThreadCount; i++)
            {
                PWCUtility.StartPointperThread[i] = threadRanges[i].StartIndex;
                PWCUtility.PointsperThread[i] = threadRanges[i].Length;
            }

            // Each process handles more than one thread. So set the process ranges accordingly
            Range[] processRanges = new Range[PWCUtility.MPI_Size];
            int start, end;

            // At each iteration i is set to a multiple of PWCUtility.ThreadCount, i.e. i % PWCUtility.ThreadCount is 0
            for (int i = 0; i < totalThreads; )
            {
                start = threadRanges[i].StartIndex;
                end = threadRanges[i + PWCUtility.ThreadCount - 1].EndIndex;
                processRanges[i / PWCUtility.ThreadCount] = new Range(start, end);
                i += PWCUtility.ThreadCount;
            }

            // The range of data points handled by this process
            Range processRange = processRanges[PWCUtility.MPI_Rank];
            PWCUtility.PointCount_Process = processRange.Length;
            PWCUtility.PointStart_Process = processRange.StartIndex;

            // Find the largest number of points handled by a process
            PWCUtility.PointCount_Largest = int.MinValue;
            foreach (Range r in processRanges)
            {
                PWCUtility.PointCount_Largest = Math.Max(r.Length, PWCUtility.PointCount_Largest);
            }

            // Create a reader to read only the number of rows handled by this process. Number of cols is
            // equal to the total number of datapoints.
            MatrixBinaryReader reader = new MatrixBinaryReader(processRange.Length, PWCUtility.PointCount_Global);

            start = PWCUtility.MPI_Rank * PWCUtility.ThreadCount;
            end = start + PWCUtility.ThreadCount - 1;

            // todo: saliya - fix Read methods for other types
#if USE_UINT16 
            PWCUtility.PointDistances = reader.ReadUInt16(fname, processRange.StartIndex, processRange.Length);
#elif USE_INT16
            PWCUtility.PointDistances = reader.ReadInt16(path, prefix, ext, start, end, threadRanges);
#else
            PWCUtility.PointDistances = reader.ReadDouble(fname, processRange.StartIndex, processRange.Length);
#endif

        }

        // read data from file to memory
        // set starting position and end number for data points assigned to each thread. 
        // These are in startindex and lenindex respectively
        // Return data stored as
        // Distance [ClusterCenter,ClusterIndex] with all ClusterIndex<=ClusterCenter stored in order for each ClusterCenter. 
        // Diagonal values are stored (as zero)
        // Total space m(m+1)/2 if lower triangular store (checkerboard =1)
        public static void ReadDataFromFile(string fname)
        {
            //	First divide points among processes
            Range[] processRanges = RangePartitioner.Partition(PWCUtility.PointCount_Global, PWCUtility.MPI_Size);
            Range processRange = processRanges[PWCUtility.MPI_Rank];

            PWCUtility.PointCount_Process = processRange.Length;
            PWCUtility.PointStart_Process = processRange.StartIndex;
            PWCUtility.PointCount_Largest = int.MinValue;

            foreach (Range r in processRanges)
            {
                PWCUtility.PointCount_Largest = Math.Max(r.Length, PWCUtility.PointCount_Largest);
            }

            //	Now divide points among threads
            Range[] threadRanges = RangePartitioner.Partition(processRange, PWCUtility.ThreadCount);
            PWCUtility.StartPointperThread = new int[PWCUtility.ThreadCount];
            PWCUtility.PointsperThread = new int[PWCUtility.ThreadCount];

            for (int i = 0; i < PWCUtility.ThreadCount; i++)
            {
                PWCUtility.StartPointperThread[i] = threadRanges[i].StartIndex;
                PWCUtility.PointsperThread[i] = threadRanges[i].Length;
            }

            MatrixBinaryReader reader = new MatrixBinaryReader(PWCUtility.PointCount_Global, PWCUtility.PointCount_Global);

#if USE_UINT16 
            PWCUtility.PointDistances = reader.ReadUInt16(fname, processRange.StartIndex, processRange.Length);
#elif USE_INT16
            PWCUtility.PointDistances = reader.ReadInt16(fname, processRange.StartIndex, processRange.Length);
#else
            PWCUtility.PointDistances = reader.ReadDouble(fname, processRange.StartIndex, processRange.Length);
#endif

        }   //  End routine controlling reading of data

        public static double getDistanceValue(int row, int col)
        {
            row -= PWCUtility.PointStart_Process;

            try
            {
#if USE_UINT16 || USE_INT16
                return (PWCUtility.PointDistances[row, col] / (TDistance.MaxValue * 1.0));
#else
                return PWCUtility.PointDistances[row, col];
#endif
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error reading distances[{0},{1}]", row, col);
                throw ex;
            }

        }

        public static int OwnerforThisPoint(int GlobalPointIndex)
        {   // Return process number for GlobalPointIndex

            int startpoint = 0;
            for (int mpiloop = 0; mpiloop < PWCUtility.MPI_Size; mpiloop++)
            {
                int endpoint = startpoint + PWCUtility.PointsperProcess[mpiloop];
                if (GlobalPointIndex < endpoint)
                    return mpiloop;
                startpoint = endpoint;
            }
            PWCUtility.SALSAError(" Illegal Point in No Process " + GlobalPointIndex.ToString());
            return -1;

        }   // End InThisProcess(int GlobalPointIndex)

    }
}
