﻿//smbeason
#if USE_UINT16
using TDistance = System.UInt16;
#elif USE_INT16
using System;
using MPI;
using Manxcat;
using Salsa.Core;
using Salsa.Core.Blas;
using Environment = MPI.Environment;
using TDistance = System.Int16;

#else
using TDistance = System.Double;
#endif

//smbeason

namespace SALSALibrary
{
    public class SALSAParallelism
    {
        // Set Parallelism related Parameters

//  Distance related variables
        public static TDistance[][] buffer;

        public static void SetupParallelism(ref string[] args)
        {
            //  Set up MPI
            SALSAUtility.MPI_Environment = new Environment(ref args);
            SALSAUtility.MPI_communicator = Communicator.world; //initializing MPI world communicator
            SALSAUtility.MPI_Rank = SALSAUtility.MPI_communicator.Rank; // Rank of this process
            SALSAUtility.MPI_Size = SALSAUtility.MPI_communicator.Size; // Number of MPI Processes

            // Set up MPI
            SALSAUtility.MPIperNodeCount = SALSAUtility.MPI_Size/SALSAUtility.NodeCount;

            //smbeason
            ManxcatCentral.Configuration.MPIperNodeCount = SALSAUtility.MPIperNodeCount;
            //smbeason

            if ((SALSAUtility.MPIperNodeCount*SALSAUtility.NodeCount) != SALSAUtility.MPI_Size)
            {
                Exception e = SALSAUtility.SALSAError("Inconsistent MPI counts Nodes "
                                                      + SALSAUtility.NodeCount.ToString() + " Size " +
                                                      SALSAUtility.MPI_Size.ToString());

                throw (e);
            }

            SALSAUtility.ParallelPattern = "Machine:" + Environment.ProcessorName + " " +
                                           SALSAUtility.ThreadCount.ToString() + "x" +
                                           SALSAUtility.MPIperNodeCount.ToString() + "x" +
                                           SALSAUtility.NodeCount.ToString();
            if (SALSAUtility.MPI_Rank == 0)
            {
                SALSAUtility.SALSAPrint(0, " Distance Data Type: " + typeof (TDistance));
                SALSAUtility.SALSAPrint(0, SALSAUtility.ParallelPattern);
            }
        }

        public static void TearDownParallelism()
        {
            if (SALSAUtility.MPI_Environment != null)
            {
                SALSAUtility.MPI_Environment.Dispose();
                SALSAUtility.MPI_Environment = null;
            }
        }

        public static void SetParallelDecomposition()
        {
            //	First divide points among processes
            Range[] processRanges = RangePartitioner.Partition(SALSAUtility.PointCount_Global, SALSAUtility.MPI_Size);
            Range processRange = processRanges[SALSAUtility.MPI_Rank]; // The answer for this process

            SALSAUtility.PointStart_Process = processRange.StartIndex;
            SALSAUtility.PointCount_Process = processRange.Length;
            SALSAUtility.VariedPointStart_Process = SALSAUtility.PointStart_Process;
            SALSAUtility.VariedPointCount_Process = SALSAUtility.PointCount_Process;
            SALSAUtility.PointCount_Largest = int.MinValue;

            foreach (Range r in processRanges)
            {
                SALSAUtility.PointCount_Largest = Math.Max(r.Length, SALSAUtility.PointCount_Largest);
            }

            // We need points per process for all processes as used by load balancing algorithm and then to reorder distance matrix consistently
            SALSAUtility.PointsperProcess = new int[SALSAUtility.MPI_Size];
            SALSAUtility.PointsperThreadperProcess = new int[SALSAUtility.MPI_Size][];

            for (int i = 0; i < SALSAUtility.MPI_Size; i++)
            {
                SALSAUtility.PointsperProcess[i] = processRanges[i].Length;
                Range[] threadprocessRanges = RangePartitioner.Partition(processRanges[i], SALSAUtility.ThreadCount);
                SALSAUtility.PointsperThreadperProcess[i] = new int[SALSAUtility.ThreadCount];
                for (int j = 0; j < SALSAUtility.ThreadCount; j++)
                {
                    SALSAUtility.PointsperThreadperProcess[i][j] = threadprocessRanges[j].Length;
                }
            }

            //	Now divide points among threads for this process
            Range[] threadRanges = RangePartitioner.Partition(processRanges[SALSAUtility.MPI_Rank],
                                                              SALSAUtility.ThreadCount);
            SALSAUtility.PointsperThread = new int[SALSAUtility.ThreadCount];
            SALSAUtility.StartPointperThread = new int[SALSAUtility.ThreadCount];

            for (int j = 0; j < SALSAUtility.ThreadCount; j++)
            {
                SALSAUtility.PointsperThread[j] = threadRanges[j].Length;
                SALSAUtility.StartPointperThread[j] = threadRanges[j].StartIndex;
            }
        }

        //  Read Distance Data
        public static void ReadDataFromFile(string fname)
        {
            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
            {
                SALSAUtility.SALSAPrint(1, "Starting to read data: " +
                                           " Distance Cut " + SALSAUtility.DistanceCut.ToString("F3"));
            }
            double countremoveddistances = 0.0;
            double counttotaldistances = 0.0;

            SALSAUtility.StoredDistanceOption = Math.Max(2, SALSAUtility.StoredDistanceOption);
            if (SALSAUtility.PointCount_Global == SALSAUtility.NumberOriginalPoints)
                SALSAUtility.DiskDistanceOption = Math.Max(2, SALSAUtility.DiskDistanceOption);

            // Remove unsupported options
            if (SALSAUtility.DiskDistanceOption == 3)
                SALSAUtility.StoredDistanceOption = 3;
            if (SALSAUtility.StoredDistanceOption == 3)
                SALSAUtility.CalcFixedCrossFixed = false;

// Set sizes of matrix on disk
            int rowcount = SALSAUtility.PointCount_Global;
            int colcount = SALSAUtility.PointCount_Global;
            if (SALSAUtility.DiskDistanceOption == 1)
            {
                rowcount = SALSAUtility.NumberOriginalPoints;
                colcount = SALSAUtility.NumberOriginalPoints;
            }
            if (SALSAUtility.DiskDistanceOption == 3) rowcount = SALSAUtility.NumberVariedPoints;

            //          MatrixTextReader reader = new MatrixTextReader(rowcount,colcount);
            var reader = new MatrixBinaryReader(rowcount, colcount);


            bool Oneread = true;
            if (SALSAUtility.StoredDistanceOption != SALSAUtility.DiskDistanceOption)
                Oneread = false;
            if (SALSAUtility.Usedreordered)
                Oneread = false;

            if (Oneread)
            {
#if USE_UINT16 
            SALSAUtility.PointDistances = reader.ReadUInt16(fname, SALSAUtility.PointStart_Process, SALSAUtility.PointCount_Process);
#elif USE_INT16
                SALSAUtility.PointDistances = reader.ReadInt16(fname, SALSAUtility.PointStart_Process,
                                                               SALSAUtility.PointCount_Process);
#else
            SALSAUtility.PointDistances = reader.ReadDouble(fname, SALSAUtility.PointStart_Process, SALSAUtility.PointCount_Process);
#endif
                int numberofcolumns = SALSAUtility.PointCount_Global;
                for (int GlobalPointIndex = SALSAUtility.PointStart_Process;
                     GlobalPointIndex < SALSAUtility.PointStart_Process + SALSAUtility.PointCount_Process;
                     GlobalPointIndex++)
                {
                    int rowindex = GlobalPointIndex;
                    if (SALSAUtility.StoredDistanceOption == 2)
                    {
                        rowindex = rowindex - SALSAUtility.PointStart_Process;
                    }
                    if (SALSAUtility.StoredDistanceOption == 3)
                    {
                        int originalpoint = SALSAUtility.UsedPointtoOriginalPointMap[rowindex];
                        int variedpoint = SALSAUtility.OriginalPointDisposition[originalpoint] - SALSAUtility.SALSASHIFT;
                        if (variedpoint < 0)
                        {
                            Exception e =
                                SALSAUtility.SALSAError(" Illegal Distance Request Used Point " + rowindex.ToString() +
                                                        " Original " + originalpoint.ToString());
                            throw (e);
                        }
                        rowindex = variedpoint - SALSAUtility.VariedPointStart_Process;
                    }
                    for (int columnindex = 0; columnindex < numberofcolumns; columnindex++)
                    {
                        TDistance Temp = SALSAUtility.PointDistances[rowindex][columnindex];
                        counttotaldistances = counttotaldistances + 1.0;

                        if (SALSAUtility.DistanceCut > 0.0)
                        {
                            double distancevalue = (SALSAUtility.PointDistances[rowindex][columnindex]/
                                                    (TDistance.MaxValue*1.0));
                            if (distancevalue > SALSAUtility.DistanceCut)
                            {
                                SALSAUtility.PointDistances[rowindex][columnindex] = TDistance.MaxValue;
                                countremoveddistances = countremoveddistances + 1.0;
                            }
                        }
                    }
                }
            }
            else
            {
                int colsread = SALSAUtility.PointCount_Global;
                if (SALSAUtility.DiskDistanceOption == 1)
                    colsread = SALSAUtility.NumberOriginalPoints;

                if (SALSAUtility.StoredDistanceOption == 2)
                    SALSAUtility.PointDistances = new TDistance[SALSAUtility.PointCount_Process][];
                if (SALSAUtility.StoredDistanceOption == 3)
                    SALSAUtility.PointDistances = new TDistance[SALSAUtility.VariedPointCount_Process][];

                for (int GlobalPointIndex = SALSAUtility.PointStart_Process;
                     GlobalPointIndex < SALSAUtility.PointStart_Process + SALSAUtility.PointCount_Process;
                     GlobalPointIndex++)
                {
                    int rowtostore = GlobalPointIndex;
                    if (SALSAUtility.StoredDistanceOption == 2)
                        rowtostore = GlobalPointIndex - SALSAUtility.PointStart_Process;
                    if (SALSAUtility.StoredDistanceOption == 3)
                    {
                        int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                        rowtostore = SALSAUtility.OriginalPointDisposition[OriginalIndex] - SALSAUtility.SALSASHIFT;
                        if (rowtostore < 0)
                            continue;
                        rowtostore = rowtostore - SALSAUtility.VariedPointStart_Process;
                    }


                    int rowtoread = -1;
                    if (SALSAUtility.DiskDistanceOption == 1)
                        rowtoread = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                    if (SALSAUtility.DiskDistanceOption == 2)
                    {
                        rowtoread = SALSAUtility.ActualtoNaiveUsedOrder[GlobalPointIndex];
                    }
                    if (SALSAUtility.DiskDistanceOption == 3)
                    {
                        int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                        rowtoread = SALSAUtility.OriginalPointDisposition[OriginalIndex] - SALSAUtility.SALSASHIFT;
                        if (rowtoread < 0)
                            continue;
                        rowtoread = SALSAUtility.ActualtoNaiveUsedOrder[rowtoread];
                    }
                    SALSAUtility.PointDistances[rowtostore] = new TDistance[colcount];
#if USE_UINT16 
                    buffer = reader.ReadUInt16(fname, rowtoread, 1);
#elif USE_INT16
                    buffer = reader.ReadInt16(fname, rowtoread, 1);
#else
                    buffer = reader.ReadDouble(fname, rowtoread, 1);
#endif
                    // Store buffer in PointDistances
                    for (int colIndex = 0; colIndex < colsread; colIndex++)
                    {
                        int coltostore = colIndex;
                        if (SALSAUtility.DiskDistanceOption == 1)
                        {
                            coltostore = SALSAUtility.OriginalPointtoUsedPointMap[colIndex];
                            if (coltostore < 0) continue;
                        }
                        else if (SALSAUtility.DiskDistanceOption > 1)
                            coltostore = SALSAUtility.NaivetoActualUsedOrder[colIndex];
                        SALSAUtility.PointDistances[rowtostore][coltostore] = buffer[0][colIndex];
                        counttotaldistances = counttotaldistances + 1.0;

                        if (SALSAUtility.DistanceCut > 0.0)
                        {
                            double distancevalue = (SALSAUtility.PointDistances[rowtostore][coltostore]/
                                                    (TDistance.MaxValue*1.0));
                            if (distancevalue > SALSAUtility.DistanceCut)
                            {
                                SALSAUtility.PointDistances[rowtostore][coltostore] = TDistance.MaxValue;
                                countremoveddistances = countremoveddistances + 1.0;
                            }
                        }
                    }
                }
            }
            SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming);
            counttotaldistances = SALSAUtility.MPI_communicator.Allreduce(counttotaldistances, Operation<double>.Add);
            countremoveddistances = SALSAUtility.MPI_communicator.Allreduce(countremoveddistances, Operation<double>.Add);
            SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming);
            double fractionleft = 1.0 - countremoveddistances/counttotaldistances;
            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
            {
                SALSAUtility.SALSAPrint(1,
                                        "Total Distances " + counttotaldistances.ToString("F0") +
                                        " Distances Removed on Input " + countremoveddistances.ToString("F0") +
                                        " Fraction Left " + fractionleft.ToString("F5"));
            }
        }

        // End ReadDataFromFile(string fname)

        // Input row and col are Used values
        // Return -1 if distance undefined
        public static double getDistanceValue(int row, int col)
        {
            if (SALSAUtility.StoredDistanceOption == 2)
            {
                row = row - SALSAUtility.PointStart_Process;
            }
            if (SALSAUtility.StoredDistanceOption == 3)
            {
                int originalpoint = SALSAUtility.UsedPointtoOriginalPointMap[row];
                int variedpoint = SALSAUtility.OriginalPointDisposition[originalpoint] - SALSAUtility.SALSASHIFT;
                if (variedpoint < 0)
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Illegal Distance Request Used Point " + row.ToString() + " Original " +
                                                originalpoint.ToString());
                    throw (e);
                }
                row = variedpoint - SALSAUtility.VariedPointStart_Process;
            }

            TDistance Temp = SALSAUtility.PointDistances[row][col];

            // Return -1.0 if input distance < 0 or if equals TDistance.MaxValue
            if (Temp == TDistance.MaxValue)
            {
                if (SALSAUtility.DistanceCut > 0.0)
                {
                    if (SALSAUtility.UndefinedDistanceValue < 0.0)
                        return -1.0;
                    if (SALSAUtility.DistanceProcessingOption == 2)
                        return SALSAUtility.UndefinedDistanceValue*SALSAUtility.UndefinedDistanceValue;
                    else
                        return SALSAUtility.UndefinedDistanceValue;
                }
                else
                    return 1.0;
            }


#if USE_UINT16

#endif
#if USE_INT16
            if (Temp < 0)
            {
                if (SALSAUtility.UndefinedDistanceValue < 0.0)
                    return -1.0;
                if (SALSAUtility.DistanceProcessingOption == 2)
                    return SALSAUtility.UndefinedDistanceValue*SALSAUtility.UndefinedDistanceValue;
                else
                    return SALSAUtility.UndefinedDistanceValue;
            }
#endif

#if USE_UINT16 || USE_INT16
            double distancevalue = (Temp/(TDistance.MaxValue*1.0));
#else
            if( Temp < 0 )
                return -1.0;
            distancevalue =  Temp ;
#endif

            if (SALSAUtility.DistanceProcessingOption == 2)
                distancevalue = distancevalue*distancevalue;
            return distancevalue;
        }

        // End getDistanceValue

        // Input row and col are Used values
        public static void putDistanceValue(int row, int col, double value)
        {
            if (SALSAUtility.StoredDistanceOption == 2)
            {
                row = row - SALSAUtility.PointStart_Process;
            }
            if (SALSAUtility.StoredDistanceOption == 3)
            {
                int originalpoint = SALSAUtility.UsedPointtoOriginalPointMap[row];
                int variedpoint = SALSAUtility.OriginalPointDisposition[originalpoint] - SALSAUtility.SALSASHIFT;
                if (variedpoint < 0)
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Illegal Distance Put Request Used Point " + row.ToString() +
                                                " Original " + originalpoint.ToString());
                    throw (e);
                }
                row = variedpoint - SALSAUtility.VariedPointStart_Process;
            }

            TDistance Temp;

#if USE_UINT16 || USE_INT16
            if (value > 1.0)
            {
                Exception e =
                    SALSAUtility.SALSAError(" Illegal Distance value Put Request Used Point " + value.ToString("F4") +
                                            " Coordinates " + row.ToString() + " " + col.ToString());
                throw (e);
            }
            if ((value == 1.0) && (SALSAUtility.DistanceCut > 0.0) && (SALSAUtility.UndefinedDistanceValue < 0.0))
            {
                // Inconsistent value
                Exception e =
                    SALSAUtility.SALSAError(" Illegal Distance value 1.0 Put Request Used Point " + value.ToString("F4") +
                                            " Coordinates " + row.ToString() + " " + col.ToString());
                throw (e);
            }
            if (value < 0.0)
            {
                if (SALSAUtility.DistanceCut < 0.0)
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Illegal Distance value < 0.0 Put Request Used Point " +
                                                value.ToString("F4") + " Coordinates " + row.ToString() + " " +
                                                col.ToString());
                    throw (e);
                }

                Temp = TDistance.MaxValue;
            }
            else
            {
#if USE_INT16
                Temp = Convert.ToInt16(value*TDistance.MaxValue);
#endif
#if USE_UINT16
                Temp =  Convert.ToUInt16(value * TDistance.MaxValue);
#endif
            }
#else
            Temp = value;
#endif
            SALSAUtility.PointDistances[row][col] = Temp;
            return;
        }

        // End putDistanceValue
    }

    // end SALSAParallelism
}

// End namespace SALSALibrary