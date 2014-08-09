using System;
using System.Threading.Tasks;
using MPI;
using Salsa.Core;

namespace SALSALibrary
{
    public class GlobalReductions
    {
        public static void FindMinimumSet(double newvalue, int newindex, ref int currentcut, double[] SmallValues,
                                          int[] SmallIndices, int NumberSmallones)
        {
            if (currentcut < 0)
            {
                currentcut = 0;
                SmallValues[0] = newvalue;
                SmallIndices[0] = newindex;
                return;
            }
            if (SmallIndices[NumberSmallones - 1] < 0)
            {
                // Not all positions are filled so add at next available
                // Reset currentcut if worst
                for (int ivalue = 0; ivalue < NumberSmallones; ivalue++)
                {
                    if (SmallIndices[ivalue] < 0)
                    {
                        SmallValues[ivalue] = newvalue;
                        SmallIndices[ivalue] = newindex;
                        if (SmallValues[ivalue] > SmallValues[currentcut])
                        {
                            currentcut = ivalue;
                        }
                        return;
                    }
                }
            }
            if (newvalue >= SmallValues[currentcut])
                return;

            // Replace currentcut position with new values and Reset new worst position
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
        }

        #region Nested type: Find2DDoubleArraySum

        public class Find2DDoubleArraySum
        {
            public int InnerDimension;
            public double[] NumberOfPoints;
            public int NumberOfThreads;
            public int OuterDimension;
            public double[][][] Sum; // NumberOfThreads x OuterDimension x InnerDimension
            public double TotalNumberofPoints;
            public double[][] TotalSum; // OuterDimension x InnerDimension

            public Find2DDoubleArraySum(int numThreads, int outerDimension, int innerDimension)
            {
                NumberOfThreads = numThreads;
                OuterDimension = outerDimension;
                InnerDimension = innerDimension;

                NumberOfPoints = new double[numThreads];
                TotalSum = new double[outerDimension][];
                Sum = new double[numThreads][][];

                TotalNumberofPoints = 0.0;
                for (int i = 0; i < outerDimension; ++i)
                {
                    TotalSum[i] = new double[innerDimension];
                    for (int j = 0; j < innerDimension; ++j)
                    {
                        TotalSum[i][j] = 0.0;
                    }
                }
            }

            public void startthread(int threadNo)
            {
                NumberOfPoints[threadNo] = 0.0;
                Sum[threadNo] = new double[OuterDimension][];
                for (int i = 0; i < OuterDimension; ++i)
                {
                    Sum[threadNo][i] = new double[InnerDimension];
                    for (int j = 0; j < InnerDimension; ++j)
                    {
                        Sum[threadNo][i][j] = 0.0;
                    }
                }
            }

            public void addapoint(int threadNo, int row, int col)
            {
                if (row < 0 || row >= OuterDimension || col < 0 || col >= InnerDimension)
                    return;
                NumberOfPoints[threadNo] += 1.0;
                Sum[threadNo][row][col] += 1.0;
            }

            public void sumoverthreadsandmpi()
            {
                for (int threadNo = 0; threadNo < NumberOfThreads; threadNo++)
                {
                    TotalNumberofPoints += NumberOfPoints[threadNo];
                    for (int i = 0; i < OuterDimension; ++i)
                    {
                        for (int j = 0; j < InnerDimension; ++j)
                        {
                            TotalSum[i][j] += Sum[threadNo][i][j];
                        }
                    }
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    for (int i = 0; i < OuterDimension; ++i)
                    {
                        TotalSum[i] = SALSAUtility.MPI_communicator.Allreduce(TotalSum[i], Operation<double>.Add);
                    }
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        #region Nested type: FindArrayMean

        public class FindArrayMean
        {
            public int ArraySize;
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double TotalNumberofPoints;
            public double[] Totalmean;
            public double[][] mean;

            public FindArrayMean(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                Totalmean = new double[ArraySize];
                mean = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    mean[ThreadNo] = new double[ArraySize];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    {
                        mean[ThreadNo][ArrayLoop] = 0.0;
                    }
                }
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    Totalmean[ArrayLoop] = 0.0;
            }

            public void addapoint(int ThreadNo, double[] value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    mean[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                        Totalmean[ArrayLoop] += mean[ThreadNo][ArrayLoop];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    Totalmean = SALSAUtility.MPI_communicator.Allreduce(Totalmean, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    Totalmean[ArrayLoop] = Totalmean[ArrayLoop]/TotalNumberofPoints;
            }
        }

        #endregion

        #region Nested type: FindBoolOr

        public class FindBoolOr
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public bool[] Orvalue;
            public double TotalNumberofPoints;
            public bool TotalOr;

            public FindBoolOr(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                Orvalue = new bool[NumThreads];

                TotalNumberofPoints = 0.0;
                TotalOr = false;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    Orvalue[ThreadNo] = false;
                }
            }

            public void addapoint(int ThreadNo, bool value)
            {
                NumberofPoints[ThreadNo] += 1.0;
                Orvalue[ThreadNo] = Orvalue[ThreadNo] || value;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    TotalOr = Orvalue[ThreadNo] || TotalOr;
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalOr = SALSAUtility.MPI_communicator.Allreduce(TotalOr, Operation<bool>.LogicalOr);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        #region Nested type: FindCorrelation

        public class FindCorrelation
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double TotalNumberofPoints;
            public double Totalcross12;
            public double Totalmean1;
            public double Totalmean2;
            public double Totalsigma1;
            public double Totalsigma2;
            public double Totalsquare1;
            public double Totalsquare2;
            public double[] cross12;
            public double[] mean1;
            public double[] mean2;
            public double[] square1;
            public double[] square2;

            public FindCorrelation(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                mean1 = new double[NumThreads];
                mean2 = new double[NumThreads];
                square1 = new double[NumThreads];
                square2 = new double[NumThreads];
                cross12 = new double[NumThreads];

                TotalNumberofPoints = 0.0;
                Totalmean1 = 0.0;
                Totalmean2 = 0.0;
                Totalsquare1 = 0.0;
                Totalsquare2 = 0.0;
                Totalcross12 = 0.0;
                Totalsigma1 = 0.0;
                Totalsigma2 = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    mean1[ThreadNo] = 0.0;
                    mean2[ThreadNo] = 0.0;
                    square1[ThreadNo] = 0.0;
                    square2[ThreadNo] = 0.0;
                    cross12[ThreadNo] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double value1, double value2)
            {
                NumberofPoints[ThreadNo] += 1.0;
                mean1[ThreadNo] += value1;
                mean2[ThreadNo] += value2;
                square1[ThreadNo] += value1*value1;
                square2[ThreadNo] += value2*value2;
                cross12[ThreadNo] += value1*value2;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    Totalmean1 += mean1[ThreadNo];
                    Totalmean2 += mean2[ThreadNo];
                    Totalsquare1 += square1[ThreadNo];
                    Totalsquare2 += square2[ThreadNo];
                    Totalcross12 += cross12[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    Totalmean1 = SALSAUtility.MPI_communicator.Allreduce(Totalmean1, Operation<double>.Add);
                    Totalmean2 = SALSAUtility.MPI_communicator.Allreduce(Totalmean2, Operation<double>.Add);
                    Totalsquare1 = SALSAUtility.MPI_communicator.Allreduce(Totalsquare1, Operation<double>.Add);
                    Totalsquare2 = SALSAUtility.MPI_communicator.Allreduce(Totalsquare2, Operation<double>.Add);
                    Totalcross12 = SALSAUtility.MPI_communicator.Allreduce(Totalcross12, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;

                Totalmean1 = Totalmean1/TotalNumberofPoints;
                Totalmean2 = Totalmean2/TotalNumberofPoints;
                Totalsquare1 = (Totalsquare1/TotalNumberofPoints) - Totalmean1*Totalmean1;
                Totalsquare2 = (Totalsquare2/TotalNumberofPoints) - Totalmean2*Totalmean2;
                Totalcross12 = (Totalcross12/TotalNumberofPoints) - Totalmean1*Totalmean2;
                Totalsigma1 = Math.Sqrt(Totalsquare1);
                Totalsigma2 = Math.Sqrt(Totalsquare2);
                Totalcross12 = Totalcross12/(Totalsigma1*Totalsigma2);
            }

            public void print(string label, string FPformat)
            {
                if ((SALSAUtility.DebugPrintOption == 0) || (SALSAUtility.MPI_Rank != 0))
                    return;
                SALSAUtility.SALSAPrint(1,
                                        label + " means " + Totalmean1.ToString(FPformat) + " " +
                                        Totalmean2.ToString(FPformat)
                                        + " sigmas " + Totalsigma1.ToString(FPformat) + " " +
                                        Totalsigma2.ToString(FPformat) + " correl " + Totalcross12.ToString("F4"));
            }
        }

        #endregion

        // End FindBoolOr

        // End FindIntSum

        #region Nested type: FindDoubleArraySum

        public class FindDoubleArraySum
        {
            // Used to do histograms
            // Must call startthread method at start of threads
            public int NumberinSum;
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[][] Sum;
            public double TotalNumberofPoints;
            public double[] TotalSum;

            public FindDoubleArraySum(int NumThreads, int ArraySize)
            {
                NumberofThreads = NumThreads;
                NumberinSum = ArraySize;

                NumberofPoints = new double[NumThreads];
                TotalSum = new double[ArraySize];
                Sum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                for (int loop = 0; loop < ArraySize; loop++)
                {
                    TotalSum[loop] = 0.0;
                }
            }

            public void startthread(int ThreadNo)
            {
                NumberofPoints[ThreadNo] = 0.0;
                Sum[ThreadNo] = new double[NumberinSum];
                for (int loop = 0; loop < NumberinSum; loop++)
                {
                    Sum[ThreadNo][loop] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, int loopvalue)
            {
                if ((loopvalue < 0) || (loopvalue >= NumberinSum))
                    return;
                NumberofPoints[ThreadNo] += 1.0;
                Sum[ThreadNo][loopvalue] += 1.0;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int loop = 0; loop < NumberinSum; loop++)
                    {
                        TotalSum[loop] += Sum[ThreadNo][loop];
                    }
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalSum = SALSAUtility.MPI_communicator.Allreduce(TotalSum, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        // End FindDoubleArraySum

        #region Nested type: FindDoubleMax

        public class FindDoubleMax
        {
            public double[] Maxvalue;
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double TotalMax;
            public double TotalNumberofPoints;

            public FindDoubleMax(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                Maxvalue = new double[NumThreads];

                TotalNumberofPoints = 0.0;
                TotalMax = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    Maxvalue[ThreadNo] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double value)
            {
                NumberofPoints[ThreadNo] += 1.0;
                Maxvalue[ThreadNo] = Math.Max(Maxvalue[ThreadNo], value);
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    TotalMax = Math.Max(TotalMax, Maxvalue[ThreadNo]);
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalMax = SALSAUtility.MPI_communicator.Allreduce(TotalMax, Operation<double>.Max);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        // End FindDoubleMax

        // End FindDoubleSum

        #region Nested type: FindDoubleMean

        public class FindDoubleMean
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double TotalNumberofPoints;
            public double Totalmean;
            public double[] mean;

            public FindDoubleMean(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                mean = new double[NumThreads];

                TotalNumberofPoints = 0.0;
                Totalmean = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    mean[ThreadNo] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                mean[ThreadNo] += value1;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    Totalmean += mean[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    Totalmean = SALSAUtility.MPI_communicator.Allreduce(Totalmean, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;

                Totalmean = Totalmean/TotalNumberofPoints;
            }
        }

        #endregion

        #region Nested type: FindDoubleSum

        public class FindDoubleSum
        {
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly double[] TotalinThread;
            public double Total;
            public double TotalNumberofPoints;

            public FindDoubleSum(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                TotalinThread = new double[NumThreads];

                TotalNumberofPoints = 0.0;
                Total = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    TotalinThread[ThreadNo] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                TotalinThread[ThreadNo] += value1;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    Total += TotalinThread[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    Total = SALSAUtility.MPI_communicator.Allreduce(Total, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        #region Nested type: FindIndirectVectorDoubleSum

        public class FindIndirectVectorDoubleSum
        {
            private readonly int ArraySize;
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly Range[] ParallelArrayRanges;
            private readonly double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;

            public FindIndirectVectorDoubleSum(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];
                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;

                ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, SALSAUtility.ThreadCount);
            }

            public void startthread(int ThreadNo)
            {
                NumberofPoints[ThreadNo] = 0.0;
                VectorSum[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                {
                    VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, int NumberLocations, int[] location, double[] value1)
            {
                if (NumberLocations <= 0)
                    return;
                NumberofPoints[ThreadNo] += 1.0;
                for (int ArrayLoop = 0; ArrayLoop < NumberLocations; ArrayLoop++)
                    VectorSum[ThreadNo][location[ArrayLoop]] += value1[ArrayLoop];
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                    TotalNumberofPoints += NumberofPoints[ThreadNo];

                Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                             (ArrayThreadNo) =>
                                 {
                                     int beginindex = ParallelArrayRanges[ArrayThreadNo].StartIndex;
                                     int indexlength = ParallelArrayRanges[ArrayThreadNo].Length;
                                     for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++)
                                     {
                                         TotalVectorSum[ArrayLoop] = 0.0;
                                         for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                                             TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                                     }
                                 });

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalVectorSum = SALSAUtility.MPI_communicator.Allreduce(TotalVectorSum, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        #region Nested type: FindIntSum

        public class FindIntSum
        {
            private readonly int[] Intvalue;
            private readonly int[] NumberofPoints;
            private readonly int NumberofThreads;
            public int TotalInt;
            public int TotalNumberofPoints;

            public FindIntSum(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new int[NumThreads];
                Intvalue = new int[NumThreads];

                TotalNumberofPoints = 0;
                TotalInt = 0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0;
                    Intvalue[ThreadNo] = 0;
                }
            }

            public void addapoint(int ThreadNo, int value)
            {
                NumberofPoints[ThreadNo] += 1;
                Intvalue[ThreadNo] += value;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    TotalInt += Intvalue[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<int>.Add);
                    TotalInt = SALSAUtility.MPI_communicator.Allreduce(TotalInt, Operation<int>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        #region Nested type: FindManyMinValuewithIndex

        public class FindManyMinValuewithIndex
        {
            // Finds top LimitNumberStored points by minimizing value given in addapoint
            // Rsults returned in order of figure of merit
            // Store values and indices
            // Uses FindMinimumSet
            // Results TotalNumberofPoints OrderedMinValue OrderedIndexValue

            public int[] CurrentWorstbythread;
            public int[][] IndexValuebythread;
            public double[][] MinValuebythread;
            public double[] NumberofPoints;
            public int NumberofThreads;
            public int Numbertofind;
            public int[] OrderedIndexValue;
            public double[] OrderedMinValue;
            public int[] TotalIndexValue;
            public double[] TotalMinValue;
            public double TotalNumberofPoints;
            public int TotalWorst;

            public FindManyMinValuewithIndex(int NumThreads, int LimitNumberStored)
            {
                NumberofThreads = NumThreads;
                Numbertofind = LimitNumberStored;

                NumberofPoints = new double[NumThreads];
                MinValuebythread = new double[NumThreads][];
                IndexValuebythread = new int[NumThreads][];
                CurrentWorstbythread = new int[NumThreads];

                TotalNumberofPoints = 0.0;
                TotalMinValue = new double[LimitNumberStored];
                TotalIndexValue = new int[LimitNumberStored];
                OrderedMinValue = new double[LimitNumberStored];
                OrderedIndexValue = new int[LimitNumberStored];

                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    CurrentWorstbythread[ThreadNo] = -1;
                    MinValuebythread[ThreadNo] = new double[LimitNumberStored];
                    IndexValuebythread[ThreadNo] = new int[LimitNumberStored];
                    for (int storeloop = 0; storeloop < LimitNumberStored; storeloop++)
                    {
                        MinValuebythread[ThreadNo][storeloop] = -1.0;
                        IndexValuebythread[ThreadNo][storeloop] = -1;
                    }
                }
            }

            public void addapoint(int ThreadNo, int indexposition, double value)
            {
                NumberofPoints[ThreadNo] += 1.0;
                FindMinimumSet(value, indexposition, ref CurrentWorstbythread[ThreadNo], MinValuebythread[ThreadNo],
                               IndexValuebythread[ThreadNo], Numbertofind);
            }

            public void sumoverthreadsandmpi()
            {
                for (int storeloop = 0; storeloop < Numbertofind; storeloop++)
                {
                    TotalMinValue[storeloop] = -1.0;
                    TotalIndexValue[storeloop] = -1;
                }
                TotalWorst = -1;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    for (int storeloop = 0; storeloop < Numbertofind; storeloop++)
                    {
                        if (IndexValuebythread[ThreadNo][storeloop] < 0)
                            continue; // End this thread
                        FindMinimumSet(MinValuebythread[ThreadNo][storeloop], IndexValuebythread[ThreadNo][storeloop],
                                       ref TotalWorst, TotalMinValue, TotalIndexValue, Numbertofind);
                    }
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                // Sort in absolute order and accumulate over processes. This takes Numbertofindsteps
                for (int OrderLoop = 0; OrderLoop < Numbertofind; OrderLoop++)
                {
                    int localindex = -1; // unset
                    double localvalue = -1.0;
                    int loopused = -1;
                    for (int internalloop = 0; internalloop < Numbertofind; internalloop++)
                    {
                        // Find minimum
                        if (TotalIndexValue[internalloop] < 0)
                            continue;
                        if ((localindex < 0) || (TotalMinValue[internalloop] < localvalue))
                        {
                            localindex = TotalIndexValue[internalloop];
                            localvalue = TotalMinValue[internalloop];
                            loopused = internalloop;
                        }
                    }
                    int oldlocalindex = localindex;
                    if (SALSAUtility.MPI_Size > 1)
                    {
                        SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                        SALSAUtility.AllReduceMinWithIndex(ref localvalue, ref localindex);
                        SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                    }

                    OrderedMinValue[OrderLoop] = localvalue;
                    OrderedIndexValue[OrderLoop] = localindex;
                    if ((oldlocalindex >= 0) && (OrderedIndexValue[OrderLoop] == oldlocalindex))
                    {
                        TotalIndexValue[loopused] = -1;
                        TotalMinValue[loopused] = -1.0;
                    }
                } // Loop over Order Loop

                return;
            }
        }

        #endregion

        #region Nested type: FindMeanSigma

        public class FindMeanSigma
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double TotalNumberofPoints;
            public double Totalmean;
            public double Totalsigma;
            public double Totalsquare;
            public double[] mean;
            public double[] square;

            public FindMeanSigma(int NumThreads)
            {
                NumberofThreads = NumThreads;

                NumberofPoints = new double[NumThreads];
                mean = new double[NumThreads];
                square = new double[NumThreads];

                TotalNumberofPoints = 0.0;
                Totalmean = 0.0;
                Totalsquare = 0.0;
                Totalsigma = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    mean[ThreadNo] = 0.0;
                    square[ThreadNo] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                mean[ThreadNo] += value1;
                square[ThreadNo] += value1*value1;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    Totalmean += mean[ThreadNo];
                    Totalsquare += square[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    Totalmean = SALSAUtility.MPI_communicator.Allreduce(Totalmean, Operation<double>.Add);
                    Totalsquare = SALSAUtility.MPI_communicator.Allreduce(Totalsquare, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;

                Totalmean = Totalmean/TotalNumberofPoints;
                Totalsquare = (Totalsquare/TotalNumberofPoints) - Totalmean*Totalmean;
                Totalsigma = Math.Sqrt(Math.Max(0.0, Totalsquare));
            }
        }

        #endregion

        #region Nested type: FindMinorMaxValuewithIndex

        public class FindMinorMaxValuewithIndex
        {
            public int[] IndexValue;
            public double[] MaxOrMinvalue;
            public int MinMaxPointer; // =0 Min = 1 Max
            public double[] NumberofPoints;
            public int NumberofThreads;
            public int TotalIndexValue;
            public double TotalMaxOrMin;
            public double TotalNumberofPoints;

            public FindMinorMaxValuewithIndex(int NumThreads, int UseMax)
            {
                NumberofThreads = NumThreads;
                MinMaxPointer = UseMax;

                NumberofPoints = new double[NumThreads];
                MaxOrMinvalue = new double[NumThreads];
                IndexValue = new int[NumThreads];

                TotalNumberofPoints = 0.0;
                TotalMaxOrMin = 0.0;
                TotalIndexValue = -1;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    MaxOrMinvalue[ThreadNo] = 0.0;
                    IndexValue[ThreadNo] = -1;
                }
            }

            public void addapoint(int ThreadNo, int indexposition, double value)
            {
                NumberofPoints[ThreadNo] += 1.0;
                if (MinMaxPointer != 0)
                {
                    // Max
                    if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] > value))
                        return;
                }
                else
                {
                    // Min
                    if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] <= value))
                        return;
                }
                MaxOrMinvalue[ThreadNo] = value;
                IndexValue[ThreadNo] = indexposition;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    if (IndexValue[ThreadNo] < 0) continue;

                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    if (MinMaxPointer != 0)
                    {
                        if ((TotalIndexValue >= 0) && (TotalMaxOrMin > MaxOrMinvalue[ThreadNo]))
                            continue;
                    }
                    else
                    {
                        if ((TotalIndexValue >= 0) && (TotalMaxOrMin <= MaxOrMinvalue[ThreadNo]))
                            continue;
                    }

                    TotalMaxOrMin = MaxOrMinvalue[ThreadNo];
                    TotalIndexValue = IndexValue[ThreadNo];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    if (MinMaxPointer != 0)
                        SALSAUtility.AllReduceMaxWithIndex(ref TotalMaxOrMin, ref TotalIndexValue);
                    else
                        SALSAUtility.AllReduceMinWithIndex(ref TotalMaxOrMin, ref TotalIndexValue);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
                return;
            }
        }

        #endregion

        // End FindDoubleMean


        // End FindVectorIntSum

        #region Nested type: FindVectorDoubleMax

        public class FindVectorDoubleMax
        {
            private readonly int ArraySize;
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly double[][] VectorMax;
            public double TotalNumberofPoints;
            public double[] TotalVectorMax;

            public FindVectorDoubleMax(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                TotalVectorMax = new double[ArraySize];
                VectorMax = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    VectorMax[ThreadNo] = new double[ArraySize];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    {
                        VectorMax[ThreadNo][ArrayLoop] = -1.0E10;
                    }
                }
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    TotalVectorMax[ArrayLoop] = -1.0E10;
            }

            public void addapoint(int ThreadNo, double[] value)
            {
                NumberofPoints[ThreadNo] += 1.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    VectorMax[ThreadNo][ArrayLoop] = Math.Max(VectorMax[ThreadNo][ArrayLoop], value[ArrayLoop]);
            }

            public void addapoint(int ThreadNo, double value1, int position)
            {
                NumberofPoints[ThreadNo] += 1.0;
                VectorMax[ThreadNo][position] = Math.Max(VectorMax[ThreadNo][position], value1);
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                        TotalVectorMax[ArrayLoop] = Math.Max(TotalVectorMax[ArrayLoop], VectorMax[ThreadNo][ArrayLoop]);
                }

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalVectorMax = SALSAUtility.MPI_communicator.Allreduce(TotalVectorMax, Operation<double>.Max);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        // End FindVectorDoubleMax

        #region Nested type: FindVectorDoubleSum

        public class FindVectorDoubleSum
        {
            private readonly int ArraySize;
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;

            public FindVectorDoubleSum(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];
                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    NumberofPoints[ThreadNo] = 0.0;
                    VectorSum[ThreadNo] = new double[ArraySize];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    {
                        VectorSum[ThreadNo][ArrayLoop] = 0.0;
                    }
                }
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    TotalVectorSum[ArrayLoop] = 0.0;
            }

            public void addapoint(int ThreadNo, double[] value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }

            public void addapoint(int ThreadNo, double value1, int position)
            {
                NumberofPoints[ThreadNo] += 1.0;
                VectorSum[ThreadNo][position] += value1;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                        TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                }

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalVectorSum = SALSAUtility.MPI_communicator.Allreduce(TotalVectorSum, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        // End FindVectorDoubleSum

        #region Nested type: FindVectorDoubleSum2

        public class FindVectorDoubleSum2
        {
            private readonly int ArraySize;
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly Range[] ParallelArrayRanges;
            private readonly double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;

            public FindVectorDoubleSum2(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];
                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    TotalVectorSum[ArrayLoop] = 0.0;

                ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, SALSAUtility.ThreadCount);
            }

            public void startthread(int ThreadNo)
            {
                NumberofPoints[ThreadNo] = 0.0;
                VectorSum[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                {
                    VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
            }

            public void addapoint(int ThreadNo, double[] value1)
            {
                NumberofPoints[ThreadNo] += 1.0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }

            public void sumoverthreadsandmpi()
            {
                Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                             (ArrayThreadNo) =>
                                 {
                                     int beginindex = ParallelArrayRanges[ArrayThreadNo].StartIndex;
                                     int indexlength = ParallelArrayRanges[ArrayThreadNo].Length;
                                     for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++)
                                     {
                                         TotalVectorSum[ArrayLoop] = 0.0;
                                         for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                                             TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                                     }
                                 });

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    TotalVectorSum = SALSAUtility.MPI_communicator.Allreduce(TotalVectorSum, Operation<double>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        // End FindVectorDoubleSum2

        #region Nested type: FindVectorDoubleSum3

        public class FindVectorDoubleSum3
        {
            private readonly int ArraySize;
            private readonly int ArraySize1;
            private readonly double[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly Range[] ParallelArrayRanges;
            private readonly double[][] VectorSum;
            private int ArraySize2;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;

            public FindVectorDoubleSum3(int NumThreads, int NumberinArray1, int NumberinArray2)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray1*NumberinArray2;
                ArraySize1 = NumberinArray1;
                ArraySize2 = NumberinArray2;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];

                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;

                ParallelArrayRanges = RangePartitioner.Partition(ArraySize, SALSAUtility.ThreadCount);
            }

            public void startthread(int ThreadNo)
            {
                NumberofPoints[ThreadNo] = 0.0;
                VectorSum[ThreadNo] = new double[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    VectorSum[ThreadNo][ArrayLoop] = 0.0;
            }

            public void addapoint(int ThreadNo, double[] value1, int index2)
            {
                NumberofPoints[ThreadNo] += 1.0;
                int additive = index2*ArraySize1;
                for (int ArrayLoop1 = 0; ArrayLoop1 < ArraySize1; ArrayLoop1++)
                    VectorSum[ThreadNo][ArrayLoop1 + additive] += value1[ArrayLoop1];
            }

            public void sumoverthreadsandmpi()
            {
                SALSAUtility.StartSubTimer(SALSAUtility.ThreadTiming);
                Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                             (ArrayThreadNo) =>
                                 {
                                     int beginindex = ParallelArrayRanges[ArrayThreadNo].StartIndex;
                                     int indexlength = ParallelArrayRanges[ArrayThreadNo].Length;
                                     for (int ArrayLoop = beginindex; ArrayLoop < beginindex + indexlength; ArrayLoop++)
                                     {
                                         double tmp = 0.0;
                                         for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                                             tmp += VectorSum[ThreadNo][ArrayLoop];
                                         TotalVectorSum[ArrayLoop] = tmp;
                                     }
                                 });
                SALSAUtility.StopSubTimer(SALSAUtility.ThreadTiming);

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<double>.Add);
                    int bigsize = TotalVectorSum.Length;
                    if (bigsize <= 4096)
                    {
                        TotalVectorSum = SALSAUtility.MPI_communicator.Allreduce(TotalVectorSum, Operation<double>.Add);
                    }
                    else
                    {
                        var buffer = new double[4096];
                        int start = 0;
                        while (start < bigsize)
                        {
                            int whatsleft = Math.Min(bigsize - start, 4096);
                            for (int innerloop = 0; innerloop < whatsleft; innerloop++)
                                buffer[innerloop] = TotalVectorSum[start + innerloop];
                            buffer = SALSAUtility.MPI_communicator.Allreduce(buffer, Operation<double>.Add);
                            for (int innerloop = 0; innerloop < whatsleft; innerloop++)
                                TotalVectorSum[start + innerloop] = buffer[innerloop];
                            start += whatsleft;
                        }
                    }
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        #region Nested type: FindVectorIntSum

        public class FindVectorIntSum
        {
            private readonly int ArraySize;
            private readonly int[] NumberofPoints;
            private readonly int NumberofThreads;
            private readonly int[][] VectorSum;
            public int TotalNumberofPoints;
            public int[] TotalVectorSum;

            public FindVectorIntSum(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new int[NumThreads];
                TotalVectorSum = new int[ArraySize];
                VectorSum = new int[NumThreads][];

                TotalNumberofPoints = 0;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    TotalVectorSum[ArrayLoop] = 0;
            }

            public void startthread(int ThreadNo)
            {
                NumberofPoints[ThreadNo] = 0;
                VectorSum[ThreadNo] = new int[ArraySize];
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                {
                    VectorSum[ThreadNo][ArrayLoop] = 0;
                }
            }

            public void addapoint(int ThreadNo, int[] value1)
            {
                NumberofPoints[ThreadNo] += 1;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    VectorSum[ThreadNo][ArrayLoop] += value1[ArrayLoop];
            }

            public void addapoint(int ThreadNo, int value1, int position)
            {
                NumberofPoints[ThreadNo] += 1;
                VectorSum[ThreadNo][position] += value1;
            }

            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                        TotalVectorSum[ArrayLoop] += VectorSum[ThreadNo][ArrayLoop];
                }
                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = SALSAUtility.MPI_communicator.Allreduce(TotalNumberofPoints,
                                                                                  Operation<int>.Add);
                    TotalVectorSum = SALSAUtility.MPI_communicator.Allreduce(TotalVectorSum, Operation<int>.Add);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming1);
                }
            }
        }

        #endregion

        // End FindVectorDoubleSum3

        // End Find2DDoubleArraySum
    }

    // End class GlobalReductions
}