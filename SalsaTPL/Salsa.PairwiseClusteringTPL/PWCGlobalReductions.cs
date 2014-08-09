using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MPI;
using Salsa.Core;

using SALSALibrary;
using Salsa.PairwiseClusteringTPL;

namespace SALSALibrary
{
    public class GlobalReductions
    {

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalOr = PWCUtility.MPI_communicator.Allreduce<bool>(TotalOr, Operation<bool>.LogicalOr);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                return;
            }
        }   // End FindBoolOr

        public class FindIntSum
        {
            private int[] NumberofPoints;
            private int NumberofThreads;
            private int[] Intvalue;
            public int TotalNumberofPoints;
            public int TotalInt;

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<int>(TotalNumberofPoints, Operation<int>.Add);
                    TotalInt = PWCUtility.MPI_communicator.Allreduce<int>(TotalInt, Operation<int>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                return;
            }
        }   // End FindIntSum
        
        public class FindDoubleArraySum
        {   // Used to do histograms
            // Must call startthread method at start of threads
            public double[] NumberofPoints;
            public int NumberofThreads;
            public int NumberinSum;
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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalSum = PWCUtility.MPI_communicator.Allreduce<double>(TotalSum, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                return;
            }
        }   // End FindDoubleArraySum

        public class FindDoubleMax
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[] Maxvalue;
            public double TotalNumberofPoints;
            public double TotalMax;

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalMax = PWCUtility.MPI_communicator.Allreduce<double>(TotalMax, Operation<double>.Max);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                return;
            }
        }   // End FindDoubleMax

        public class FindMeanSigma
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[] mean;
            public double[] square;
            public double TotalNumberofPoints;
            public double Totalmean;
            public double Totalsquare;
            public double Totalsigma;

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
                square[ThreadNo] += value1 * value1;
            }
            public void sumoverthreadsandmpi()
            {
                for (int ThreadNo = 0; ThreadNo < NumberofThreads; ThreadNo++)
                {
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    Totalmean += mean[ThreadNo];
                    Totalsquare += square[ThreadNo];
                }
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    Totalmean = PWCUtility.MPI_communicator.Allreduce<double>(Totalmean, Operation<double>.Add);
                    Totalsquare = PWCUtility.MPI_communicator.Allreduce<double>(Totalsquare, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;

                Totalmean = Totalmean / TotalNumberofPoints;
                Totalsquare = (Totalsquare / TotalNumberofPoints) - Totalmean * Totalmean;
                Totalsigma = Math.Sqrt(Math.Max(0.0, Totalsquare));
            }
        }   // End FindMeanSigma

        public class FindDoubleSum
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[] TotalinThread;
            public double TotalNumberofPoints;
            public double Total;

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

            public void zero()
            {
                this.TotalNumberofPoints = 0.0;
                this.Total = 0.0;
                for (int ThreadNo = 0; ThreadNo < this.NumberofThreads; ThreadNo++)
                {
                    this.NumberofPoints[ThreadNo] = 0.0;
                    this.TotalinThread[ThreadNo] = 0.0;
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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    Total = PWCUtility.MPI_communicator.Allreduce<double>(Total, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

            }
        }   // End FindDoubleSum

        public class FindDoubleMean
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[] mean;
            public double TotalNumberofPoints;
            public double Totalmean;

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    Totalmean = PWCUtility.MPI_communicator.Allreduce<double>(Totalmean, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;

                Totalmean = Totalmean / TotalNumberofPoints;
            }
        }   // End FindDoubleMean


        public class FindVectorIntSum
        {
            private int[] NumberofPoints;
            private int NumberofThreads;
            private int[][] VectorSum;
            public int TotalNumberofPoints;
            public int[] TotalVectorSum;
            private int ArraySize;

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<int>(TotalNumberofPoints, Operation<int>.Add);
                    TotalVectorSum = PWCUtility.MPI_communicator.Allreduce<int>(TotalVectorSum, Operation<int>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
            }

        }   // End FindVectorIntSum

        public class FindVectorDoubleMax
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[][] VectorMax;
            public double TotalNumberofPoints;
            public double[] TotalVectorMax;
            private int ArraySize;

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
                    VectorMax[ThreadNo][ArrayLoop] = Math.Max( VectorMax[ThreadNo][ArrayLoop], value[ArrayLoop]);
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
                        TotalVectorMax[ArrayLoop] = Math.Max( TotalVectorMax[ArrayLoop],  VectorMax[ThreadNo][ArrayLoop]);
                }

                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalVectorMax = PWCUtility.MPI_communicator.Allreduce<double>(TotalVectorMax, Operation<double>.Max);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

            }

        }   // End FindVectorDoubleMax

        public class FindVectorDoubleSum
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;
            private int ArraySize;

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
                        VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    TotalVectorSum[ArrayLoop] = 0.0;
            }

            public void zero()
            {
                this.TotalNumberofPoints = 0.0;
                for (int ThreadNo = 0; ThreadNo < this.NumberofThreads; ThreadNo++)
                {
                    this.NumberofPoints[ThreadNo] = 0.0;
                    for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                        this.VectorSum[ThreadNo][ArrayLoop] = 0.0;
                }
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    this.TotalVectorSum[ArrayLoop] = 0.0;
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

                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalVectorSum = PWCUtility.MPI_communicator.Allreduce<double>(TotalVectorSum, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                
            }

        }   // End FindVectorDoubleSum

        public class FindVectorDoubleSum2
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;
            private int ArraySize;
            private Range[] ParallelArrayRanges;

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

                ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, PWCUtility.ThreadCount);
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

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ArrayThreadNo) =>
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

                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalVectorSum = PWCUtility.MPI_communicator.Allreduce<double>(TotalVectorSum, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

            }

        }   // End FindVectorDoubleSum2

        public class FindVectorDoubleSum3
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;
            private int ArraySize;
            private int ArraySize1;
            private int ArraySize2;
            private Range[] ParallelArrayRanges;

            public FindVectorDoubleSum3(int NumThreads, int NumberinArray1, int NumberinArray2)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray1 * NumberinArray2;
                ArraySize1 = NumberinArray1;
                ArraySize2 = NumberinArray2;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];

                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;

                ParallelArrayRanges = RangePartitioner.Partition(ArraySize, PWCUtility.ThreadCount);
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
                int additive = index2 * ArraySize1;
                for (int ArrayLoop1 = 0; ArrayLoop1 < ArraySize1; ArrayLoop1++)
                    VectorSum[ThreadNo][ArrayLoop1 + additive] += value1[ArrayLoop1];
            }

            public void sumoverthreadsandmpi()
            {
                PWCUtility.StartSubTimer(PWCUtility.ThreadTiming);
                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ArrayThreadNo) =>
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
                PWCUtility.StopSubTimer(PWCUtility.ThreadTiming);

                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    int bigsize = TotalVectorSum.Length;
                    if (bigsize <= 4096)
                    {
                        TotalVectorSum = PWCUtility.MPI_communicator.Allreduce<double>(TotalVectorSum, Operation<double>.Add);
                    }
                    else
                    {
                        double[] buffer = new double[4096];
                        int start = 0;
                        while (start < bigsize)
                        {
                            int whatsleft = Math.Min(bigsize-start, 4096);
                            for (int innerloop = 0; innerloop < whatsleft; innerloop++)
                                buffer[innerloop] = TotalVectorSum[start + innerloop];
                            buffer = PWCUtility.MPI_communicator.Allreduce<double>(buffer, Operation<double>.Add);
                            for (int innerloop = 0; innerloop < whatsleft; innerloop++)
                                TotalVectorSum[start + innerloop] = buffer[innerloop];
                            start += whatsleft;
                        }

                    }
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
            }

        }   // End FindVectorDoubleSum3

        public class FindIndirectVectorDoubleSum
        {
            private double[] NumberofPoints;
            private int NumberofThreads;
            private double[][] VectorSum;
            public double TotalNumberofPoints;
            public double[] TotalVectorSum;
            private int ArraySize;
            private Range[] ParallelArrayRanges;

            public FindIndirectVectorDoubleSum(int NumThreads, int NumberinArray)
            {
                NumberofThreads = NumThreads;
                ArraySize = NumberinArray;

                NumberofPoints = new double[NumThreads];
                TotalVectorSum = new double[ArraySize];
                VectorSum = new double[NumThreads][];

                TotalNumberofPoints = 0.0;
                
                ParallelArrayRanges = RangePartitioner.Partition(NumberinArray, PWCUtility.ThreadCount);
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

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ArrayThreadNo) =>
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

                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    TotalVectorSum = PWCUtility.MPI_communicator.Allreduce<double>(TotalVectorSum, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

            }

        }   // End FindIndirectVectorDoubleSum

        public class FindArrayMean
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[][] mean;
            public double TotalNumberofPoints;
            public double[] Totalmean;
            public int ArraySize;

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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    Totalmean = PWCUtility.MPI_communicator.Allreduce<double>(Totalmean, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }

                if (TotalNumberofPoints < 0.5)
                    return;
                for (int ArrayLoop = 0; ArrayLoop < ArraySize; ArrayLoop++)
                    Totalmean[ArrayLoop]  = Totalmean[ArrayLoop] / TotalNumberofPoints;
            }
        }   // End FindArrayMean

        public class FindMinorMaxValuewithIndex
        {
            public double[] NumberofPoints;
            public int NumberofThreads;
            public double[] MaxOrMinvalue;
            public double TotalNumberofPoints;
            public double TotalMaxOrMin;
            public int TotalIndexValue;
            public int[] IndexValue;
            public int MinMaxPointer;   // =0 Min = 1 Max

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
                {   // Max
                    if ((IndexValue[ThreadNo] >= 0) && (MaxOrMinvalue[ThreadNo] > value))
                        return;
                }
                else
                {   // Min
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
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    if (MinMaxPointer != 0)
                        PWCUtility.AllReduceMaxWithIndex(ref TotalMaxOrMin, ref TotalIndexValue);
                    else
                        PWCUtility.AllReduceMinWithIndex(ref TotalMaxOrMin, ref TotalIndexValue);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                return;
            }
        }   // End FindMinorMaxValuewithIndex

        public class FindManyMinValuewithIndex
        {   // Finds top LimitNumberStored points by minimizing value given in addapoint
            // Rsults returned in order of figure of merit
            // Store values and indices
            // Uses FindMinimumSet
            // Results TotalNumberofPoints OrderedMinValue OrderedIndexValue

            public double[] NumberofPoints;
            public int NumberofThreads;
            public int Numbertofind;
            public double[][] MinValuebythread;
            public int[][] IndexValuebythread;
            public int[] CurrentWorstbythread;

            public double TotalNumberofPoints;
            public double[] TotalMinValue;
            public int[] TotalIndexValue;
            public int TotalWorst;
            public double[] OrderedMinValue;
            public int[] OrderedIndexValue;

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
                FindMinimumSet(value, indexposition, ref CurrentWorstbythread[ThreadNo], MinValuebythread[ThreadNo], IndexValuebythread[ThreadNo], Numbertofind);
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
                    TotalNumberofPoints += NumberofPoints[ThreadNo];
                    for (int storeloop = 0; storeloop < Numbertofind; storeloop++)
                    {
                        if (IndexValuebythread[ThreadNo][storeloop] < 0)
                            continue;   // End this thread
                        FindMinimumSet(MinValuebythread[ThreadNo][storeloop], IndexValuebythread[ThreadNo][storeloop], ref TotalWorst, TotalMinValue, TotalIndexValue, Numbertofind);
                    }
                }
                if (PWCUtility.MPI_Size > 1)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                    TotalNumberofPoints = PWCUtility.MPI_communicator.Allreduce<double>(TotalNumberofPoints, Operation<double>.Add);
                    PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                }
                // Sort in absolute order and accumulate over processes. This takes Numbertofindsteps
                for (int OrderLoop = 0; OrderLoop < Numbertofind; OrderLoop++)
                {
                    int localindex = -1;    // unset
                    double localvalue = -1.0;
                    int loopused = -1;
                    for (int internalloop = 0; internalloop < Numbertofind; internalloop++) 
                    {   // Find minimum
                        if (TotalIndexValue[internalloop] < 0)
                            continue;
                        if ((localindex < 0) || (TotalMinValue[internalloop] < localvalue) )
                        {
                            localindex = TotalIndexValue[internalloop];
                            localvalue = TotalMinValue[internalloop];
                            loopused = internalloop;
                        }
                    }
                    int oldlocalindex = localindex;
                    if (PWCUtility.MPI_Size > 1)
                    {
                        PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming1);
                        PWCUtility.AllReduceMinWithIndex(ref localvalue, ref localindex);
                        PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming1);
                    }
                    
                    OrderedMinValue[OrderLoop] = localvalue;
                    OrderedIndexValue[OrderLoop] = localindex;
                    if ((oldlocalindex >= 0) && (OrderedIndexValue[OrderLoop] == oldlocalindex))
                    {
                        TotalIndexValue[loopused] = -1;
                        TotalMinValue[loopused] = -1.0;
                    }
                }   // Loop over Order Loop

                return;
            }
        }   // End FindManyMinValuewithIndex

        //  Support finding list of minimum values by inserting new point with value newvalue and index newindex into lists SmallValues SmallIndices
        //  In SmallIndices negative values correspond to unset values
        //  NumberSmallOnes is total number wanted
        //  Currentcut is position in SmallValues, SmallIndices of largest min value
        public static void FindMinimumSet(double newvalue, int newindex, ref int currentcut, double[] SmallValues, int[] SmallIndices, int NumberSmallones)
        {
            if (currentcut < 0)
            {
                currentcut = 0;
                SmallValues[0] = newvalue;
                SmallIndices[0] = newindex;
                return;
            }
            if (SmallIndices[NumberSmallones - 1] < 0)
            {   // Not all positions are filled so add at next available
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
                if(SmallIndices[ivalue] < 0 )
                    continue;
                if (SmallValues[ivalue] > maxvalue)
                {
                    currentcut = ivalue;
                    maxvalue = SmallValues[ivalue];
                }
            }
            return;

        }   // End FindMinimumSet


    }   // End class GlobalReductions
}
