using System;
using MPI;
using MathNet.Numerics.LinearAlgebra;
using SALSALibrary;
using LMatrix = MathNet.Numerics.LinearAlgebra.Matrix;
using LU = MathNet.Numerics.LinearAlgebra.LUDecomposition;

namespace Manxcat
{
    //  Routines for Generic Manxcat with parallel computation of components for small number of parameters
    public class GenericManxcat
    {
        public static ChisqFirstandSecond[] LocalAccum;
        public static ChisqFirstandSecond SummedoverThreadsAccum, SummedoverProcessesAccum;

        public static void SetupGenericManxcat()
        {
            LocalAccum = new ChisqFirstandSecond[SALSAUtility.ThreadCount];
            int NumberParms = Hotsun.npar;
            for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++)
            {
                LocalAccum[ThreadNo] = new ChisqFirstandSecond(NumberParms);
            }
            SummedoverThreadsAccum = new ChisqFirstandSecond(NumberParms);
            SummedoverProcessesAccum = new ChisqFirstandSecond(NumberParms);

            Hotsun.Number_VectorParameters = Hotsun.npar;
            Hotsun.ParameterVectorDimension = 1;
            Hotsun.DecomposeParameters = false;
            SALSAUtility.sequentialBLAS = true;
            Hotsun.fullmatrixset = true;
        }

        // End SetupGenericManxcat()

        public static void Accum(int ThreadNo, double value, double[] deriv)
        {
            // Well Loved Function accum for threading
            // Accum called for each point in Chisq Sum with normation Chisq =  sum(points) value**2
            // deriv[i] is d(value)/d(parameter i)

            LocalAccum[ThreadNo].chisq += value*value;
            int linearcount = 0;
            for (int iparm1 = 0; iparm1 < Hotsun.npar; iparm1++)
            {
                LocalAccum[ThreadNo].first[iparm1] += value*deriv[iparm1];
                for (int iparm2 = iparm1; iparm2 < Hotsun.npar; iparm2++)
                {
                    LocalAccum[ThreadNo].second[linearcount] += deriv[iparm1]*deriv[iparm2];
                    ++linearcount;
                }
            }
            return;
        }

        // End Accum(int ThreadNo, double value, double[] deriv)

        public static void AccumDebug(int index)
        {
            // SALSAUtility.SALSAPrint(0, index.ToString() + " Last two second derivs " + LocalAccum[0].second[26].ToString("e4") + " " + LocalAccum[0].second[27].ToString("e4"));
        }


        public static void ReInitializeAccum()
        {
            for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++)
            {
                ChisqFirstandSecond.Zeromember(ref LocalAccum[ThreadNo]);
            }
            ChisqFirstandSecond.Zeromember(ref SummedoverThreadsAccum);
            ChisqFirstandSecond.Zeromember(ref SummedoverProcessesAccum);
        }

        // End ReInitializeAccum()

        public static void AddupChisqContributions(Desertwind TotalSolution)
        {
            // Sum over Threads

            int NumberParms = SummedoverThreadsAccum.VectorSize;
            int SecondSize = SummedoverThreadsAccum.SecondSize;
            for (int ThreadNo = 0; ThreadNo < SALSAUtility.ThreadCount; ThreadNo++)
            {
                SummedoverThreadsAccum.chisq += LocalAccum[ThreadNo].chisq;
                for (int iparm = 0; iparm < NumberParms; iparm++)
                {
                    SummedoverThreadsAccum.first[iparm] += LocalAccum[ThreadNo].first[iparm];
                }
                for (int iparm = 0; iparm < SecondSize; iparm++)
                {
                    SummedoverThreadsAccum.second[iparm] += LocalAccum[ThreadNo].second[iparm];
                }
            }

            //  Sum over MPI Processes
            ChisqFirstandSecond UsethisAccum = SummedoverThreadsAccum;
            if (SALSAUtility.MPI_Size > 1)
            {
                UsethisAccum = SummedoverProcessesAccum;
                SummedoverProcessesAccum.chisq = SALSAUtility.MPI_communicator.Allreduce(SummedoverThreadsAccum.chisq,
                                                                                         Operation<double>.Add);
                SummedoverProcessesAccum.first = SALSAUtility.MPI_communicator.Allreduce(SummedoverThreadsAccum.first,
                                                                                         Operation<double>.Add);
                SummedoverProcessesAccum.second = SALSAUtility.MPI_communicator.Allreduce(
                    SummedoverThreadsAccum.second, Operation<double>.Add);
            }

            //  Put into Hotsun Form
            TotalSolution.Chisquared = UsethisAccum.chisq;
            Hotsun.zerocr = UsethisAccum.chisq;
            int linearcount = 0;
            for (int iparm1 = 0; iparm1 < NumberParms; iparm1++)
            {
                TotalSolution.first[iparm1][0] = UsethisAccum.first[iparm1];
                for (int iparm2 = iparm1; iparm2 < NumberParms; iparm2++)
                {
                    double tmp = UsethisAccum.second[linearcount];
                    linearcount++;
                    TotalSolution.FullMatrix[iparm1, iparm2][0, 0] = tmp;
                    TotalSolution.FullMatrix[iparm2, iparm1][0, 0] = tmp;
                    TotalSolution.ExactFullMatrix[iparm1, iparm2][0, 0] = tmp;
                    TotalSolution.ExactFullMatrix[iparm2, iparm1][0, 0] = tmp;
                }
            }
            // SALSAUtility.SALSAPrint(0, "Second Deriv " + TotalSolution.FullMatrix[5, 6][0, 0].ToString("E4"));
        }

        // End AddupChisqContributions(ChisqFirstandSecond[] SubTotal, Desertwind TotalSolution)

        public static void FindQlimits(Desertwind Solution, ref double Qhigh, ref double Qlow, ref int ReasontoStop1,
                                       ref int ReasontoStop2)
        {
            if (Hotsun.FullSecondDerivative)
                FindTraceandNorm(Solution.ExactFullMatrix, ref Hotsun.ChisqMatrixTrace, ref Hotsun.ChisqMatrixNorm);
            else
                FindTraceandNorm(Solution.FullMatrix, ref Hotsun.ChisqMatrixTrace, ref Hotsun.ChisqMatrixNorm);

            var ConventionalMatrix = new double[Hotsun.npar,Hotsun.npar];

            // Set Up Matrix to find eigenvalues
            // Scale but do NOT add Q
            double[,][,] Matrix;
            if (Hotsun.FullSecondDerivative)
                Matrix = Solution.ExactFullMatrix;
            else
                Matrix = Solution.FullMatrix;

            //  Set actual matrix
            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++)
                {
                    double MatrixElement = Matrix[GlobalIndex1, GlobalIndex2][0, 0];
                    if (Hotsun.UseDiagonalScaling)
                        MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex2][0];
                    ConventionalMatrix[GlobalIndex1, GlobalIndex2] = MatrixElement;
                } // End GlobalIndex2
            } // End GlobalIndex1
            // Find Minimum and Maximum eigenvalue of ConventionalMatrix


            // Begin Added by smbeason 5/17/2009
            LMatrix lmatrix = LMatrix.Create(ConventionalMatrix);
            EigenvalueDecomposition eigenValueDecomp = lmatrix.EigenvalueDecomposition;

            double minEigenValue = double.MaxValue;
            double maxEigenValue = double.MinValue;

            //Assuming you want on the real part...
            for (int i = 0; i < eigenValueDecomp.RealEigenvalues.Length; i++)
            {
                minEigenValue = Math.Min(minEigenValue, eigenValueDecomp.RealEigenvalues[i]);
                maxEigenValue = Math.Max(maxEigenValue, eigenValueDecomp.RealEigenvalues[i]);
            }
            // End Added by smbeason 5/17/2009

            ReasontoStop1 = 1;
            ReasontoStop2 = 1;
            Qlow = minEigenValue;
            Qhigh = maxEigenValue;
            return;
        }

        // End FindQlimits(Desertwind Solution, ref double Qhigh, ref double Qlow,ref int ReasontoStop1, ref int ReasontoStop2)

        public static void FindTraceandNorm(double[,][,] Matrix, ref double Trace, ref double Norm)
        {
            Trace = 0.0;
            Norm = 0.0;
            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                if (Hotsun.FixedParameter[GlobalIndex1][0])
                    continue;
                double RowNorm = 0.0;
                for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++)
                {
                    if (Hotsun.FixedParameter[GlobalIndex2][0])
                        continue;
                    double MatrixElement = Matrix[GlobalIndex1, GlobalIndex2][0, 0];
                    if (Hotsun.UseDiagonalScaling)
                        MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex2][0];
                    if (GlobalIndex1 == GlobalIndex2)
                    {
                        if (Hotsun.AddMarquardtQDynamically)
                            MatrixElement += Hotsun.Q;
                        Trace += Math.Abs(MatrixElement);
                    }
                    RowNorm += Math.Abs(MatrixElement);
                } // End GlobalIndex2

                if (Norm < RowNorm)
                    Norm = RowNorm;
            } // End GlobalIndex1
            return;
        }

        // End FindTraceandNorm(double[,][,] Matrix, double[][] GlobalxVector, ref double Trace, ref double Norm)

        public static bool SolveMatrix(double[][] Answer, Desertwind Solution)
        {
            var ConventionalMatrix = new double[Hotsun.npar,Hotsun.npar];
            var ConventionalFirst = new double[Hotsun.npar,1];
            var ConventionalAnswer = new double[Hotsun.npar][];
            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                ConventionalAnswer[GlobalIndex1] = new double[1];
            }

            double[,][,] Matrix;
            if (Hotsun.FullSecondDerivative)
                Matrix = Solution.ExactFullMatrix;
            else
                Matrix = Solution.FullMatrix;

            //  Set actual matrix
            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                if (Hotsun.FixedParameter[GlobalIndex1][0])
                    ConventionalFirst[GlobalIndex1, 0] = 0.0;
                else
                    ConventionalFirst[GlobalIndex1, 0] = Solution.first[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex1][0];

                for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++)
                {
                    double MatrixElement = Matrix[GlobalIndex1, GlobalIndex2][0, 0];
                    if (Hotsun.UseDiagonalScaling)
                        MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex2][0];
                    if (GlobalIndex1 == GlobalIndex2)
                    {
                        if (Hotsun.AddMarquardtQDynamically)
                            MatrixElement += Hotsun.Q;
                    }
                    ConventionalMatrix[GlobalIndex1, GlobalIndex2] = MatrixElement;
                } // End GlobalIndex2
            } // End GlobalIndex1

            // Form ConventionalMatrix-1 * ConventionalFirst

            LMatrix cMatrix = LMatrix.Create(ConventionalMatrix);
            LMatrix RightHandSide = LMatrix.Create(ConventionalFirst);
            var Fred = new LU(cMatrix);
            LMatrix MatrixAnswer = Fred.Solve(RightHandSide);
            ConventionalAnswer = MatrixAnswer;


            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                Answer[GlobalIndex1][0] = ConventionalAnswer[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex1][0];
            }
            bool success = true;

            return success;
        }

        // end SolveMatrix(double[][] Answer, Desertwind Solution)

        public static void GlobalMatrixVectorProduct(double[][] DistributedVector, Desertwind Solution, bool useexact,
                                                     double[][] GlobalxVector, double[][] GlobalVectoronRight)
        {
            // DistributedVector = Matrix (calculated using GlobalxVector) . GlobalVectoronRight
            // In generic case, Matrix is fully calculated and does not use GlobalxVector

            double[,][,] Matrix;
            if (useexact)
                Matrix = Solution.ExactFullMatrix;
            else
                Matrix = Solution.FullMatrix;

            for (int GlobalIndex1 = 0; GlobalIndex1 < Hotsun.npar; GlobalIndex1++)
            {
                if (Hotsun.FixedParameter[GlobalIndex1][0])
                    DistributedVector[GlobalIndex1][0] = 0.0;
                else
                {
                    double tmp = 0.0;
                    for (int GlobalIndex2 = 0; GlobalIndex2 < Hotsun.npar; GlobalIndex2++)
                    {
                        if (Hotsun.FixedParameter[GlobalIndex2][0])
                            continue;
                        double MatrixElement = Matrix[GlobalIndex1, GlobalIndex2][0, 0];
                        if (Hotsun.UseDiagonalScaling)
                            MatrixElement *= Hotsun.sqdginv[GlobalIndex1][0]*Hotsun.sqdginv[GlobalIndex2][0];
                        if (Hotsun.AddMarquardtQDynamically && (GlobalIndex1 == GlobalIndex2))
                            MatrixElement += Hotsun.Q;
                        tmp += MatrixElement*GlobalVectoronRight[GlobalIndex2][0];
                    } // End GlobalIndex2
                    DistributedVector[GlobalIndex1][0] = tmp;
                } // End varied parameter
            } // End GlobalIndex1
        }

        // End Generic GlobalMatrixVectorProduct
    }

    // End GenericManxcat

    [Serializable]
    public class ChisqFirstandSecond
    {
        public int SecondSize;
        public int VectorSize;
        public double chisq;
        public double[] first;
        public double[] second;

        public ChisqFirstandSecond(int NumberParms)
        {
            VectorSize = NumberParms;
            first = new double[NumberParms];
            SecondSize = (NumberParms*(1 + NumberParms))/2;
            second = new double[SecondSize];
            chisq = 0.0;
            for (int iparm = 0; iparm < NumberParms; iparm++)
            {
                first[iparm] = 0.0;
            }
            for (int iparm = 0; iparm < SecondSize; iparm++)
            {
                second[iparm] = 0.0;
            }
        }

        public static void Zeromember(ref ChisqFirstandSecond Instance)
        {
            Instance.chisq = 0.0;
            for (int iparm = 0; iparm < Instance.VectorSize; iparm++)
            {
                Instance.first[iparm] = 0.0;
            }
            for (int iparm = 0; iparm < Instance.SecondSize; iparm++)
            {
                Instance.second[iparm] = 0.0;
            }
        }
    }
}