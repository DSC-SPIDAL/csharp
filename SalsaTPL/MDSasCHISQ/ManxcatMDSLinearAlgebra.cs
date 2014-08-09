using System;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using Manxcat;
using SALSALibrary;
using UserRoutines = MDS.ManxcatMDS;

namespace MDS
{
    public class MDSLinearAlgebra
    {
        public static double[][] DistributedNewIteratedVector = null; // LHS of Iteration Equation
        public static double[][] DistributedOldIteratedVector; // RHS of Iteration Equation
        public static double[][] GlobalOldIteratedVector; // Global version of RHS of Iteration Equation
        public static double[][] DistributedCGVector_R; // Vector r in Conjugate Gradient Method
        public static double[][] DistributedCGVector_Q; // Vector q in Conjugate Gradient Method
        public static double[][] DistributedCGVector_P; // Vector p in Conjugate Gradient Method
        public static double[][] GlobalCGVector_P; // Global version of Vector p in Conjugate Gradient Method

        public static double[][] GlobalCGVector_R;
                                 // Global version of Vector R in Conjugate Gradient Method (only used for debug)

        public static void Initialize()
        {
            int LocalNumberofPoints = SALSAUtility.PointCount_Process;
            DistributedOldIteratedVector = new double[LocalNumberofPoints][];
            DistributedNewIteratedVector = new double[LocalNumberofPoints][];
            DistributedCGVector_R = new double[LocalNumberofPoints][];

            int GlobalNumberofPoints = SALSAUtility.PointCount_Global;
            GlobalOldIteratedVector = new double[GlobalNumberofPoints][];
            GlobalCGVector_R = new double[GlobalNumberofPoints][];

            for (int LongIndex = 0; LongIndex < LocalNumberofPoints; LongIndex++)
            {
                DistributedOldIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
                DistributedNewIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
                DistributedCGVector_R[LongIndex] = new double[Hotsun.ParameterVectorDimension];
            }
            for (int LongIndex = 0; LongIndex < GlobalNumberofPoints; LongIndex++)
            {
                GlobalOldIteratedVector[LongIndex] = new double[Hotsun.ParameterVectorDimension];
                GlobalCGVector_R[LongIndex] = new double[Hotsun.ParameterVectorDimension];
            }

            // re-use Initialized arrays for Conjugate Gradient
            DistributedCGVector_Q = DistributedNewIteratedVector;
            DistributedCGVector_P = DistributedOldIteratedVector;
            GlobalCGVector_P = GlobalOldIteratedVector;
        }

        // End Initialize() in MDSLinearAlgebra

        //  MaxIndicator = 0 Find Maximum Eigenvalue of Matrix
        //  MaxIndicator = 1 Find Maximum Eigenvalue of (Maximizr - Matrix)
        public static int PowerIterate(Desertwind Solution, int MaxIndicator, double Maximizer,
                                       out double PowerEigenvalue)
        {
            if (DistributedNewIteratedVector == null)
                Initialize();
            Hotsun.UseDiagonalScaling = true;
            Hotsun.AddMarquardtQDynamically = false;

            //  Set up Initial Power Vectors
            SetInitialPowerVector(DistributedNewIteratedVector);
            double AdotA = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector, DistributedNewIteratedVector);

            int somethingtodo = 0;
            int PowerIterationCount = 0;
            PowerEigenvalue = -1.0;
            while (true)
            {
                // Iterate over Power Multiplications by Chisq Matrix

                //  Normalize Current A and move from New to Old

                double OldNorm = 1.0/Math.Sqrt(AdotA);
                SALSABLAS.LinearCombineVector(DistributedOldIteratedVector, OldNorm,
                                              DistributedNewIteratedVector, 0.0, DistributedNewIteratedVector);

                //  Make a Global Vector of DistributedOldIteratedVector
                ManxcatCentral.MakeVectorGlobal(DistributedOldIteratedVector, GlobalOldIteratedVector);

                //  Form Chisq Matrix Product with Old Vector
                if (Hotsun.FullSecondDerivative)
                    UserRoutines.GlobalMatrixVectorProduct(DistributedNewIteratedVector, Solution, true,
                                                           Hotsun.GlobalParameter, GlobalOldIteratedVector);
                else
                    UserRoutines.GlobalMatrixVectorProduct(DistributedNewIteratedVector, Solution, false,
                                                           Hotsun.GlobalParameter, GlobalOldIteratedVector);

                //  Correct case MaxIndicator = 1
                if (MaxIndicator > 0)
                {
                    SALSABLAS.LinearCombineVector(DistributedNewIteratedVector, -1.0, DistributedNewIteratedVector,
                                                  Maximizer, DistributedOldIteratedVector);
                }
                SALSABLAS.LinearCombineVector(DistributedNewIteratedVector, 1.0, DistributedNewIteratedVector,
                                              Hotsun.addonforQcomputation, DistributedOldIteratedVector);

                //  Form Scalar Products
                double NewEigenvalue = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector,
                                                                     DistributedOldIteratedVector);
                AdotA = SALSABLAS.VectorScalarProduct(DistributedNewIteratedVector, DistributedNewIteratedVector);

                ++PowerIterationCount;
                somethingtodo = -1;
                if (SALSAUtility.MPI_Rank == 0)
                {
                    if (PowerIterationCount > 10 && (NewEigenvalue > 0.0))
                    {
                        // Arbitary criteria for starting +tests
                        somethingtodo = 0;
                        double scaleit = 1.0;
                        if (MaxIndicator > 0)
                        {
                            scaleit = Hotsun.extraprecision;
                        }
                        if (Math.Abs(NewEigenvalue - PowerEigenvalue) > PowerEigenvalue*scaleit*Hotsun.eigenvaluechange)
                            ++somethingtodo;
                        double delta = AdotA - NewEigenvalue*NewEigenvalue; // (Ax- Eigenvalue*Axold)**2
                        if (Math.Abs(delta) > NewEigenvalue*NewEigenvalue*scaleit*Hotsun.eigenvectorchange)
                            ++somethingtodo;
                    }
                }
                PowerEigenvalue = NewEigenvalue;

                if (SALSAUtility.MPI_Size > 1)
                {
                    SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
                    SALSAUtility.MPI_communicator.Broadcast(ref somethingtodo, 0);
                    SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
                }
                if (PowerIterationCount >= Hotsun.PowerIterationLimit)
                {
                    somethingtodo = -2;
                    break;
                }
                if (somethingtodo == 0)
                {
                    somethingtodo = PowerIterationCount;
                    break;
                }
            } // End while over PowerIterationCounts

            return somethingtodo;
        }

        // End PowerIterate()

        //  Solve for Distributed Answer = (Matrix-1) First 
        //  Note RealMatrixSize takes account of forced zero parameters which does not crop up in explicit
        //  algorithm as implemented in Matrix Vector Multiplier
        public static bool ConjugateGradientSolver(double[][] Answer, Desertwind Solution, bool useexact,
                                                   double[][] GlobalxVector, double[][] DistributedRHS,
                                                   ref int NumberofIterations, int RealMatrixSize, double LimitonNormofR)
        {
            bool matrixsuccess = true;

            //  Initialize
            SALSABLAS.zrword(Answer); // Zero Solution x called xshift in Manxcat and stored in Answer
            SALSABLAS.CopyVector(DistributedCGVector_R, DistributedRHS, 0, SALSAUtility.PointCount_Process);
                // Set R(0) = RHS
            double InitialRNorm = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);
            double RNormTest = InitialRNorm*LimitonNormofR; // Limit for Test on Iteration of Norm of R
            var SaveRNorm = new double[RealMatrixSize + 1];
            SaveRNorm[0] = InitialRNorm;
            int CountSteps = 0;

            double lastrho = 1.0;
            double currentrho = 1.0;

            //  Loop over Conjugate Gradient Steps
            while (true)
            {
                // Set value of rho
                lastrho = currentrho;
                currentrho = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);

                // Set Vector P
                ++CountSteps;
                if (CountSteps == 1)
                {
                    SALSABLAS.CopyVector(DistributedCGVector_P, DistributedCGVector_R, 0,
                                         SALSAUtility.PointCount_Process);
                }
                else
                {
                    SALSABLAS.LinearCombineVector(DistributedCGVector_P, currentrho/lastrho, DistributedCGVector_P,
                                                  1.0, DistributedCGVector_R);
                }

                //  Make a Global Vector of DistributedCGVector_P
                ManxcatCentral.MakeVectorGlobal(DistributedCGVector_P, GlobalCGVector_P);

                //  Distributed Q = Matrix . Global P
                UserRoutines.GlobalMatrixVectorProduct(DistributedCGVector_Q, Solution, useexact,
                                                       Hotsun.GlobalParameter, GlobalCGVector_P);

                //  New Answer is Old answer + (Current Rho / (P dot Q)) Vector P
                double PdotQ = SALSABLAS.VectorScalarProduct(DistributedCGVector_P, DistributedCGVector_Q);
                double alpha = currentrho/PdotQ;
                SALSABLAS.LinearCombineVector(Answer, alpha, DistributedCGVector_P, 1.0, Answer);

                // New residual R = Old Residual - (Current Rho / (P dot Q)) Vector Q
                SALSABLAS.LinearCombineVector(DistributedCGVector_R, -alpha, DistributedCGVector_Q, 1.0,
                                              DistributedCGVector_R);

                //  See if we can or should End
                double CurrentRNorm = SALSABLAS.VectorScalarProduct(DistributedCGVector_R, DistributedCGVector_R);
                SaveRNorm[CountSteps] = CurrentRNorm;
                bool TestRNorm = CurrentRNorm <= RNormTest;
                SALSAUtility.SynchronizeMPIvariable(ref TestRNorm);
                if (TestRNorm)
                {
                    matrixsuccess = true;
                    break;
                }

                // Singular Matrix
                if (CountSteps >= RealMatrixSize)
                {
                    matrixsuccess = false;
                    ManxcatCentral.MakeVectorGlobal(DistributedCGVector_R, GlobalCGVector_R);
                    if (SALSAUtility.MPI_Rank == 0)
                    {
                        SALSAUtility.SALSAPrint(0, " CG Failure after " + RealMatrixSize.ToString() + " Steps");
                        string ListofNorms = "";
                        int Normindex = CountSteps;
                        for (int inorm = 0; inorm < 10; inorm++)
                        {
                            ListofNorms += " " + SaveRNorm[Normindex].ToString("E4");
                            --Normindex;
                            if (Normindex < 0)
                                break;
                        }
                        SALSAUtility.SALSAPrint(0, "Last 10 Norms " + ListofNorms);

                        string fname = ManxcatCentral.ResultDirectoryName + "\\BadCGVector" +
                                       Hotsun.TotalCGFailures.ToString();
                        var sw = new StreamWriter(fname, false, Encoding.UTF8);
                        double fractionnorm = 1.0/Math.Sqrt(CurrentRNorm);
                        try
                        {
                            for (int GlobalPointIndex = 0;
                                 GlobalPointIndex < SALSAUtility.PointCount_Global;
                                 GlobalPointIndex++)
                            {
                                string Coordinates = "";
                                int UsedPointIndex = SALSAUtility.NaivetoActualUsedOrder[GlobalPointIndex];
                                double pointsize = 0.0;
                                for (int LocalVectorIndex = 0;
                                     LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                     LocalVectorIndex++)
                                {
                                    Coordinates += GlobalCGVector_R[UsedPointIndex][LocalVectorIndex].ToString("E4") +
                                                   "\t";
                                    pointsize += GlobalCGVector_R[UsedPointIndex][LocalVectorIndex]*
                                                 GlobalCGVector_R[UsedPointIndex][LocalVectorIndex];
                                }
                                pointsize = Math.Sqrt(pointsize)*fractionnorm;
                                sw.WriteLine(
                                    String.Format((GlobalPointIndex).ToString() + "\t" + Coordinates + " " +
                                                  pointsize.ToString("E4")));
                            }

                            sw.Flush();
                            sw.Close();
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine("Failed writing data in CG Solver " + e);
                            throw (e);
                        }
                    }
                    break;
                }
            } // End Loop over Conjugate Gradient Steps

            NumberofIterations = CountSteps;
            return matrixsuccess;
        }

        // End ConjugateGradientSolver(Desertwind Solution)

        // Scale Distributed Vector by a Global Vector
        public static void DiagScaleVector(double[][] ScaledDistributedVector, double[][] OldDistributedVector,
                                           double[][] GlobalScalingVector)
        {
            int LocalVectorDimension = ScaledDistributedVector[0].GetLength(0);

            // Parallel Setting
            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;
                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < LocalVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
                                         {
                                             ScaledDistributedVector[LongIndex][LocalVectorIndex] = 0.0;
                                         }
                                         else
                                         {
                                             ScaledDistributedVector[LongIndex][LocalVectorIndex] = OldDistributedVector
                                                                                                        [LongIndex][
                                                                                                            LocalVectorIndex
                                                                                                        ]
                                                                                                    *
                                                                                                    GlobalScalingVector[
                                                                                                        GlobalIndex][
                                                                                                            LocalVectorIndex
                                                                                                        ];
                                         }
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // // End DiagScaleVector    


        // Set Initial Power Vector to be a random number with normalization 1
        // Be sure to set fixed variables to be zero
        public static void SetInitialPowerVector(double[][] TobeSet)
        {
            int LocalVectorDimension = TobeSet[0].GetLength(0);
            var Randobject = new Random();

            // Parallel Setting
            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;
                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < LocalVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
                                         {
                                             TobeSet[LongIndex][LocalVectorIndex] = 0.0;
                                         }
                                         else
                                         {
                                             TobeSet[LongIndex][LocalVectorIndex] = Randobject.NextDouble();
                                         }
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // // End SetInitialPowerVector(double[][] TobeSet)    


        //  Calculate Matrix Global Vector product storing as a distributed vector
        public static void FindTraceandNorm(double[][,] MatrixDiagonals, double[][] GlobalxVector,
                                            ref double Trace, ref double Norm)
        {
            Trace = 0.0;
            Norm = 0.0;
            var FindTrace = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
            var FindNorm = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
//	Start Code to calculate Trace and Norm
                                 double LocalTrace = 0.0;
                                 double LocalNorm = 0.0;
                                 double WeightFunction1, WeightFunction2;

                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;
                                 for (int LongIndex1 = beginpoint; LongIndex1 < indexlen + beginpoint; LongIndex1++)
                                 {
                                     int GlobalIndex1 = LongIndex1 + SALSAUtility.PointStart_Process;
                                     for (int LocalVectorIndex1 = 0;
                                          LocalVectorIndex1 < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex1++)
                                     {
                                         if (Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1])
                                             continue;

                                         double RowNorm = 0.0;
                                         for (int GlobalIndex2 = 0;
                                              GlobalIndex2 < SALSAUtility.PointCount_Global;
                                              GlobalIndex2++)
                                         {
                                             if (
                                                 !ManxcatMDS.SetChisqWeights(GlobalIndex1, GlobalIndex2,
                                                                             out WeightFunction1, out WeightFunction2))
                                                 continue;

                                             double DistanceFudge1 = 1.0;
                                             double DistanceFudge2 = 1.0;
                                             double SquaredDistance = 0.0;
                                             double funcl = 0.0;
                                             if (GlobalIndex1 != GlobalIndex2)
                                             {
                                                 for (int LocalVectorIndex2 = 0;
                                                      LocalVectorIndex2 < Hotsun.ParameterVectorDimension;
                                                      LocalVectorIndex2++)
                                                 {
                                                     double tmp1 = GlobalxVector[GlobalIndex1][LocalVectorIndex2] -
                                                                   GlobalxVector[GlobalIndex2][LocalVectorIndex2];
                                                     SquaredDistance += tmp1*tmp1;
                                                 }
                                                 double ActualDistance = SquaredDistance;
                                                 if (UserRoutines.DistanceFormula == 1)
                                                 {
                                                     SquaredDistance += UserRoutines.SQUAREMinimumDistance;
                                                     ActualDistance = Math.Sqrt(SquaredDistance);
                                                     DistanceFudge1 = 0.5/ActualDistance;
                                                     DistanceFudge2 = DistanceFudge1/SquaredDistance;
                                                 }
                                                 funcl = WeightFunction1 - WeightFunction2*ActualDistance;
                                             }
                                             for (int LocalVectorIndex2 = 0;
                                                  LocalVectorIndex2 < Hotsun.ParameterVectorDimension;
                                                  LocalVectorIndex2++)
                                             {
                                                 if (Hotsun.FixedParameter[GlobalIndex2][LocalVectorIndex2])
                                                     continue;
                                                 double MatrixElement;
                                                 if (GlobalIndex1 == GlobalIndex2)
                                                 {
                                                     // Diagonal Term
                                                     MatrixElement =
                                                         MatrixDiagonals[LongIndex1][
                                                             LocalVectorIndex1, LocalVectorIndex2];
                                                     if (Hotsun.UseDiagonalScaling)
                                                         MatrixElement *=
                                                             Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1]*
                                                             Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex2];
                                                     if (Hotsun.AddMarquardtQDynamically &&
                                                         (LocalVectorIndex1 == LocalVectorIndex2))
                                                         MatrixElement += Hotsun.Q;
                                                     if (LocalVectorIndex1 == LocalVectorIndex2)
                                                         LocalTrace += Math.Abs(MatrixElement);
                                                 }
                                                 else
                                                 {
                                                     // Off Diagonal Term
                                                     double correction = 0.0;
                                                     double VectorCrossProduct =
                                                         (GlobalxVector[GlobalIndex1][LocalVectorIndex1] -
                                                          GlobalxVector[GlobalIndex2][LocalVectorIndex1])
                                                         *
                                                         (GlobalxVector[GlobalIndex1][LocalVectorIndex2] -
                                                          GlobalxVector[GlobalIndex2][LocalVectorIndex2]);
                                                     MatrixElement = -8.0*VectorCrossProduct*DistanceFudge1*
                                                                     DistanceFudge1*WeightFunction2*WeightFunction2;
                                                     if (Hotsun.FullSecondDerivative)
                                                     {
                                                         if ((UserRoutines.DistanceFormula == 2) &&
                                                             (LocalVectorIndex1 == LocalVectorIndex2))
                                                             correction = 4.0*funcl*WeightFunction2;
                                                         if ((UserRoutines.DistanceFormula == 1) &&
                                                             (LocalVectorIndex1 == LocalVectorIndex2))
                                                             correction = 4.0*funcl*WeightFunction2*DistanceFudge1;
                                                         if (UserRoutines.DistanceFormula == 1)
                                                             correction += - 4.0*funcl*VectorCrossProduct*
                                                                           WeightFunction2*DistanceFudge2;
                                                         MatrixElement += correction;
                                                     }
                                                     if (Hotsun.UseDiagonalScaling)
                                                         MatrixElement *=
                                                             Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1]*
                                                             Hotsun.sqdginv[GlobalIndex2][LocalVectorIndex2];
                                                 }
                                                 RowNorm += Math.Abs(MatrixElement);
                                             } // End Loop over LocalVectorIndex2
                                         } // End loop over GlobalIndex2

                                         if (LocalNorm < RowNorm)
                                             LocalNorm = RowNorm;
                                     } //  End Loop over LocalVectorIndex1
                                 } //  // End loop over LongIndex1

                                 FindTrace.addapoint(ThreadNo, LocalTrace);
                                 FindNorm.addapoint(ThreadNo, LocalNorm);
                             }); // End loop over Point dependent quantities

            FindTrace.sumoverthreadsandmpi();
            FindNorm.sumoverthreadsandmpi();

            Trace = FindTrace.Total;
            Norm = FindNorm.TotalMax;
            return;
        }

        // End FindTraceandNorm(double[][] MatrixDiagonals,double[][] GlobalxVector, out double Trace, out double Norm)
    }
}

// End Namespace MDS