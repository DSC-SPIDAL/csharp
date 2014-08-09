using System;
using SALSALibrary;

namespace Manxcat
{
    public class Hotsun
    {
        /*  Delegate Functions needed for each application */

        #region Delegates

        public delegate bool CalcfgSignature(Desertwind Solution);

        public delegate void FindQlimitsSignature(
            Desertwind Solution, ref double Qhigh, ref double Qlow, ref int ReasontoStop1, ref int ReasontoStop2);

        public delegate void GlobalMatrixVectorProductSignature(
            double[][] DistributedVector, Desertwind Solution, bool usexact, double[][] GlobalxVector,
            double[][] GlobalVectoronRight);

        public delegate void InitializeParametersSignature(Desertwind Solution, int CountStartingPoints);

        public delegate void IntializeSignature();

        public delegate void SequelSignature();

        public delegate bool SolveMatrixSignature(double[][] Answer, Desertwind Solution);

        public delegate void WriteSolutionSignature(string fname, int OutputOption, double[][] param, double[][] perr);

        #endregion

        // To do any processing after doing the main chore of MDS


        /* Input Parameters */
        public static int npar = 0; // Number of Parameters = Number_Vectors*VectorDimension
        public static int Number_VectorParameters = 0; // Number of Vectors
        public static int ParameterVectorDimension = 1; // Vector Dimension
        public static long ndata = 0; // Number of Data Points

        public static bool funcset = false; // If true the residuals are set in Desertwind
        public static bool fullmatrixset = false; // if true set full matrices in Desertwind

        public static int nbadgo = 6; // Stop if nbadgo tries produce less change than dellim
        public static int maxit = 20; // Maxiumum number of iterations

        public static double ChisqChangePerPoint = 0.001;
                             // Stop when predicted change in Chisq per point <= ChisqChangePerPoint

        public static double rho = 0.25; // Fletcher's parameter rho to define "good enough" changes in Chisq
        public static double sigma = 0.75; // Fletcher's parameter sigma to define "excellent" changes in Chisq
        public static double Omega = 0.25; // Parameter Omega to define "Excellent" changes in Chisq
        public static int OmegaOption = 1; // Control line search if Actual Change large

        public static double Qscale = 1.5; //   Force all  Q changes to be at least a factor(>1) of Sqrt(Qscale)
        public static double QHighInitialFactor = 0.01; // Initial call uses this fraction of Qhigh for Q (added 2009)
        public static double QgoodReductionFactor = 0.5; // Set to zero to make Q = 0. for good solutions

        public static int QLimitscalculationInterval = 6;
                          //  recalculate Qhigh Qlow Trace and Norm after this number of iterations

        public static int InitialSteepestDescents = 0;
                          // Start minimization with this number of Steepest Descent Iterations

        public static int DecisionMethod = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_1 = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_2 = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_3 = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_4 = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_5 = 0; // Specify why Manxcat made its choice
        public static int DecisionMethod_LineSearch = 0; // Specify why LineSearch made its choice

        public static bool UseDiagonalScalinginSolvers = true; // Surprisingly Original Manxcat did NOT Scale in Solvers
        //  The scaling is not applied in matrices but applied in calculation (this is a fixed input variable)
        // Original Manxcat did scale matrix to estimate Q and then added Q * diagonal element to matrix before solving equations

        public static bool AddMarquardtQExplicitly = false;
                           // If true add Marquardt explicitly in solvers (input variable)

        public static bool derivtest = false; // If true do a derivative test
        public static bool doingderivtest = false; // If true you are doing a derivative test
        public static bool FullSecondDerivative = false; // If true use full second derivative matrix

        public static double timmax = -1.0; // if positive, Stop when time used in milliseconds greater than timmax

        public static double CGResidualLimit = 0.00001; // Stop when Norm of Residual is less than this
        public static int PowerIterationLimit = 200; // Limit on Power Iterations
        public static double eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
        public static double eigenvectorchange = 0.001; // convergence test in power method for eigenvector
        public static double addonforQcomputation = 2.0; // Add this to make Power Method Safer
        public static double extraprecision = 0.05; //  Extra precision for Qlow

        public static int errcal = 1; // On input or output errcal = 0 means do not calculate errors
        // On input errcal = 1 means calculate errors
        // On output errcal = 1 means errors calculated successfully
        // On output errcal = -1 means error calculation tried failed

        public static string HotsunComment = ""; // Comment on Result from Hotsun

        //  Internal and Output Variables
        public static double zeromn = 0.0; // Minimum Chisq
        public static double zerocr = 0.0; // Current Chisq

        public static int idata = 0; // Count data points during update
        public static int numit = 0; // Count Iterations
        public static int bnderr = 0; // Number of Consequitive Boundary errors reported by Calcfg
        public static int bnderrLimit = 3; // Limit on Number of Consequitive Boundary errors reported by Calcfg

        public static int endres = 0; // Reason Fit finished
        // endres =1 Iteration Limit
        // endres =2 Converged Correctly
        // endres =3 Time Cut
        // endres =4 No Progress in nbadgo iterations
        // endres =5 Boundary Violation
        // endres =6 Matrix Singular even though Q added

        public static double tcalcfg = 0.0; // Total time in user routine in milliseconds
        public static double tsolve = 0.0; // Total time in fitting routine in milliseconds
        public static double teigen = 0.0; // Total time in Q range (eigen) routine in milliseconds
        public static double TotalTimeUsed = 0.0; // Total Time Used in milliseconds

        public static bool UseDiagonalScaling = false;
                           // if true, scale diagonals in this computation (this is dynamic variable)

        public static Boolean AddMarquardtQDynamically = false;
                              // If true  add in Marquardt Factor in current computation (this is dynamic variable)

        // In original Manxcat, the factor Q was added explicitly but given we don't multiply in scaling factors
        //  it is cleaner code not to add in Q but do both in Solvers as these are specialized anyway

        public static int NumberofCGIterations = -2; // Number of Iterations used in Conjugate Gradient Solver
        public static int CGIterationLimit = 200; // Number of Iterations used in Conjugate Gradient Solver

        public static int EigenvalueIndicator1 = 0; // Success Indicator in finding Qhigh
        public static int EigenvalueIndicator2 = 0; // Success Indicator in finding Qlow

        public static int TotalCGIterations = 0; // Accumulate Number of Conjugate Gradient Solver Iterations
        public static int TotalCGFailures = 0; // Accumulate Number of Conjugate Gradient Solver Failures
        public static int TotalPowerIterations = 0; // Accumulate Number of Eigenvalue Solver Iterations
        public static int TotalPowerFailures = 0; // Accumulate Number of Eigenvalue Solver Failures
        public static int TotalSearchIterations = 0; // Accumulate Number of Linear Search Solver Iterations

        public static bool DecomposeParameters = true;
        public static Desertwind CurrentSolution; // DECOMPOSED on rows
        public static Desertwind BestSolution; // DECOMPOSED on rows
        public static Desertwind PreviousSolution; // DECOMPOSED on rows
        public static Desertwind SearchSolution1; // DECOMPOSED on rows
        public static Desertwind SearchSolution2; // DECOMPOSED on rows
        public static Desertwind SearchSolution3; // DECOMPOSED on rows
        public static Desertwind SearchSolution4; // DECOMPOSED on rows

        public static Desertwind BeginningLinePositionSolution;
                                 // Set to starting point of search (holds param of starting point and value of first derivative there)

        public static Desertwind EndingLinePositionSolution;
                                 // Set to starting point of search (holds xshift from starting point and value and first derivative of end of line )


        //  The following arrays are decomposed on rows
        public static double[][] perr; // Errors in parameters
        // dgsave removed as Marquardt Q always added dynamically
        public static double[][] UtilityLocalVector1; // Utility Decomposed Vector
        public static double[][] UtilityLocalVector2; // Utility Decomposed Vector

        //  The following arrays are NOT decomposed on rows
        public static double[][] diag; // True diagonals
        //      public static double[][] sqdiag;  Square root of true diagonals removed in this version
        //      public static double[][] dginv;   // Inverse of Diagonals removed in this version
        public static double[][] sqdginv; // Square root of Inverse of Diagonals

        public static double[][] GlobalParameter; // Current Solution for all points
        public static double[][] UtilityGlobalVector1; // General non decomposed Vector used in utility fashion

        public static bool[][] FixedParameter;
                               // If true fix this parameter at initial value -- often zero -- Added 2009

        //  Local variables in original manxcat.f
        //  Means not needed in Chisq/derivative calculation or in initialization/final return values
        public static int materr = 0; // If nonzero specifies error code in matrix solver
        public static int isaved = 0; // Specify where best solution is
        // isaved = 0 best solution is in CurrentSolution
        // isaved = 1 best solution is in BestSolution
        // isaved = 2 best solution is in CurrentSolution except for residuals which are in BestSolution.func
        public static bool succ = false; // If true current fit reduced Chisq; if false it increased it
        public static int igood = 0; // Counts number of consequitive good iterations
        public static int CountToSD = 0; // Counts number of consequitive modest iterations
        public static int CountQSmall = 0; // Counts number of consequitive small Q's
        public static int CountQLarge = 0; // Counts number of consequitive large Q's
        public static int ResetQLimit = 10; // Reset Q if stays in same region this length of time
        public static int SDLimit = 10; // Choose Steepest Descent if CountToSD reaches this limit
        public static int RandomLimit = 15; // Choose Random move if irandomize reaches this limit
        public static int CountToRandom = 0; // Count towards explosive choice
        public static int CountExercusion = 0; // Count Random Exercusion
        public static int ExercusionLimit = 3; // Count Random Exercusion
        public static double QgoodFactor = 1.0; // Decrease Q by this factor when lowering below Qlow estimate
        public static int IterationforQlimit = 0; // Iteration at which Qlimits found

        public static double Q = 0.0; // Marquardt Lambda parameter

        public static double Qlow = 1.0;
                             //     Lower limit on Q (Estimate of minimum eigenvalue of scaled second derivative matrix)

        public static double Qhigh = 1.0;
                             //    Upper limit on Q (Estimate of maximum eigenvalue of scaled second derivative matrix)

        public static double Qgood = -1.0; //   Value of Q on last successful iteration
        public static double Qsd = 0.0; // Value of Q used in Steepest Descent Try
        public static double Qbest = 0.0; // Value of Q used in Best Solution
        public static double Qrange = 0.001; // Do not allow Qlow to be smaller than Qrange * Qhigh Added 2009
        public static double ChisqMatrixTrace = 0.0; // Trace of Chisq Matrix
        public static double ChisqMatrixNorm = 0.0; // Norm of Chisq Matrix

        public static double xnorm = 0.0; // Sum over points of Diag[point][vector[ * xshift[point][vector]**2
        public static int isdtry = 0; // Specify Steepest descent options
        // isdtry = 0 Steepest Descent not tried
        // isdtry = 1 Steepest Descent  to be tried next iteration
        // isdtry = 2 Steepest Descent just tried
        // isdtry = 3 Steepest Descent to be tried next iteration after Matrix Failure
        // isdtry = 4 Steepest Descent after Matrix Failure just tried

        public static double dellim = 0.0;
                             // Stop when predicted change in Chisq  <= dellim (calculated from ChisqChangePerPoint)

        public static double expchg = 0.0; //  Expected chisq change at next iteration
        public static double delchi = 0.0; //  Actual chisq change at last iteration
        public static double pred1 = 0.0; //  Expected chisq change at next iteration from linear term in expansion

        public static double pred2 = 0.0;
                             //  Expected chisq change at next iteration from quadratic term in expansion using Chisq Approximation

        public static double pred3 = 0.0;
                             //  Expected chisq change at next iteration from quadratic term in expansion using FULL formula

        public static double[] chsave = new double[nbadgo]; // Circular array to save previous nbadgo values of chisq
        public static int ichsav = 0; // Counts from 0 to nbadgo-1 in chsave positions

        public static int InitializationLoops = 1; // Number of Initial Conditions to be looped over
        public static int InitializationLoopCount = 0; // Count number of Initial Conditions to be looped over
        public static Desertwind BestLoopedSolution; // Solution gotten from previous loops (DECOMPOSED on rows)
        public static double[] InitLoopChisq; // Store Chisq on each Initialization Loop
        public static int BestChisqLoop = -1; // Loop on which best solution calculated

        public static int CoordinateWriteFrequency = 80;
                          // Frequency, in iterations, with which to write out MDS coordinates.

        public static void SetupManxcat()
        {
            //  Set Hotsun Input parameters
            maxit = ManxcatCentral.Configuration.Maxit;
            nbadgo = ManxcatCentral.Configuration.Nbadgo;
            ChisqChangePerPoint = ManxcatCentral.Configuration.ChisqChangePerPoint;
            rho = ManxcatCentral.Configuration.FletcherRho;
            Omega = ManxcatCentral.Configuration.Omega;
            OmegaOption = ManxcatCentral.Configuration.OmegaOption;
            sigma = ManxcatCentral.Configuration.FletcherSigma;
            QHighInitialFactor = ManxcatCentral.Configuration.QHighInitialFactor;
            QgoodReductionFactor = ManxcatCentral.Configuration.QgoodReductionFactor;
            QLimitscalculationInterval = ManxcatCentral.Configuration.QLimitscalculationInterval;
            InitialSteepestDescents = ManxcatCentral.Configuration.InitialSteepestDescents;
            timmax = ManxcatCentral.Configuration.TimeCutmillisec;
            CGResidualLimit = ManxcatCentral.Configuration.CGResidualLimit;
            PowerIterationLimit = ManxcatCentral.Configuration.PowerIterationLimit;
            eigenvaluechange = ManxcatCentral.Configuration.Eigenvaluechange;
            eigenvectorchange = ManxcatCentral.Configuration.Eigenvectorchange;
            extraprecision = ManxcatCentral.Configuration.Extraprecision;
            addonforQcomputation = ManxcatCentral.Configuration.AddonforQcomputation;
            derivtest = ManxcatCentral.Configuration.Derivtest;
            dellim = ChisqChangePerPoint; // Set dellim for current size problem (original Manxcat didn't do this)
            CoordinateWriteFrequency = ManxcatCentral.Configuration.CoordinateWriteFrequency;

            //  Parameters for Loop over Choices

            int NumberofPoints = Number_VectorParameters;
            InitializationLoops = ManxcatCentral.Configuration.InitializationLoops;
            InitLoopChisq = new double[InitializationLoops];

            diag = new double[NumberofPoints][];
            sqdginv = new double[NumberofPoints][];
            GlobalParameter = new double[NumberofPoints][];
            UtilityGlobalVector1 = new double[NumberofPoints][];
            FixedParameter = new bool[NumberofPoints][];

            for (int LongIndex = 0; LongIndex < NumberofPoints; LongIndex++)
            {
                diag[LongIndex] = new double[ParameterVectorDimension];
                sqdginv[LongIndex] = new double[ParameterVectorDimension];
                ;
                GlobalParameter[LongIndex] = new double[ParameterVectorDimension];
                UtilityGlobalVector1[LongIndex] = new double[ParameterVectorDimension];
                FixedParameter[LongIndex] = new bool[ParameterVectorDimension];

                for (int LocalVectorIndex = 0; LocalVectorIndex < ParameterVectorDimension; LocalVectorIndex++)
                {
                    FixedParameter[LongIndex][LocalVectorIndex] = false;
                }
            }

            int LocalNumberofPoints = Number_VectorParameters;

            if (DecomposeParameters)
                LocalNumberofPoints = SALSAUtility.PointCount_Process;

            perr = new double[LocalNumberofPoints][];
            UtilityLocalVector1 = new double[LocalNumberofPoints][];
            UtilityLocalVector2 = new double[LocalNumberofPoints][];

            for (int LongIndex = 0; LongIndex < LocalNumberofPoints; LongIndex++)
            {
                perr[LongIndex] = new double[ParameterVectorDimension];
                UtilityLocalVector1[LongIndex] = new double[ParameterVectorDimension];
                UtilityLocalVector2[LongIndex] = new double[ParameterVectorDimension];
            }

            CurrentSolution = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            BestSolution = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            PreviousSolution = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            BestLoopedSolution = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            SearchSolution1 = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            SearchSolution2 = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            SearchSolution3 = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
            SearchSolution4 = new Desertwind(LocalNumberofPoints, ParameterVectorDimension);
        }
    }

    // End Hotsun

    // Basic data for vector based parameters
    // func not stored as in general too big
    public class Desertwind
    {
        public double Chisquared; // Value of Chisq

        public double[][,] DiagonalofMatrix;
                           // Diagonal of Matrix as a matrix in VectorDimension for usual Chisq approximation; HALF Second Derivative

        public double[][,] ExactDiagonalofMatrix;
                           // Diagonal of Matrix as a matrix in VectorDimension without usual Chisq approximation; HALF Second Derivative

        public double[,][,] ExactFullMatrix; // Full matrix
        public double[,][,] FullMatrix; // Full matrix

        public int IterationCalculated; // Iteration solution calculated on
        public double[][] first; // HALF First Derivatives
        public double[][] func = null; // Chisq residuals

        public double[][] param;
                          // Parameter values (See xshift. These are incremented by MINUS xshift from previous solution)

        public double[][] xshift;
                          // Minus predicted change in parameter values -- This is Stored with Solution that was made by it

        public Desertwind(int NumberofPoints, int VectorDimension)
        {
            param = new double[NumberofPoints][];
            first = new double[NumberofPoints][];
            xshift = new double[NumberofPoints][];

            if (Hotsun.fullmatrixset)
            {
                FullMatrix = new double[NumberofPoints,NumberofPoints][,];
                ExactFullMatrix = new double[NumberofPoints,NumberofPoints][,];
            }
            else
            {
                DiagonalofMatrix = new double[NumberofPoints][,];
                ExactDiagonalofMatrix = new double[NumberofPoints][,];
            }

            for (int LongIndex = 0; LongIndex < NumberofPoints; LongIndex++)
            {
                if (Hotsun.fullmatrixset)
                {
                    for (int OtherIndex = 0; OtherIndex < NumberofPoints; OtherIndex++)
                    {
                        FullMatrix[LongIndex, OtherIndex] = new double[VectorDimension,VectorDimension];
                        ExactFullMatrix[LongIndex, OtherIndex] = new double[VectorDimension,VectorDimension];
                    }
                }
                else
                {
                    DiagonalofMatrix[LongIndex] = new double[VectorDimension,VectorDimension];
                    ExactDiagonalofMatrix[LongIndex] = new double[VectorDimension,VectorDimension];
                }
                param[LongIndex] = new double[VectorDimension];
                first[LongIndex] = new double[VectorDimension];
                xshift[LongIndex] = new double[VectorDimension];
            }

            if (Hotsun.funcset == false)
                func = null;
            Chisquared = 0.0;
        }
    }

    // end Desertwind
}

// end Namespace Manxcat