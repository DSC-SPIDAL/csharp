using System;
using System.Threading;
using System.Threading.Tasks;
using MPI;

using SALSALibrary;

namespace Salsa.PairwiseClusteringTPL
{
    //	Class Dist **********************************************************
    // ------------------------------------------------------------------------------
    //	getDist controls complete pairwise computation
    //	It sets up problem and calls parallel threads PairwiseThread. 

    //	PairwiseThread does major pairwise computation looping over all EM iterations

    //	In first iteration it is assumed that Malpha(k) is set  and epsilonalpha(k) is calculated
    //	Other iterations use value of epsilonalpha(k) to first calculate a new value of Malpha(k)
    //	One calls calculatEpsi to find A(k) Balpha(k) from C(k) Malpha(k) and distance matrix
    //	Then Program.convergenceTest in Program decides if current EM iteration converged
    //	Program.ShouldWeSplit is called to see if clusters should be split
    //	Program.DoTheSplit is called to split chosen cluster
    //	Dist.Temperature is reduced and EM Iteration continues
    //	
    //	calculateEpsi finds A(k) Balpha(k) from C(k) Malpha(k) and distance matrix
    //	These give the new value of epsilonalpha(k)

    // ------------------------------------------------------------------------------
    // public void getDist()
    // Implements logic of EM and Temperature loop and decision on splitting
    // 
    // ------------------------------------------------------------------------------
    // public static void PrintIteration(string looptype)
    // Print summary of EM iteration with a given temperature value
    // 
    // ------------------------------------------------------------------------------
    // public void PairwiseThread()
    // Control calculation of key functions starting by calculating M and then calling
    // calculateEpsi(Malpha_k_, A_k_, Balpha_k_, C_k_, Epsilonalpha_k_, Dist.Ncent)
    // 
    // ------------------------------------------------------------------------------
    // public void calculateEpsi(double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, double[][] localepsi, int localNcent)
    //	Perform multiple parallel steps calculating A_k_, Balpha_k_, epsi and differences
    //	Value of M is assumed
    // Complicated MPI parallelism broadcasting M in parts with a ring
    // 
    // ------------------------------------------------------------------------------
    // public static int convergenceTest(double ChangeLimit)
    // The change of sum(points) Delta(epsi)/# Points should be less than argument ChangeLimit for each cluster
    // epsiDiff calculated in CalculateEpsi while Epsi are updated
    //  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
    // 
    // Test if iteration has converged for given temperature
    // 
    // ------------------------------------------------------------------------------
    // public void SaveCurrentTask()
    // {   // Save current task that should be split but we need to converge first
    // public void RestorePreviousTask()
    // {   // Restore previous task that should be split but we needed to converge current cluster configuration first
    // 
    // Currently before splitting we converge current cluster configuration to zero Temperature before saving
    // These routines allow one to save configuration before split and return to it
    // 
    // ------------------------------------------------------------------------------
    // JiggleClusters() // Jiggle Values of Epsilon -- leaving other Cluster Parameters
    // 
    // Either do totally randomly or along most negative eigenvector Based on Program.JigglePerturbation
    // Eigenvalue analysis called from JiggleClusters()
    // 
    // Set new values ofd epsilon(alpha, k).
    // Calculate change in C and M for output only
    // Set flag that M has to be calculated
    // 
    // ------------------------------------------------------------------------------
    // public bool shouldweSplit()
    // Decide if to split based on negative eigenvalue of second derivative matrix
    // MinimumEigenvalue is MINIMUM eigenvalue
    // ClustertoSplit is cluster number of this
    // D is distance matrix
    //  Skip first time this happens as initialize to two clusters
    // 
    // Control by Program.Eigenvalue_Methodology
    // =0 Do not look at second derivative for splitting
    // =1 original simple method using unaltered second derivative matrix (use for Continuous Clustering)
    // =2 quick and dirty estimate of Duplicated second derivative matrix
    // =3 EM estimate of Duplicated second derivative matrix (DEFAULT). Each cluster has weight increased to 2 and system iterated to convergence
    // =4 Sophisticated estimate of Duplicated second derivative matrix (used in test only as 3 faster and same)
    // 
    // Program.PerformEigenTest is defaulted to False.
    // If true perform a test of approach
    // 
    // shouldweSplit() calls eigenvalue routine
    // 
    // ------------------------------------------------------------------------------
    // public void dothesplit()
    //  Do a split on the identified cluster found in shouldweSplit()
    //  Malpha_k_ is value of M indexed by (point,cluster)
    //  New clusters stored in last cluster and old position
    //  Cluster count is incremented by one
    //  ClustertoSplit is cluster number to split
    //	Embarassingly Parallel over Processes -- no MPI calls needed
    // 
    // public static int PerturbationVehicle = 0; // = 0 Perturb Epsilon ; if 1 Perturb Malpha_k_
    // public static double SplitPerturbationFactor = 1.0; // Normalization of Perturbation for Split
    // 
    // M is halved and either M or epsilon perturbed
    // epsilon uses eigenvector for perturbation direction (assumes currently that calculated)
    // 
    // -------------------------------------------------------------------------------
    // public static void initializeTemperature()
    //	Find initial Temperature that can be used to start with one cluster
    // Just need an upper bound as whips through this part ofd program
    // 
    // ------------------------------------------------------------------------------
    // Calculate Distance Correlation Function
    // public void CorrelationCalculation()
    // Use at end of run for final clusters

    // Calculate Pairwise Algorithm

    public class ClusteringSolution
    {
        public double[][] Old_Epsilonalpha_k_;      // Previous value of Epsilon
        public double[][] Epsilonalpha_k_;          // Epsilon
        public double[][] Best_Epsilonalpha_k_;     // best Epsilon in loop over Clusters Duplicated
        public double[][] Master_Epsilonalpha_k_;   // Master Epsilon in loop over Clusters Duplicated
        public double[][] Malpha_k_;                // Probability that point in Cluster k 
        public double[][] Previous_Malpha_k_;       // Previous value of Malpha_k_ used to test converge of p(k) Malpha(k) loop
        public double[][] Balpha_k_;                // B(alpha,k) = sum_N d(ClusterCenter, a) M_i(k) / C(k)

        public double[] A_k_;                       // A(k) = -0.5 sum_N Balpha(k) M_a(k) / C(k)
        public double[] C_k_;                       // Summation of C(k) over parallel threads
        public double[] Weight_k_;                  // weight -- typically 1 -- of each cluster
        public double[] FreezingMeasure_k_;         // Freezing Measure for Cluster
        public double[] P_k_;                       // P(k) used in Continuous Clustering

        public double[][] Saved_Ax;                 // New power vector
        public double[][] Saved_oldAx;              // Old power vector

        public double PairwiseHammy = 0.0;                // Value of Hamiltonian
        public double OldHammy = 0.0;                     //  Previous value of Hamiltonian
        public int Ncent = 1;                           //  the current number of clusters
        public int ClustertoSplit = -1;                  //  the cluster label that will be split
        public double MinimumEigenvalue = 0.0;            // Minimum full eigenvalue associated with ClustertoSplit
        public double ActualCoolingFactor = Program.InitialCoolingFactor;          // Actual Cooling Factor
        public double Temperature = Dist.Tinitial;                  // The current temperature
        public int IterationSetAt = 0;                      // Iteration set at
        public bool Axset = false;                          // If True Saved_Ax and Saved_oldAx set
        public bool SolutionSet = false;                    // If True Solution Set

        public static int NumberofPointsinProcess = -1;  // Number of Points in Process
        public static int MaximumNumberClusters = 0;     // Maximum Number of Centers
        public static int cachelinesize = 0;            // Cacheline increment
       
        public static void SetParameters(int NumberofPointsinProcessINPUT, int MaximumNumberClustersINPUT, int cachelinesizeINPUT)
        {
            if (NumberofPointsinProcess > 0)
                return;
            NumberofPointsinProcess = NumberofPointsinProcessINPUT;
            MaximumNumberClusters = MaximumNumberClustersINPUT;
            cachelinesize = cachelinesizeINPUT;
        }

        public ClusteringSolution()
        {
            if (NumberofPointsinProcess < 0)
            {
                Exception e = PWCUtility.SALSAError("NumberofPointsinProcess Unset");
                throw (e);
            }
            Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
            Old_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
            Best_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
            Master_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
            Malpha_k_ = new double[NumberofPointsinProcess][];
            Previous_Malpha_k_ = new double[NumberofPointsinProcess][];
            Balpha_k_ = new double[NumberofPointsinProcess][];

            Saved_Ax = new double[NumberofPointsinProcess][];
            Saved_oldAx = new double[NumberofPointsinProcess][];

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Old_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Best_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Master_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Malpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Previous_Malpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Balpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];  //Balpha[k] = Sum(ClusterCenter) D[ClusterCenter,alpha]Malpha[k]/C[k]

                    Saved_Ax[ProcessPointIndex] = new double[MaximumNumberClusters + cachelinesize];
                    Saved_oldAx[ProcessPointIndex] = new double[MaximumNumberClusters + cachelinesize];

                    for (int ClusterIndex = 0; ClusterIndex < Program.maxNcent; ClusterIndex++)
                    {
                        Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = 0.0;
                        Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = 0.0;
                    }
                }
            }); // End loop initialing Point dependent quantities

            C_k_ = new double[MaximumNumberClusters + cachelinesize]; // Final value of C(k)
            A_k_ = new double[MaximumNumberClusters + cachelinesize];
            P_k_ = new double[MaximumNumberClusters + cachelinesize];

            FreezingMeasure_k_ = new double[MaximumNumberClusters + cachelinesize];
            Weight_k_ = new double[MaximumNumberClusters + cachelinesize];

            Axset = false;
        }

        public static void CopySolution(ClusteringSolution From, ClusteringSolution To)
        {
            To.PairwiseHammy = From.PairwiseHammy; // Value of Hamiltonian
            To.OldHammy = From.OldHammy;  //  Previous value of Hamiltonian
            To.Ncent = From.Ncent;    //the current number of clusters
            To.ClustertoSplit = From.ClustertoSplit;   //the cluster label that will be split
            To.MinimumEigenvalue = From.MinimumEigenvalue; // Minimum full eigenvalue associated with ClustertoSplit
            To.ActualCoolingFactor = From.ActualCoolingFactor;  // Actual Cooling Factor
            To.Temperature = From.Temperature; //the current temperature
            To.Axset = From.Axset;
            To.SolutionSet = From.SolutionSet;
            To.IterationSetAt = From.IterationSetAt;

            int NumberClusters = From.Ncent;
            for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
            {
                To.C_k_[ClusterIndex] = From.C_k_[ClusterIndex];
                To.A_k_[ClusterIndex] = From.A_k_[ClusterIndex];
                To.P_k_[ClusterIndex] = From.P_k_[ClusterIndex];

                To.FreezingMeasure_k_[ClusterIndex] = From.FreezingMeasure_k_[ClusterIndex];
                To.Weight_k_[ClusterIndex] = From.Weight_k_[ClusterIndex]; 
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {

                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                    {
                        To.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = From.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = From.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = From.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = From.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Malpha_k_[ProcessPointIndex][ClusterIndex] = From.Malpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] = From.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Balpha_k_[ProcessPointIndex][ClusterIndex] = From.Balpha_k_[ProcessPointIndex][ClusterIndex];


                        if (From.Axset)
                        {
                            To.Saved_Ax[ProcessPointIndex][ClusterIndex] = From.Saved_Ax[ProcessPointIndex][ClusterIndex];
                            To.Saved_oldAx[ProcessPointIndex][ClusterIndex] = From.Saved_oldAx[ProcessPointIndex][ClusterIndex];
                        }
                    }
                }
            }); // End loop copying Point dependent quantities
        }   // End CopySolution

        public static void RemoveCluster(ClusteringSolution Changing, int RemovedIndex)
        {
            --Changing.Ncent;
            if (RemovedIndex == Changing.Ncent)
                return;

            for (int ClusterIndex = RemovedIndex; ClusterIndex < Changing.Ncent; ClusterIndex++)
            {
                Changing.C_k_[ClusterIndex] = Changing.C_k_[ClusterIndex+1];
                Changing.A_k_[ClusterIndex] = Changing.A_k_[ClusterIndex+1];
                Changing.P_k_[ClusterIndex] = Changing.P_k_[ClusterIndex+1];

                Changing.FreezingMeasure_k_[ClusterIndex] = Changing.FreezingMeasure_k_[ClusterIndex+1];
                Changing.Weight_k_[ClusterIndex] = Changing.Weight_k_[ClusterIndex+1];
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {

                    for (int ClusterIndex = 0; ClusterIndex < Changing.Ncent; ClusterIndex++)
                    {
                        Changing.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex]+1;
                        Changing.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex+1];
                        Changing.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex+1];
                        Changing.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex+1];
                        Changing.Malpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Malpha_k_[ProcessPointIndex][ClusterIndex+1];
                        Changing.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex+1];
                        Changing.Balpha_k_[ProcessPointIndex][ClusterIndex] = Changing.Balpha_k_[ProcessPointIndex][ClusterIndex+1];


                        if (Changing.Axset)
                        {
                            Changing.Saved_Ax[ProcessPointIndex][ClusterIndex] = Changing.Saved_Ax[ProcessPointIndex][ClusterIndex+1];
                            Changing.Saved_oldAx[ProcessPointIndex][ClusterIndex] = Changing.Saved_oldAx[ProcessPointIndex][ClusterIndex+1];
                        }
                    }
                }
            }); // End loop copying Point dependent quantities
        }   // End RemoveCluster

        public static void SetAxinSolution(ClusteringSolution Solution)
        {
            Solution.Axset = true;
            int NumberClusters = Solution.Ncent;
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                    {
                        Solution.Saved_Ax[ProcessPointIndex][ClusterIndex] = vectorclass.Ax[ProcessPointIndex][ClusterIndex];
                        Solution.Saved_oldAx[ProcessPointIndex][ClusterIndex] = vectorclass.oldAx[ProcessPointIndex][ClusterIndex];
                    }
                }
            }); // End loop copying Point dependent quantities
        }

        public static void RestoreAxfromSolution(ClusteringSolution Solution)
        {
            if (Solution.Axset == false)
            {
                Exception e = PWCUtility.SALSAError("Axset false but restore requested");
                throw (e);
            }

            int NumberClusters = Solution.Ncent;
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++)
                    {
                        vectorclass.Ax[ProcessPointIndex][ClusterIndex] = Solution.Saved_Ax[ProcessPointIndex][ClusterIndex];
                        vectorclass.oldAx[ProcessPointIndex][ClusterIndex] = Solution.Saved_oldAx[ProcessPointIndex][ClusterIndex];
                    }
                }
            }); // End loop copying Point dependent quantities
        }

    }   // End ClusteringSolution

    class Dist
    {
        public static double[] Diff_Epsilon_k_ = null;

        public static double[,] DistanceCorrelation;    // Distance Correlation
        public static int[] ClusterSelected; // NONZERO values are used to select Eigenvalues to Process

        public static ClusteringSolution RunningPWC;    // The main Solution
        public static ClusteringSolution SavedPWC;      // Saved Solution
        public static ClusteringSolution BestPWC;       // Current Best Solution

        public static int PlaceforEigsforMultipleCluster = -1;  // Where eigenvalues for complex situations stored
        public static int[] ListofClusterstoSplit;  // List of Clusters to Split
        public static double[] EigsofClusterstoSplit;  // List of Eigenvalues of Clusters to Split
        public static int Numberthatcanbesplit = 0; // Length of entries in ListofClusterstoSplit which is at most Program.MaxNumberSplitClusters

        public static double Tmin;  //minimal temperature 
        public static double Tinitial;  //Initial temperature
        public static int ActualMaxNcent;   // Maximum reduced if small clusters
        public static int CountValidityFailures = 0;    // Count failures in validity checks of final solution
        public static bool OnLastleg = false;   // If True we are on final convergence step with full number of clusters

        public static int EMIterationCount = 0; // iteration number in EM method
        public static int EMIterationStepCount = 0;
        public static int Extra_EMIterationCount = 0; // Extra iterations in EM method
        public static int ActualWaititerations = 0;
        public static int countAfterFixingClusterCount = 0;   // Count loops after end on number of centers
        public static int PMLoopUsed = 0;   // Number of PM Loops
        public static int IterationNumberPrinted = -1;             // Iteration Number Printed

        public static int oldepsiset = 0;	// oldepsiset =1 if oldepsi set
        public static int needtocalculateMalpha_k_ = -1;	// =0 Already Set; = 1 Calculate from EM; = -1 Initialize

        public static bool HitConvergenceLoopLimit = false;
        public static bool justconverging = false;  // If true just converging -- No splits
        public static bool taskswaiting = false;  // If true just converging -- No splits
        public static int SplitFailures = 0;    // Counts splitting failures

        public static bool initialized = false;
        public static bool HammyNotSet = true;  // If True Hamiltonian Set


        /*
         * The basic function in the program. It implements all steps of the pairwiseDA algorithm.
         *  
         */
        public void getDist()
        {

            //allocate memory on first and indeed only call
            if (!initialized)
            {
                initialized = true;

                ClusteringSolution.SetParameters(PWCUtility.PointCount_Process, Program.maxNcent, Program.cachelinesize);
                RunningPWC = new ClusteringSolution();
                SavedPWC = new ClusteringSolution();
                BestPWC = new ClusteringSolution();

                Program.MaxNumberSplitClusters = Math.Max(Program.MaxNumberSplitClusters, 1);
                Dist.ListofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
                Dist.EigsofClusterstoSplit = new double[Program.MaxNumberSplitClusters];
                
                int cachelinesize = Program.cachelinesize;

                ClusterSelected = new int[Program.maxNcent + cachelinesize];


            }	//end Initialization

            //  Do EM calculation
            Dist.ActualMaxNcent = Program.maxNcent;
            Dist.OnLastleg = false;

            // Set Initial Dist.Temperature
            initializeTemperature();
            Dist.RunningPWC.Temperature = Dist.Tinitial;
            Dist.Tmin = Dist.RunningPWC.Temperature / 100000.0;
            Dist.RunningPWC.ActualCoolingFactor = Program.InitialCoolingFactor;

            Dist.RunningPWC.Ncent = Program.InitialNcent;

            Dist.ActualWaititerations = Program.Waititerations; // Wait this number of Temperature iterations before splitting
            // This is changed upto 8 times input value if no cluster to split

            EMIterationCount = 0; // Initialize iteration count -- there is no limit
            EMIterationStepCount = 0;   // This is number of counts in current step converging epsilon at given temperature
            Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached (this decreases temperature)
            int HammyViolations = 0;    // Counts number of (illegal) consequitive INCREASES in Hamiltonian
            int CountBetweenJiggles = 0;    // Count iterations in Jiggle with maximum JiggleOption (Only used if JiggleOption > 0)
            int CountBetweenSplits = 0; // Count iterations between splits. Limit i ActualWaitIterations

            //  Variables to control pending task spun off when converging current case
            taskswaiting = false;   // If true we saved a Solution in SavedPWC and evolved that solution to define output at a particular cluster count
            bool currenttaskfinished = false;
            bool TakePreviousSplitDecision = false; // If true, we decided on a split but spun off a convergence task

            Dist.needtocalculateMalpha_k_ = -1; // =0 Malpha_k_ Already Set; = 1 Calculate Malpha_k_ from EM (usual); = -1 Initialize Malpha_k_

            //  Variables tracking convergence
            HitConvergenceLoopLimit = false;// If true, EMIterationStepCount has hit limit
            SplitFailures = 0;
            justconverging = false; // If true we are in a special loop freezing current case
            Dist.RunningPWC.OldHammy = Tinitial * PWCUtility.PointCount_Global * PWCUtility.PointCount_Global;  // Set Old value of Hamiltonian

            //	Loop over EM calculations
            while (true)
            {

                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                {
                    Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
                }

                Dist.EMIterationCount++; // Increment EM Loop Count
                Dist.EMIterationStepCount++;    // Iterate within a fixed temperature converging epilson
                if (Program.EMIterationStepCountMax >= 0)
                    Program.EMIterationStepCountMax = Math.Max(Program.EMIterationStepCountMax, Dist.EMIterationStepCount);
                else
                    Program.EMIterationStepCountMax = Dist.EMIterationStepCount;

                PairwiseThread();
                Dist.RunningPWC.SolutionSet = true;
                Dist.RunningPWC.IterationSetAt = Dist.EMIterationCount;

                //	Now see if we are done -- take results from Rank 0
                int convergence = 0;
                convergence = Dist.convergenceTest(Program.Epsi_max_change);
                PWCUtility.SynchronizeMPIvariable(ref convergence);

                //  If not converged at this temperature and cluster number just proceed with while(true) iteration
                if (convergence == 0)
                    continue;

// Case we are converged = 1 in middle or = 2 at end of refinement
                Dist.EMIterationStepCount = 0;

// Calculate measures of Clustering and save features for later printing
                double[] Save_ToPrint_C_k_ = new double[Dist.RunningPWC.Ncent];
                double[] Save_ToPrint_FreezingMeasure_k_ = new double[Dist.RunningPWC.Ncent];
                double[] Save_ToPrint_A_k_ = new double[Dist.RunningPWC.Ncent];
                double[] Save_ToPrint_P_k_ = new double[Dist.RunningPWC.Ncent];

                Dist.RunningPWC.PairwiseHammy = 0.0;
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                {
                    Dist.RunningPWC.PairwiseHammy += -2.0 * Dist.RunningPWC.A_k_[ClusterIndex] * Dist.RunningPWC.C_k_[ClusterIndex];
                    Save_ToPrint_C_k_[ClusterIndex] = Dist.RunningPWC.C_k_[ClusterIndex];
                    Save_ToPrint_FreezingMeasure_k_[ClusterIndex] = Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex];
                    Save_ToPrint_A_k_[ClusterIndex] = Dist.RunningPWC.A_k_[ClusterIndex];
                    Save_ToPrint_P_k_[ClusterIndex] = Dist.RunningPWC.P_k_[ClusterIndex];
                }

                //  Set Best Solution and progress flag decreasing
                bool decreasing;
                if (HammyNotSet)
                {
                    HammyNotSet = false;
                    decreasing = true;
                }
                else
                    decreasing = Dist.RunningPWC.PairwiseHammy <= Dist.RunningPWC.OldHammy;
                PWCUtility.SynchronizeMPIvariable(ref decreasing);
                if (decreasing)
                    ClusteringSolution.CopySolution(Dist.RunningPWC, Dist.BestPWC);

                double ChangeinHammy = Dist.RunningPWC.PairwiseHammy - Dist.RunningPWC.OldHammy;
                Dist.RunningPWC.OldHammy = Dist.RunningPWC.PairwiseHammy;

                //  Case when at end of stage (e.g. given number of clusters or of whole job). This still needs to be iterated to low Temperatures
                if (convergence == 2 || justconverging)
                {
                    if (Dist.countAfterFixingClusterCount < Program.Iterationatend)         //    do Iterationatend iterations after reaching the maximum cluster#
                    {
                        if (Dist.countAfterFixingClusterCount == 0)
                        {   // First step of "justconverging stage"
                            HammyViolations = 0;
                            Dist.RunningPWC.ActualCoolingFactor = Program.ConvergingCoolingFactor;
                        } 
                        int toobigfreezing = 0;
                        for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                        {
                            if (Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex] > Program.FreezingLimit)
                                ++toobigfreezing;
                        }
                        PWCUtility.SynchronizeMPIvariable(ref toobigfreezing);
                        Dist.countAfterFixingClusterCount++;

                        string looptype = "Final Clean Up";
                        if (taskswaiting)
                            looptype = "Intermediate Cluster Convergence Steps";
                        int printtype = -1;
                        if ((Dist.countAfterFixingClusterCount == (Program.Iterationatend - 1))
                            || (toobigfreezing == 0) || (!decreasing))
                            printtype = Dist.RunningPWC.IterationSetAt;
                        PrintIteration(looptype, printtype);

                        if (decreasing)
                        {
                            HammyViolations = 0;
                        }
                        else
                        {
                            ++HammyViolations;
                            if (HammyViolations < 5)
                                decreasing = true;
                        }
                        PWCUtility.SynchronizeMPIvariable(ref decreasing);
                        if (PWCUtility.MPI_Rank == 0)
                            currenttaskfinished = !((toobigfreezing > 0) && decreasing);
                        PWCUtility.SynchronizeMPIvariable(ref currenttaskfinished);

                        // Iteration at end still going Reduce Temperature and proceed
                        if (!currenttaskfinished)
                        {   
                            Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
                            ++CountBetweenJiggles;
                            if (CountBetweenJiggles > Program.JiggleOption)
                                CountBetweenJiggles = 0;
                            if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
                            {   // Jiggle Parameters
                                JiggleClusters();
                            }
                            continue;   // Continue while over EM Loop
                        }

                        // We have finished this task
                        if (!decreasing)
                            PWCUtility.SALSAPrint(1, "Stop as Hamiltonian Increasing with Change " + ChangeinHammy.ToString("E4"));
                        else
                            PWCUtility.SALSAPrint(1, "Stop as Freezing Measures smaller than " + Program.FreezingLimit.ToString());

                    }   // End Convergence=2  or justconverging doing Final Program.Iterationatend iterations counted by countAfterFixingClusterCount

                    else
                    {   // Case when Program.Iterationatend iterations counted by countAfterFixingClusterCount exceeded
                        PWCUtility.SALSAPrint(1, " Stop as Final Iteration Count larger than " + Program.Iterationatend.ToString());
                        currenttaskfinished = true;
                    }   // End Convergence=2 or justconverging case where extra iteration count completed
                    

                    // This completes processing for convergence=2 and justconverging=true case
                    if (currenttaskfinished)
                    {
                        if (Dist.RunningPWC.Temperature < Dist.Tmin)
                            PWCUtility.SALSAPrint(1, "Tmin Reached " + Dist.Tmin.ToString("E4"));
                        if (Dist.HitConvergenceLoopLimit == true)
                            PWCUtility.SALSAPrint(1, "EM Convergence Loop Reached for fixed Temperature -- Terminate Run");
                        if (Dist.SplitFailures > 1)
                            PWCUtility.SALSAPrint(1, "Consequitive Split Failures");
                        if (!taskswaiting)
                        {
                            if (!Dist.CheckValidSolution())
                            {   // Restart with small clusters removed
                                OnLastleg = true;
                                Dist.needtocalculateMalpha_k_ = 1;
                                Dist.countAfterFixingClusterCount = 0;
                                CountBetweenJiggles = 0;
                                CountBetweenSplits = 0;
                                HammyNotSet = true;
                                justconverging = true;
                                taskswaiting = false;
                                currenttaskfinished = false;
                                TakePreviousSplitDecision = false;
                                continue;
                            }
                            break;
                        }

                        //  Resume Previous Task
                        PWCUtility.SALSAPrint(1, "\nResuming Previous Task -- Implement Cluster Splitting");

                        //  Output Converged Results
                        if (Program.ClusterCountOutput > 0)
                        {
                            int[] counts = new int[Dist.RunningPWC.Ncent];
                            Program.OutputClusterLabels(counts);
                            PWCUtility.SALSAPrint(1, "Clusters Output " + Dist.RunningPWC.Ncent.ToString());
                        }

                        //  Restore Previous Task
                        RestorePreviousTask();
                        Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached
                        Dist.needtocalculateMalpha_k_ = 1; // =0 Already Set; = 1 Calculate from EM (usual); = -1 Initialize
                        CountBetweenJiggles = 0;
                        CountBetweenSplits = 0;
                        justconverging = false;
                        taskswaiting = false;
                        Dist.HitConvergenceLoopLimit = false;
                        currenttaskfinished = false;
                        TakePreviousSplitDecision = true;
                        convergence = 1;
                    }   // end case when task finished

                }   // End justconverging = true or convergence =2 case

                //  This section results either in EM Loop Continue (for ongoing refinement for a nonsplittable case) or EM Loop Break (to end) OR
                //  Restart previous task which has convergence=1 by definitio

                //  Convergence = 1 Case
                //	Converged so test for split -  -- take results from Rank 0 but rest do arithmetic
                //  For restarted task that decision has already been made
                bool ResultofSplittingTest = false;
                if (!TakePreviousSplitDecision)
                {   // Need to decide if split to occur
                    ++CountBetweenSplits;
                    if (CountBetweenSplits > Dist.ActualWaititerations)
                        CountBetweenSplits = 0;
                    if ((CountBetweenSplits > 0) && (Dist.ActualWaititerations > 0) && (Dist.RunningPWC.Ncent > 1))
                    {
                        PrintIteration("Ongoing Annealing", -1);
                        Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
                        ++CountBetweenJiggles;
                        if (CountBetweenJiggles > Program.JiggleOption)
                            CountBetweenJiggles = 0;
                        if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
                        {   // Jiggle Parameters
                            JiggleClusters();
                        }
                        continue;   // EM Loop for compulsory iterations between splits
                    }

                    //  Decide if to split
                    PWCUtility.StartSubTimer(0);
                    ResultofSplittingTest = shouldweSplit();
                    PWCUtility.SynchronizeMPIvariable(ref ResultofSplittingTest);
                    PWCUtility.SynchronizeMPIvariable(ref Dist.RunningPWC.ClustertoSplit);
                    PWCUtility.StopSubTimer(0);
                    PWCUtility.InterimTiming();

                    //  Diagnostic Output for splitting
                    if( ResultofSplittingTest || (Dist.EMIterationCount % Program.PrintInterval == 0) ) // *********************
                        diagnosticsplitprint(ResultofSplittingTest, Save_ToPrint_A_k_, Save_ToPrint_C_k_, Save_ToPrint_FreezingMeasure_k_, Save_ToPrint_P_k_);

                }   // End !TakePreviousSplitDecision -- case when splitting needs to be determined -- rather than taken from old task

                else
                {   // Case where splitting determined earlier
                    ResultofSplittingTest = true;
                }

                //	If split indicated perform this                          
                if (ResultofSplittingTest == true)
                {
                    if (Dist.RunningPWC.Ncent > 1)
                    {
                        if (!TakePreviousSplitDecision)
                        {
                            if (Program.ConvergeIntermediateClusters)
                            {   // Start iteration to converge current cluster output
                                SaveCurrentTask();
                                justconverging = true;
                                taskswaiting = true;
                                TakePreviousSplitDecision = true;
                                Dist.ActualWaititerations = Program.Waititerations;
                                Dist.RunningPWC.ActualCoolingFactor = Program.ConvergingCoolingFactor;
                                Dist.countAfterFixingClusterCount = 0;
                                CountBetweenSplits = 0;
                                CountBetweenJiggles = 0;
                                continue;
                            }
                            if ( (Program.ClusterCountOutput > 0) && Program.ConvergeIntermediateClusters)
                            {
                                int[] counts = new int[Dist.RunningPWC.Ncent];
                                Program.OutputClusterLabels(counts);
                            }
                        }
                    }

                    //  Need to perform already determined split
                    TakePreviousSplitDecision = false;
                    if (Dist.Numberthatcanbesplit > 1)
                    {
                        string nextline = "";
                        for (int splitlistloop = 0; splitlistloop < Dist.Numberthatcanbesplit; splitlistloop++) 
                        {
                            Dist.RunningPWC.ClustertoSplit = Dist.ListofClusterstoSplit[splitlistloop];
                            nextline += Dist.RunningPWC.ClustertoSplit.ToString() + " " + Dist.EigsofClusterstoSplit[splitlistloop].ToString("E4") + " * ";
                            dothesplit();
                        }
                        PWCUtility.SALSAPrint(1, " Multiple Clusters Split Center Eigenvalue " + nextline);
                        Dist.Numberthatcanbesplit = 0;
                    }
                    else
                        dothesplit();
                    Dist.SplitFailures = 0;
                    Dist.ActualWaititerations = Program.Waititerations;
                }   // End ResultofSplittingTest == true (either changed number of clusters or spun off a convergence task)

                //  Final portion of loop changing Temperature if needed
                if (ResultofSplittingTest == false)
                {   //	Reduce T and continue iterating if converged and no split
                    Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
                    ++CountBetweenJiggles;
                    if (CountBetweenJiggles > Program.JiggleOption)
                        CountBetweenJiggles = 0;
                    if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
                    {   // Jiggle Parameters
                        JiggleClusters();
                    }
                }
                else
                {   // Don't reduce T as clusters split
                    Dist.RunningPWC.ActualCoolingFactor = Program.FineCoolingFactor;
                    CountBetweenSplits = 0;
                    CountBetweenJiggles = 0;
                }

            }   // End while EM Loop

            //  Check if solution best
            if (!Dist.BestPWC.SolutionSet)
                return;
            if (Dist.BestPWC.Ncent != Dist.RunningPWC.Ncent)
            {
                PWCUtility.SALSAPrint(1, " Best Solution Not Used as Number of Centers " + Dist.BestPWC.Ncent + " Different from Running Solution with " + Dist.RunningPWC.Ncent);
                return;
            }
            bool changesolution = Dist.BestPWC.PairwiseHammy < Dist.RunningPWC.PairwiseHammy;
            PWCUtility.SynchronizeMPIvariable(ref changesolution);
            if (changesolution)
            {
                PWCUtility.SALSAPrint(1, " Solution at Iteration " + Dist.BestPWC.IterationSetAt.ToString() + " Chisq " + Dist.BestPWC.PairwiseHammy.ToString("E4")
                    + " Taken rather than Iteration " + Dist.RunningPWC.IterationSetAt.ToString() + " Chisq " + Dist.RunningPWC.PairwiseHammy.ToString("E4"));
                ClusteringSolution.CopySolution(Dist.BestPWC, Dist.RunningPWC);
            }
            Dist.PrintIteration(" Final Solution ", Dist.RunningPWC.IterationSetAt);

        }	// End of getDist

        public static void PrintIteration(string looptype, int linecheck)
        {
            if (linecheck < 0)
            {
                if (Dist.EMIterationCount % Program.PrintInterval != 0)
                    return;
            }
            else if (linecheck == Dist.IterationNumberPrinted)
                return;
            Dist.IterationNumberPrinted = Dist.RunningPWC.IterationSetAt;

            PWCUtility.InterimTiming();
            PWCUtility.SALSAPrint(1, "B) " + Dist.RunningPWC.Ncent.ToString() + " Iter " + Dist.EMIterationCount.ToString()
                + " T " + Dist.RunningPWC.Temperature.ToString("E4") + " PWHammy " + Dist.RunningPWC.PairwiseHammy.ToString("E4") + " PMLoop " + Dist.PMLoopUsed.ToString() + " " 
                + looptype + " Time " + PWCUtility.HPDuration.ToString("F0"));

            string nextline = "";
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                double tmp = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
                if (ClusterIndex != 0)
                    nextline += "* ";
                nextline +=  ClusterIndex.ToString() + " C " + Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1")
                    + " (Frz " + Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex].ToString("F6") + ") " + "[Wdth " + tmp.ToString("F4") + "] ";
                if (Program.ContinuousClustering)
                    nextline += "P " + Dist.RunningPWC.P_k_[ClusterIndex].ToString("F4") + " ";
            }
            PWCUtility.SALSAPrint(1, nextline);

        }   // End Print Iteration from getDist

        public void diagnosticsplitprint(bool ResultofSplittingTest, double[] Local_A_k_, double[] Local_C_k_, double[] Local_FreezingMeasure_k_, double[] Local_P_k_)
        {
            string nextline1 = "A) " + Dist.RunningPWC.Ncent.ToString() + " Iter " + Dist.EMIterationCount.ToString() + " T "
                + Dist.RunningPWC.Temperature.ToString("E4") + " PWHammy " + Dist.RunningPWC.PairwiseHammy.ToString("E4")
                + " PMLoop " + Dist.PMLoopUsed.ToString() + " Time " + PWCUtility.HPDuration.ToString("F0");
            if (Dist.RunningPWC.ClustertoSplit >= 0)
            {
                if (ResultofSplittingTest)
                    nextline1 += " C# To Split " + Dist.RunningPWC.ClustertoSplit.ToString();
                else
                    nextline1 += " No Split";
                nextline1 += " " + ResultofSplittingTest.ToString() + " Status " + vectorclass.eigenconverged[Dist.RunningPWC.ClustertoSplit].ToString()
                    + " Pass 1 " + vectorclass.Eigenvalues_Current[Dist.RunningPWC.ClustertoSplit].ToString("E4")
                    + " Pass 0 " + vectorclass.Eigenvalues_Pass0[Dist.RunningPWC.ClustertoSplit].ToString("E4");
            }
            PWCUtility.SALSAPrint(1, nextline1);
            nextline1 = "   ";
            double a_k_sum = 0.0;
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                a_k_sum += Local_A_k_[ClusterIndex];
            }
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                double tmp = Local_A_k_[ClusterIndex] / a_k_sum;
                string endofline = "]";
                if (Program.ContinuousClustering)
                    endofline = "] [P " + Local_P_k_[ClusterIndex].ToString("F4") + "]";
                if (ClusterIndex != 0)
                    nextline1 += " * ";
                nextline1 += ClusterIndex.ToString() + " " + Local_C_k_[ClusterIndex].ToString("F1") + " (Frz " + Local_FreezingMeasure_k_[ClusterIndex].ToString("F6")
                    + ") [Hfract" + tmp.ToString("F4") + endofline;
            }
            PWCUtility.SALSAPrint(1, nextline1);
            PWCUtility.SALSAPrint(1, " ");
            return;

        }   // End diagnosticsplitprint

        public void PairwiseThread()
        {
            Diff_Epsilon_k_ = new double[Dist.RunningPWC.Ncent];

            GlobalReductions.FindVectorDoubleSum Find_C_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
            GlobalReductions.FindVectorDoubleSum Find_FreezingMeasure_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
            GlobalReductions.FindDoubleSum Find_ChangeinM = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
            GlobalReductions.FindDoubleSum Find_NormalizechangeinM = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount); 

            int pmloopmax = 1;
            if (Program.ContinuousClustering && (Dist.needtocalculateMalpha_k_ == 1) && (Dist.RunningPWC.Ncent > 1) )
                pmloopmax = 40;
            double ChangeinM = 0.0;
            double NormalizechangeinM = 0.0;
            Dist.PMLoopUsed = 0;

            for (int pmloop = 0; pmloop < pmloopmax; pmloop++)
            {
                if (pmloop > 0)
                {
                    ChangeinM = 0.0;
                    Find_C_k_.zero();
                    Find_FreezingMeasure_k_.zero();
                    Find_ChangeinM.zero();
                }

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    //	Start Code setting Malpha_k_ and partialsum_C_k_
                    int localNcent = Dist.RunningPWC.Ncent;
                    double[] Accum_C_k_ = new double[localNcent];
                    double[] Accum_FreezingMeasure_k_ = new double[localNcent];
                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    {
                        Accum_C_k_[ClusterIndex] = 0.0;
                        Accum_FreezingMeasure_k_[ClusterIndex] = 0.0;
                    }

                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {   //	calculate Malpha_k_ if needed (controlled by Dist.needtocalculateMalpha_k_)
                        //                  if ((Dist.needtocalculateMalpha_k_ == -1) && (Dist.Ncent > 2)) Dist.needtocalculateMalpha_k_ = 1;
                        if (Dist.needtocalculateMalpha_k_ == -1)
                        {	// 	Initialize
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                //                          double fudge = 1.0 + Program.PerturbationFactor;
                                //                         if (PointIndex % 2 == 0) fudge = 2.0 - fudge;
                                Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = 1.0 / Dist.RunningPWC.Ncent;
                                //                          Malpha_k_[PointIndex][1] = 1.0 - 0.5 * fudge;
                            }
                        }

                        if (Dist.needtocalculateMalpha_k_ == 1)
                        {	// find Malpha_k_ from EM Method
                            double Minepsi = 0.0;
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                double tmpepsi = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                if (ClusterIndex == 0)
                                    Minepsi = tmpepsi;
                                else
                                    Minepsi = Math.Min(Minepsi, tmpepsi);
                            }
                            double tmp = 0.0;
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                double tmpepsi = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = Math.Exp(-(tmpepsi - Minepsi) / Dist.RunningPWC.Temperature); //exp(-epsi_a(k)/T)
                                if (Program.ContinuousClustering)
                                    Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] *= Dist.RunningPWC.P_k_[ClusterIndex];
                                tmp += Dist.RunningPWC.Weight_k_[ClusterIndex] * Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]; // tmp = Z_a = Sum_k(P(k)*exp(-epsi_a(k)/T))
                            }

                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            { // Normalize Malpha_k_[index][ClusterIndex] 
                                if (tmp != 0)
                                {
                                    Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] / tmp; // Malpha_k_[index][ClusterCenter] = M_a(k) = exp(-epsi_a(k)/T) / Z_a
                                }
                                else
                                {
                                    Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = 1.0 / localNcent; //  average equally over clusters 
                                }
                                if (pmloopmax > 1)
                                {
                                    if (pmloop == 0)
                                        Find_NormalizechangeinM.addapoint(ThreadNo, Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]);
                                    if (pmloop > 0)
                                        Find_ChangeinM.addapoint(ThreadNo, Math.Abs(Dist.RunningPWC.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] - Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]));
                                    Dist.RunningPWC.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex];
                                }
                            }

                        }	// End calculate Malpha_k_ from EM Method


                        for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                        {   // Find Cluster member counts C_k_ in each thread
                            double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex];
                            Accum_C_k_[ClusterIndex] += tmp;
                            Accum_FreezingMeasure_k_[ClusterIndex] += tmp * (1.0 - tmp);
                        }

                    }   // End loop over points 
                    Find_C_k_.addapoint(ThreadNo, Accum_C_k_);
                    Find_FreezingMeasure_k_.addapoint(ThreadNo, Accum_FreezingMeasure_k_);

                });	// End Code setting Malpha_k_ and partialsum_C_k_

                Find_C_k_.sumoverthreadsandmpi();
                Find_FreezingMeasure_k_.sumoverthreadsandmpi();
                if (pmloop > 0)
                {
                    Find_ChangeinM.sumoverthreadsandmpi();
                    ChangeinM = Find_ChangeinM.Total;
                }
                if ( (pmloop == 0) && (pmloopmax > 1))
                {
                    Find_NormalizechangeinM.sumoverthreadsandmpi();
                    NormalizechangeinM = Find_NormalizechangeinM.Total;
                }

                //  Form C_k_ and P_k_ from sum over threads
                for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
                {
                    Dist.RunningPWC.C_k_[ClusterCount] = Find_C_k_.TotalVectorSum[ClusterCount];
                    Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] = Find_FreezingMeasure_k_.TotalVectorSum[ClusterCount];
                }

                for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
                {
                    double tmp = Dist.RunningPWC.C_k_[ClusterCount];
                    if (tmp > 0.0)
                        Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] = Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] / tmp;
                    if (Program.ContinuousClustering)
                        Dist.RunningPWC.P_k_[ClusterCount] = Dist.RunningPWC.C_k_[ClusterCount] / PWCUtility.PointCount_Global;

                }
                //  See Shift in C after split
                if(Program.TestExpectedChange >= 0)
                {
                    double C1 = Dist.RunningPWC.C_k_[Program.TestExpectedChange];
                    double C2 = Dist.RunningPWC.C_k_[Dist.RunningPWC.Ncent - 1];
                    PWCUtility.SALSAPrint(0, "Cluster " + Program.TestExpectedChange.ToString() + " Method " + Program.ExpectedChangeMethod.ToString() + " Expected Change " + Program.ExpectedChange.ToString("E4")
                        + " Previous C " + Program.PreviousC.ToString("E4") + " New " + C1.ToString("E6") + " " + C2.ToString("E6") + " Temp " + Dist.RunningPWC.Temperature.ToString("e4") + " " + Dist.EMIterationCount.ToString());
                    Program.TestExpectedChange = -1;
                    if (Program.ExpectedChangeMethod == 1)
                    {
                        Program.deltaCoverC[0] += Math.Abs(C1 - C2) / Program.PreviousC;
                        Program.CountdeltaCoverC[0]++;
                    }
                    else
                    {
                        Program.deltaCoverC[1] += Math.Abs(C1 - C2) / Program.PreviousC;
                        Program.CountdeltaCoverC[1]++;
                    }


                }
                //  End iteration over M and p given epsilon
                if (pmloopmax > 1)
                {
                    Dist.PMLoopUsed = pmloop;
                    bool pmloopend = (pmloop > 0) && (ChangeinM < 0.001 * NormalizechangeinM);
                    PWCUtility.SynchronizeMPIvariable(ref pmloopend);
                    if (pmloopend)
                        break;
                    if (pmloop == pmloopmax - 1)
                        PWCUtility.SALSAPrint(1, " pmloop limit reached " + " Iter " + Dist.EMIterationCount.ToString()
                        + " T " + Dist.RunningPWC.Temperature.ToString("E4") + " PWHammy " + Dist.RunningPWC.PairwiseHammy.ToString("E4") + " " + ChangeinM.ToString("E4") + " " + NormalizechangeinM.ToString("E4"));
                }
            }   // End pmloop
            Program.CountPMExtraIterations += Dist.PMLoopUsed;

            calculateEpsi(Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Epsilonalpha_k_, Dist.RunningPWC.Ncent);
            Dist.needtocalculateMalpha_k_ = 1;	// as epsi are now set
            return;

        }	// End PairwiseThread 

        public int oldMaxLength = int.MinValue;
        public static MPIPacket<double> fromafar = null;
        public static MPIPacket<double> toafar = null;
        public static MPIPacket<double> myown = null;


        //	Perform multiple parallel steps calculating A_k_, Balpha_k_, epsi and differences
        public void calculateEpsi(double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, double[][] localepsi, int localNcent)
        {

            MPI.CompletedStatus MPIStatus;
            int sendtag = 0;
            int receivetag = 0;

            int Maxlength = localNcent * PWCUtility.PointCount_Largest;

            if (oldMaxLength != Maxlength)
            {
                fromafar = new MPIPacket<double>(Maxlength);
                toafar = new MPIPacket<double>(Maxlength);
                myown = new MPIPacket<double>(Maxlength);
                oldMaxLength = Maxlength;
            }
            else
            {
                fromafar.Clear();
                toafar.Clear();
                myown.Clear();
            }

            myown.FirstPoint = PWCUtility.PointStart_Process;
            myown.NumberofPoints = PWCUtility.PointCount_Process;
            if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
            {
                for (int countC = 0; countC < localNcent; countC++)
                {
                    myown.Marray[PWCUtility.PointCount_Process * localNcent + countC] = 0.0;
                }
            }
            int fromprocess = PWCUtility.MPI_Rank - 1;
            if (fromprocess < 0)
                fromprocess = PWCUtility.MPI_Size - 1;
            int toprocess = PWCUtility.MPI_Rank + 1;
            if (toprocess > PWCUtility.MPI_Size - 1)
                toprocess = 0;

            //	First communicationloop is local; then we have MPI_Size transfers of data in  a ring through processes             
            for (int MPICommunicationSteps = 0; MPICommunicationSteps < PWCUtility.MPI_Size; MPICommunicationSteps++)
            {
                if (MPICommunicationSteps == 1)
                {
                    toafar.FirstPoint = PWCUtility.PointStart_Process;
                    toafar.NumberofPoints = PWCUtility.PointCount_Process;
                    if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
                    {
                        for (int countC = 0; countC < localNcent; countC++)
                        {
                            int bigindex = PWCUtility.PointCount_Process * localNcent + countC;
                            toafar.Marray[bigindex] = myown.Marray[bigindex];
                        }
                    }
                }
                if (MPICommunicationSteps > 1)
                {
                    toafar.FirstPoint = fromafar.FirstPoint;
                    toafar.NumberofPoints = fromafar.NumberofPoints;
                    if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
                    {
                        for (int countC = 0; countC < localNcent; countC++)
                        {
                            int bigindex = PWCUtility.PointCount_Process * localNcent + countC;
                            toafar.Marray[bigindex] = fromafar.Marray[bigindex];
                        }
                    }
                }
                if (MPICommunicationSteps > 0)
                {
                    PWCUtility.StartSubTimer(PWCUtility.MPISENDRECEIVETiming);
                    PWCUtility.MPI_communicator.SendReceive<MPIPacket<double>>(toafar, toprocess, sendtag, fromprocess, receivetag, out fromafar, out MPIStatus);
                    PWCUtility.StopSubTimer(PWCUtility.MPISENDRECEIVETiming);
                }

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (item) =>
                {

                    //  Loop over Home (point) indices in thread using non-home values from MPI or locally for communicationloop == 0
                    int betastart, betatotal;
                    int indexlen = PWCUtility.PointsperThread[item];
                    int beginpoint = PWCUtility.StartPointperThread[item] - PWCUtility.PointStart_Process;

                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        if (MPICommunicationSteps == 0)
                        {   // Set up Home-Home Interaction
                            Array.Clear(localBalpha_k_[ProcessPointIndex], 0, localNcent);
                            betatotal = PWCUtility.PointCount_Process;
                            betastart = PWCUtility.PointStart_Process;
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                double tmp = localMalpha_k_[ProcessPointIndex][ClusterIndex];
                                int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                myown.Marray[bigindex] = tmp;
                                toafar.Marray[bigindex] = tmp;
                            }
                        }
                        else
                        {   // Set up Home-Away Interaction
                            betatotal = fromafar.NumberofPoints;
                            betastart = fromafar.FirstPoint;
                            if (MPICommunicationSteps != (PWCUtility.MPI_Size - 1))
                            {
                                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                                {
                                    int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                    toafar.Marray[bigindex] = fromafar.Marray[bigindex];
                                }
                            }
                        }


                        for (int betalocal = 0; betalocal < betatotal; betalocal++)  // Loop over Away indices
                        {
                            double tmp;

                            int betafull = betastart + betalocal;

                            double dijforthiscase = PWCParallelism.getDistanceValue(ProcessPointIndex + PWCUtility.PointStart_Process, betafull);

                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {

                                if (MPICommunicationSteps == 0)
                                {
                                    tmp = localMalpha_k_[betalocal][ClusterIndex];
                                }
                                else
                                {
                                    tmp = fromafar.Marray[betalocal * localNcent + ClusterIndex];
                                }

                                localBalpha_k_[ProcessPointIndex][ClusterIndex] += dijforthiscase * tmp / localC_k_[ClusterIndex];
                            }
                        }    // End Loop over Away indices

                    }	//  Loop over Home (point) indices in thread calculating all Balpha_k_
                });     // End delegate threads parallelizing home indices
            }	// End loop over communicationloop

            //  Now calculate quantities involving global sums

            //	Calculate full A(k)
            GlobalReductions.FindVectorDoubleSum Find_A_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] LocalContribution_A_k_ = new double[localNcent];
                
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    {
                        LocalContribution_A_k_[ClusterIndex] = localBalpha_k_[ProcessPointIndex][ClusterIndex] * localMalpha_k_[ProcessPointIndex][ClusterIndex];

                    }
                    Find_A_k_.addapoint(ThreadNo, LocalContribution_A_k_);
                }
            });

            Find_A_k_.sumoverthreadsandmpi();

            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                localA_k_[ClusterIndex] = -0.5 * Find_A_k_.TotalVectorSum[ClusterIndex] / localC_k_[ClusterIndex];

            // Calculate new values of epsi and do partial sums of differences   
            GlobalReductions.FindVectorDoubleSum Find_EpsiDiff = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] Local_EpsiDiff = new double[localNcent];

                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    {
                        double tmp = localBalpha_k_[ProcessPointIndex][ClusterIndex] + localA_k_[ClusterIndex];
                        localepsi[ProcessPointIndex][ClusterIndex] = tmp;
                        if (Dist.oldepsiset > 0) //the value will be used in doing merge test.
                            Local_EpsiDiff[ClusterIndex] = Math.Abs(Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] - tmp);

                        Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = tmp;
                    }
                    if (Dist.oldepsiset > 0)
                        Find_EpsiDiff.addapoint(ThreadNo, Local_EpsiDiff);
                }
            });

            if (Dist.oldepsiset > 0)
            {	//	Calculate epsidiff which is sum for each center over data points
                Find_EpsiDiff.sumoverthreadsandmpi();

                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    Diff_Epsilon_k_[ClusterIndex] = Find_EpsiDiff.TotalVectorSum[ClusterIndex];
  
            }	// End computation of epsidiff
            ++Dist.oldepsiset;

        }	// End calculateEpsi
        

        // The change of sum(points) Delta(epsi)/# Points should be less than argument ChangeLimit for each cluster
        // epsiDiff calculated in CalculateEpsi while Epsi are updated
        //  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
        public static int convergenceTest(double ChangeLimit)
        {

            Dist.HitConvergenceLoopLimit = false;
            if (Dist.oldepsiset <= 1)
                return 0; // Skip first two cases -- WHY?
            double epsidiffmax = 0;
            for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
            {
                if (Dist.Diff_Epsilon_k_[ClusterCount] > epsidiffmax)
                    epsidiffmax = Dist.Diff_Epsilon_k_[ClusterCount];
            }
            epsidiffmax = epsidiffmax / PWCUtility.PointCount_Global;
            bool epslimit = (epsidiffmax <= ChangeLimit);
            PWCUtility.SynchronizeMPIvariable(ref epslimit);

            if (epslimit)
            {
                Dist.EMIterationStepCount = 0;
                if (Dist.OnLastleg || Dist.justconverging)
                    return 2;
                bool templimit = (Dist.RunningPWC.Temperature < Dist.Tmin);
                PWCUtility.SynchronizeMPIvariable(ref templimit);
                if ((Dist.RunningPWC.Ncent >= Dist.ActualMaxNcent) || templimit || (Dist.SplitFailures > 1))
                {
                    Dist.OnLastleg = true;
                    return 2;
                }
                else
                    return 1;

            }
            if (Dist.EMIterationStepCount > Program.ConvergenceLoopLimit)
            {
                PWCUtility.SALSAPrint(0, "EM Convergence Loop Reached for fixed Temperature at Iteration " + Dist.EMIterationCount.ToString() );
                Dist.EMIterationStepCount = 0;
                Dist.HitConvergenceLoopLimit = true;
                return 2;
            }
            else
                return 0;

        }   // End convergenceTest

        public void SaveCurrentTask()
        {   // Save current task that should be split but we need to converge first

            ClusteringSolution.SetAxinSolution(Dist.RunningPWC); 
            ClusteringSolution.CopySolution(Dist.RunningPWC, Dist.SavedPWC);
            return;

        }   // End saving current task that should be split but we need to converge first

        public void RestorePreviousTask()
        {   // Restore previous task that should be split but we needed to converge current cluster configuration first

            ClusteringSolution.CopySolution(Dist.SavedPWC, Dist.RunningPWC);
            ClusteringSolution.RestoreAxfromSolution(Dist.RunningPWC);
            Dist.RunningPWC.Axset = false;
            return;

        }   // End Restore previous task that should be split but we needed to converge current cluster configuration first

        public void JiggleClusters()
        {   // Jiggle Values of Epsilon -- leaving other Cluster Parameters


            if (Dist.RunningPWC.Ncent == 1)
                return;

            // Initialize
            Random Randobject = new Random();
            vectorclass vcjiggle = new vectorclass();

            if (Program.JigglePerturbation == 2)
            {
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                {
                    Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
                    Dist.ClusterSelected[ClusterIndex] = 1;
                }
                Dist.RunningPWC.ClustertoSplit = 1;
                Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;

                // Note this calculates Full eigenvector of entire matrix -- NOT as in stability analysis, the eigenvector in each cluster sector
                vcjiggle.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
                double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

                PWCUtility.SALSAPrint(1, " Jiggle Cluster " + Dist.PlaceforEigsforMultipleCluster.ToString()
                        + " Status " + vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster].ToString() + " Pass 1 "
                        + vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster].ToString("E4")
                        + " Pass 0 " + vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster].ToString("E4"));

                if ((vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] <= 0.0) || (vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster] <= 0))
                    return;
            }   // End Initialization of JiggleOption = 2

            // NewC_k_ and AverageMalpha_k_Change calculated for output only
            GlobalReductions.FindVectorDoubleSum Find_NewC_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
            GlobalReductions.FindDoubleSum Find_AverageMalpha_k_Change = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] NewMalpha_k_ = new double[Dist.RunningPWC.Ncent];
                double[] partialsum_NewC_k_ = new double[Dist.RunningPWC.Ncent];

                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    double AverageMalpha_k_Change = 0.0;
                    double tmp = 0.0;
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        double perturb;
                        if (Program.JigglePerturbation == 1)
                        { // Perturb Epsilon randomly
                            perturb = Randobject.NextDouble();
                        }
                        else
                        {   // Perturb Epsilon with nastiest eigenvector
                            perturb = vectorclass.oldAx[ProcessPointIndex][ClusterIndex];
                        }
                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] += perturb * Program.JigglePerturbationFactor * Dist.RunningPWC.Temperature;

                        NewMalpha_k_[ClusterIndex] = Math.Exp(-Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] / Dist.RunningPWC.Temperature);
                        tmp += NewMalpha_k_[ClusterIndex];
                    }
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        double addtoC = NewMalpha_k_[ClusterIndex] / tmp;
                        partialsum_NewC_k_[ClusterIndex] += addtoC;
                        AverageMalpha_k_Change += Math.Abs(addtoC - Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]);
                    }
                    Find_NewC_k_.addapoint(ThreadNo, partialsum_NewC_k_);
                    Find_AverageMalpha_k_Change.addapoint(ThreadNo, AverageMalpha_k_Change);
                }
            });

            Find_AverageMalpha_k_Change.sumoverthreadsandmpi();
            double FullAverageMalpha_k_Change = Find_AverageMalpha_k_Change.Total / (Dist.RunningPWC.Ncent * PWCUtility.PointCount_Global);
            Dist.needtocalculateMalpha_k_ = 1;

            double[] NewC_k_ = new double[Dist.RunningPWC.Ncent];
            Find_NewC_k_.sumoverthreadsandmpi();
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                NewC_k_[ClusterIndex] = Find_NewC_k_.TotalVectorSum[ClusterIndex];

            string nextline = " Jiggle M change " + FullAverageMalpha_k_Change.ToString("F5") + "  Old(Jiggled) Sizes ";
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                nextline += Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1") + " (" + NewC_k_[ClusterIndex].ToString("F6") + ") ";
            }
            PWCUtility.SALSAPrint(1, nextline);
            return;

        }   // End JiggleCluster()

        // Decide if to split based on negative eigenvalue of second derivative matrix
        // MinimumEigenvalue is MINIMUM eigenvalue
        // ClustertoSplit is cluster number of this
        // D is distance matrix
        //  Skip first time this happens as initialize to two clusters
        public bool shouldweSplit()
        {
            //	Calculate Cluster with minimum eigenvalue -- find cluster number and eigenvalue (which could be positive)
            Dist.Numberthatcanbesplit = 0;
            int LimitonSplits = Math.Min(Program.MaxNumberSplitClusters, Dist.ActualMaxNcent - Dist.RunningPWC.Ncent);

            Dist.RunningPWC.ClustertoSplit = -1;
            Dist.RunningPWC.MinimumEigenvalue = 99999999999.0;
            int ActualMethodology = Program.Eigenvalue_Methodology;
            if (Dist.RunningPWC.Ncent == 1)
                ActualMethodology = 2;
            if (ActualMethodology == 0)
                return false;

            vectorclass vc = new vectorclass();
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                vectorclass.eigenconverged[ClusterIndex] = 0;
                Dist.ClusterSelected[ClusterIndex] = 1;
            }

            // Continuous Clustering is simple case of ActualMethodology  = 1
            if (ActualMethodology <= 2)
            {
                vc.getEigenvalues(ActualMethodology, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
                for (int ClusterToRefine = 0; ClusterToRefine < Dist.RunningPWC.Ncent; ClusterToRefine++)
                {
                    if (Dist.RunningPWC.C_k_[ClusterToRefine] <= Program.ToosmalltoSplit)
                        continue;
                    double TrueMinimumEigenvalue = vectorclass.Eigenvalues_Pass0[ClusterToRefine] - vectorclass.Eigenvalues_Current[ClusterToRefine];
                    PWCUtility.SALSAPrint(1, "       Cluster " + ClusterToRefine.ToString() + " Status "
                        + vectorclass.eigenconverged[ClusterToRefine].ToString() + " Pass 1 " + vectorclass.Eigenvalues_Current[ClusterToRefine].ToString("E4")
                        + " Pass 0 " + vectorclass.Eigenvalues_Pass0[ClusterToRefine].ToString("E4") + " Diff " + TrueMinimumEigenvalue.ToString("E2") 
                        + " Size " + Dist.RunningPWC.C_k_[ClusterToRefine].ToString("F1"));

                    if ((vectorclass.Eigenvalues_Pass0[ClusterToRefine] > 0.0) || (vectorclass.eigenconverged[ClusterToRefine] > 0))
                    {
                        if (TrueMinimumEigenvalue < Dist.RunningPWC.MinimumEigenvalue)
                        {
                            Dist.RunningPWC.MinimumEigenvalue = TrueMinimumEigenvalue;
                            Dist.RunningPWC.ClustertoSplit = ClusterToRefine;
                        }
                        if(LimitonSplits > 1) {
                            double eigtest1 = Program.MinEigtest * vectorclass.Eigenvalues_Pass0[ClusterToRefine];
                            if (TrueMinimumEigenvalue < eigtest1)
                            {   // Candidate for Split List
                                if (Dist.Numberthatcanbesplit == 0)     // Initialize Split List
                                {
                                    Dist.Numberthatcanbesplit = 1;
                                    Dist.EigsofClusterstoSplit[0] = TrueMinimumEigenvalue;
                                    Dist.ListofClusterstoSplit[0] = ClusterToRefine;
                                }
                                else     // Add to Split List
                                {
                                    int position = Dist.Numberthatcanbesplit;
                                    for (int positionloop = 0; positionloop < Dist.Numberthatcanbesplit; positionloop++)
                                    {
                                        if (TrueMinimumEigenvalue >= Dist.EigsofClusterstoSplit[positionloop])
                                            continue;
                                        position = positionloop;
                                        break;
                                    }
                                    if (position >= LimitonSplits)
                                        continue;
                                    for (int positionloop = Dist.Numberthatcanbesplit - 1; positionloop >= position; positionloop--)
                                    {
                                        if (positionloop == (LimitonSplits - 1))
                                            continue;
                                        Dist.EigsofClusterstoSplit[positionloop + 1] = Dist.EigsofClusterstoSplit[positionloop];
                                        Dist.ListofClusterstoSplit[positionloop + 1] = Dist.ListofClusterstoSplit[positionloop];
                                    }
                                    Dist.Numberthatcanbesplit = Math.Min(Dist.Numberthatcanbesplit + 1, LimitonSplits);
                                    Dist.EigsofClusterstoSplit[position] = TrueMinimumEigenvalue;
                                    Dist.ListofClusterstoSplit[position] = ClusterToRefine;
                                }
                            }
                        }
                    }
                }

                if (Program.PerformEigenTest && (Dist.RunningPWC.Ncent == 1))
                {
                    // Eigenvalue Test for Metholodology 2 -- only sensible if Ncent = 1 as otherwise Method 2 inaccurate so test will fail
                    Dist.RunningPWC.ClustertoSplit = 0;

                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {	//	Start Initialization Code
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                        for (int index = beginpoint; index < indexlen + beginpoint; index++)
                        {
                            Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.ClustertoSplit];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[index][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.ClustertoSplit];
                        }
                    });	// End 


                    Dist.RunningPWC.Ncent++;
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
                        Dist.ClusterSelected[ClusterIndex] = 0;
                    }
                    Dist.ClusterSelected[0] = 1;
                    Dist.ClusterSelected[1] = 1;
                    Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;
                    PairwiseThread();
                    vc.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
                    double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

                    string nextline = "  C)Test Cluster " + Dist.RunningPWC.ClustertoSplit + " T " + Dist.RunningPWC.Temperature.ToString("E4") + " PMLoop " + Dist.PMLoopUsed.ToString() +
                        " PWHammy " + Dist.RunningPWC.PairwiseHammy.ToString("E4") + " Status " + vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster].ToString()
                        + " Status " + vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster].ToString()
                        + " Pass 1 " + vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster].ToString("E4")
                        + " Pass 0 " + vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster].ToString("E4") + " C ";
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        nextline += Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1") + " ";
                    }
                    PWCUtility.SALSAPrint(1, nextline);

                    Dist.RunningPWC.Ncent--;
                }   // End Test of Eigenvalues for Methodology 2

            }   // End ActualMethodology <= 2

            else
            {
                // Implement exact Duplication separately for each cluster ActualMethodology = 3 or 4
                // Parallel reSetting of Epsilon
                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {	//	Start Initialization Code
                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                    for (int index = beginpoint; index < indexlen + beginpoint; index++)
                    {
                        for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                        {
                            Dist.RunningPWC.Master_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                            Dist.RunningPWC.Best_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                        }
                    }

                });	// End Initialization Code for Epsilon

                double[] Save_C_k_ = new double[Dist.RunningPWC.Ncent];
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    Save_C_k_[ClusterIndex] = Dist.RunningPWC.C_k_[ClusterIndex];

                for (int ClusterToRefine = 0; ClusterToRefine < Dist.RunningPWC.Ncent; ClusterToRefine++)
                {
                    if (Save_C_k_[ClusterToRefine] <= Program.ToosmalltoSplit)
                        continue;
                    //  Do EM calculation
                    Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached
                    Dist.needtocalculateMalpha_k_ = 1; // =0 Already Set; = 1 Calculate from EM (usual); = -1 Initialize

                    //	Loop over EM calculations
                    for (int EMLoop = 0; EMLoop < Program.EMlimit_Duplication; EMLoop++)
                    {
                        for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                        {
                            Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
                        }
                        Dist.RunningPWC.Weight_k_[ClusterToRefine] = 2.0;
                        PairwiseThread();
                        Extra_EMIterationCount++;

                        //	Now see if we are done -- take results from Rank 0
                        int convergence = 0;
                        convergence = Dist.convergenceTest(Program.Epsi_max_change_Duplication);
                        PWCUtility.SynchronizeMPIvariable(ref convergence);

                        if ( (convergence > 0) || (EMLoop == Program.EMlimit_Duplication - 1) )
                        {
                            string nextline = "     D)Special Iter " + EMLoop.ToString() + " T " + Dist.RunningPWC.Temperature.ToString("E4") + " PMLoop " + Dist.PMLoopUsed.ToString()
                                + " PWHammy " + Dist.RunningPWC.PairwiseHammy.ToString("E4") + " Cluster " + ClusterToRefine.ToString() + " Sizes ";
                            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                            {
                                double tmp = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
                                nextline += Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1") + " ("
                                    + Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex].ToString("F6") + ") [" + tmp.ToString("F1") + "] ";
                            }
                            PWCUtility.SALSAPrint(1, nextline);
                        }
                        if (convergence > 0)
                            break;

                    }   // End Loop over EMLoop

                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        vectorclass.eigenconverged[ClusterIndex] = 0;
                        Dist.ClusterSelected[ClusterIndex] = 0;
                    }
                    Dist.ClusterSelected[ClusterToRefine] = 1;
                    vc.getEigenvalues(3, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
                    double TrueMinimumEigenvalue = vectorclass.Eigenvalues_Pass0[ClusterToRefine] - vectorclass.Eigenvalues_Current[ClusterToRefine];

                    PWCUtility.SALSAPrint(1, "       Cluster " + ClusterToRefine.ToString() + " Status "
                        + vectorclass.eigenconverged[ClusterToRefine].ToString() + " Pass 1 " + vectorclass.Eigenvalues_Current[ClusterToRefine].ToString("E4")
                        + " Pass 0 " + vectorclass.Eigenvalues_Pass0[ClusterToRefine].ToString("E4"));

                    if ((vectorclass.Eigenvalues_Pass0[ClusterToRefine] > 0.0) && (vectorclass.eigenconverged[ClusterToRefine] > 0))
                    {
                        if (TrueMinimumEigenvalue < Dist.RunningPWC.MinimumEigenvalue)
                        {
                            Dist.RunningPWC.MinimumEigenvalue = TrueMinimumEigenvalue;
                            Dist.RunningPWC.ClustertoSplit = ClusterToRefine;
                        }
                    }

                    // Parallel Saving of Epsilon
                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {	//	Start Initialization Code
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                            {
                                if (Dist.RunningPWC.ClustertoSplit == ClusterToRefine)
                                    Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                if (ClusterToRefine < Dist.RunningPWC.Ncent - 1)
                                {
                                    Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                    Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                }
                                else
                                {
                                    Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                    Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                }
                            }

                        }
                    });	// End 
                }   // End case Methodology 3 or 4

                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;

                if (Program.PerformEigenTest)
                {
                    // Eigenvalue Test for Methodology 3
                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {	//	Start Initialization Code
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                        }
                    });	// End 

                    Dist.RunningPWC.Ncent++;
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
                        Dist.ClusterSelected[ClusterIndex] = 0;
                    }
                    Dist.ClusterSelected[Dist.RunningPWC.ClustertoSplit] = 1;
                    Dist.ClusterSelected[Dist.RunningPWC.Ncent - 1] = 1;

                    Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;
                    PairwiseThread();
                    vc.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
                    double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

                    string nextline = " C)Test Cluster " + Dist.RunningPWC.ClustertoSplit + " T " + Dist.RunningPWC.Temperature.ToString("E4") + " PMLoop " + Dist.PMLoopUsed.ToString() +" PWHammy "
                        + Dist.RunningPWC.PairwiseHammy.ToString("E4") + Dist.PlaceforEigsforMultipleCluster.ToString() + " Status "
                        + vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster].ToString() + " Pass 1 "
                        + vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster].ToString("E4") + " Pass 0 "
                        + vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster].ToString("E4") + " C ";
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        nextline += Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1") + " ";
                    }
                    PWCUtility.SALSAPrint(1, nextline);

                    Dist.RunningPWC.Ncent--;
                }   // End Eigenvalue Test for Metholodology 3

            }   // End case Methodology 3 or 4

            if (Dist.RunningPWC.ClustertoSplit < 0)
            {
                Dist.SplitFailures++;
                return false;
            }
            double eigtest = Program.MinEigtest * vectorclass.Eigenvalues_Pass0[Dist.RunningPWC.ClustertoSplit];
            if (Dist.RunningPWC.MinimumEigenvalue > eigtest)
            {
                return false;
            }
            return true;

        }   // End shouldweSplit

        //  Do a split on the identified cluster
        //  Malpha_k_ is value of M indexed by (point,cluster)
        //  New clusters stored in last cluster and old position
        //  Cluster count is incremented by one
        //  ClustertoSplit is cluster number to split
        //	Embarassingly Parallel over Processes -- no MPI calls needed except in normalization
        public void dothesplit()
        {
            double PerturbationNormFactor = 1.0;

            //  Calculate Normalization
            if (Program.PerturbationVehicle == 0)
            {
                GlobalReductions.FindDoubleSum SumoverShifts = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
                GlobalReductions.FindDoubleSum ShiftNorm = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
                GlobalReductions.FindDoubleSum EpsNorm = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                        tmp = tmp * vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] / Dist.RunningPWC.Temperature;
                        if (Dist.RunningPWC.Ncent == 1)
                            tmp = Math.Abs(tmp);
                        double tmp1 = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor;
                        double tmp2 = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                        SumoverShifts.addapoint(ThreadNo, tmp);
                        ShiftNorm.addapoint(ThreadNo, tmp1 * tmp1);
                        EpsNorm.addapoint(ThreadNo, tmp2 * tmp2);
                    }
                });
                SumoverShifts.sumoverthreadsandmpi();
                ShiftNorm.sumoverthreadsandmpi();
                EpsNorm.sumoverthreadsandmpi();

                Program.TestExpectedChange = Dist.RunningPWC.ClustertoSplit;
                Program.ExpectedChangeMethod = 0;
                if (Dist.RunningPWC.Ncent == 1)
                    Program.ExpectedChangeMethod = 1;

                double CShiftNorm = SumoverShifts.Total;
                CShiftNorm = Math.Abs (0.05 * Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit] / CShiftNorm);
                PerturbationNormFactor = 0.1 * Math.Sqrt(EpsNorm.Total / ShiftNorm.Total);
                PWCUtility.SALSAPrint(0, "Perturb Cluster " + Program.TestExpectedChange.ToString() + " Sum Shift " + SumoverShifts.Total.ToString("E4") + " Shift Norm " + ShiftNorm.Total.ToString("E4") + " Eps Norm " + EpsNorm.Total.ToString("E4")
                     + " CShift " + CShiftNorm.ToString("E4") + " Eps Shift " + PerturbationNormFactor.ToString("E4") + " Iter " + Dist.EMIterationCount.ToString());

                bool testchange = CShiftNorm < PerturbationNormFactor;
                PWCUtility.SynchronizeMPIvariable(ref testchange);
                Program.ExpectedChange = 0.05 * Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor;
                if (testchange)
                    PerturbationNormFactor = CShiftNorm;
                else
                {
                    if( Program.ExpectedChangeMethod == 0)
                        Program.ExpectedChangeMethod = 2;
                    Program.ExpectedChange *= PerturbationNormFactor / CShiftNorm;
                }
                if (Program.ExpectedChangeMethod == 0)
                    ++Program.Countmethodzero;
                if (Program.ExpectedChangeMethod == 2)
                    ++Program.Countmethodtwo;

                Program.PreviousC = Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit]*0.5;
                if (Program.ExpectedChangeMethod != 1)
                {
                    Program.deltaCoverC[4] += 2.0 * Math.Abs(Program.ExpectedChange) / Program.PreviousC;
                    Program.CountdeltaCoverC[4]++;
                }

                if (Program.ContinuousClustering)
                {
                    GlobalReductions.FindDoubleSum NewC1 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
                    GlobalReductions.FindDoubleSum NewC2 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                            double perturb = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor * PerturbationNormFactor;

                            double tmp2 = tmp * Math.Exp(perturb / Dist.RunningPWC.Temperature);
                            double tmp1 = tmp * Math.Exp(-perturb / Dist.RunningPWC.Temperature);
                            double NewBottom = 1.0 - 2.0 * tmp + tmp1 + tmp2;
                            NewC1.addapoint(ThreadNo, tmp1 / NewBottom);
                            NewC2.addapoint(ThreadNo, tmp2 / NewBottom);
                        }
                    });
                    NewC1.sumoverthreadsandmpi();
                    NewC2.sumoverthreadsandmpi();
                    PWCUtility.SALSAPrint(0, "Real Expectation " + NewC1.Total.ToString("E6") + " " + NewC2.Total.ToString("E6") + " Temp " + Dist.RunningPWC.Temperature.ToString("e4"));
                    if (Program.ExpectedChangeMethod == 1)
                    {
                        Program.deltaCoverC[2] += Math.Abs(NewC1.Total - NewC2.Total) / Program.PreviousC;
                        Program.CountdeltaCoverC[2]++;
                    }
                    else
                    {
                        Program.deltaCoverC[3] += Math.Abs(NewC1.Total - NewC2.Total) / Program.PreviousC;
                        Program.CountdeltaCoverC[3]++;
                    }
                }
            }
            
            Dist.RunningPWC.Ncent++;

            //  Parallel Section Splitting Cluster with delegate for action reading thread # from StartindexPort
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    if (Program.PerturbationVehicle == 1)
                    { // Perturb Malpha_k_
                        double newvalueofM = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                        double fudge = 1.0 + Program.SplitPerturbationFactor;
                        if (ProcessPointIndex % 2 == 0)
                            fudge = 2.0 - fudge;
                        Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] = newvalueofM * fudge; //split in positive direction
                        Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent - 1] = newvalueofM * (2.0 - fudge); //split in negative direction
                    }
                    else
                    { // Perturb Epsilon where we must also set P_k_ as Malpha_k_ formula involves P_k_
                        double perturb = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor * PerturbationNormFactor;
                        double original = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];

                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] = original + perturb;
                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent - 1] = original - perturb;
                    }
                }
            }); //  End Parallel Section Splitting Cluster

            if (Program.PerturbationVehicle == 1)
                Dist.needtocalculateMalpha_k_ = 0;
            else
            {
                // If continuous Clustering set P_k_ as Malpha_k_ formula involves P_k_
                // No need to perturb P_k_ -- they are just halved from last time around
                Dist.needtocalculateMalpha_k_ = 1;
                if (Program.ContinuousClustering)
                {
                    Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit] = 0.5 * Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit];
                    Dist.RunningPWC.P_k_[Dist.RunningPWC.Ncent - 1] = Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit];
                }
            }
                
            oldepsiset = 0;

        }   // End dothesplit

        //	Find initial Temperature
        public static void initializeTemperature()
        {

            GlobalReductions.FindDoubleSum Find_avg1 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
            GlobalReductions.FindDoubleSum Find_avg2 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
            GlobalReductions.FindDoubleSum Find_avg3 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadIndex) =>
            {
                double DistanceSum = 0.0;
                double NumberSum = 0.0;
                double STDSum = 0.0;
                int indexlen = PWCUtility.PointsperThread[ThreadIndex];
                int beginpoint = PWCUtility.StartPointperThread[ThreadIndex] - PWCUtility.PointStart_Process;

                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    int GlobalIndex = ProcessPointIndex + PWCUtility.PointStart_Process;

                    for (int PointIndex1 = 0; PointIndex1 < PWCUtility.PointCount_Global; PointIndex1++)
                    {
                        if (PointIndex1 == GlobalIndex)
                            continue;

                        double placevalue = PWCParallelism.getDistanceValue(GlobalIndex, PointIndex1);
                        DistanceSum += placevalue;
                        STDSum += placevalue * placevalue;
                        NumberSum += 1.0;
                    }
                }

                Find_avg1.addapoint(ThreadIndex, DistanceSum);
                Find_avg2.addapoint(ThreadIndex, NumberSum);
                Find_avg3.addapoint(ThreadIndex, STDSum);

            });  // End Parallel Section calculating average distance

            Find_avg1.sumoverthreadsandmpi();
            Find_avg2.sumoverthreadsandmpi();
            Find_avg3.sumoverthreadsandmpi();
            
            // Calculate global averages
            double[] avg_global = new double[3];
            
            avg_global[0] = Find_avg1.Total;
            avg_global[1] = Find_avg2.Total;
            avg_global[2] = Find_avg3.Total;

            //  Estimate of Initial Dist.Temperature is average distance
            //  Fudge factor of 2 over estimated critical temperature for first split
            double DistceMean = avg_global[0] / avg_global[1];
            double DistceSTD = Math.Sqrt((avg_global[2] / avg_global[1]) - DistceMean * DistceMean);
            double EstimatedDimension = 2.0 * DistceMean * DistceMean / (DistceSTD * DistceSTD);
            double IndividualSigma = Math.Sqrt(DistceMean / EstimatedDimension);
            Dist.Tinitial = 2.0 * DistceMean;

            PWCUtility.SALSAPrint(1, "Initial Temperature " + Dist.Tinitial.ToString("F4") + " # Distces " + avg_global[1].ToString()
                + " Mean " + DistceMean.ToString("F4") + " STD " + DistceSTD.ToString("F4")
                + " Estimated Dimension " + EstimatedDimension.ToString("F3") + " IndividualSigma " + IndividualSigma.ToString("F4"));
            return;

        }	// End Initialize  Dist.InitializeTemperature


        // Calculate Distance Correlation Function
        public void CorrelationCalculation()
        {
            int localNcent = Dist.RunningPWC.Ncent;
            Dist.DistanceCorrelation = new double[localNcent, localNcent];

// Summation of one row of Correlation Matrix
            GlobalReductions.FindVectorDoubleSum Find_CorrelationRow = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

            // Loop over rows over Correlation Matrix
            for (int CorrelationrowIndex = 0; CorrelationrowIndex < Dist.RunningPWC.Ncent; CorrelationrowIndex++)
            {
                if(CorrelationrowIndex != 0)
                    Find_CorrelationRow.zero();

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    //	Start Code setting partialsum_Correlation
                    double[] TempCorrel = new double[localNcent];

                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            TempCorrel[ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] * Dist.RunningPWC.Balpha_k_[ProcessPointIndex][CorrelationrowIndex];

                        Find_CorrelationRow.addapoint(ThreadNo, TempCorrel);

                    }   // End loop over points in this thread
                });	// End Code setting  partialsum_Correlation

                //  Form Correlation Row from sum over threads
                Find_CorrelationRow.sumoverthreadsandmpi();
                for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    Dist.DistanceCorrelation[CorrelationrowIndex, ClusterIndex] = Find_CorrelationRow.TotalVectorSum[ClusterIndex] / Dist.RunningPWC.C_k_[ClusterIndex];

            }   // End Loop over Rows of Correlation Matrix

            return;

        }	// End Correlation Calculation

        // Return False if current final solution has a problem
        //  Return True if there is no problem or no way of fixing problem
        public static bool CheckValidSolution()
        {
            bool[] RemoveCluster = new bool[Program.maxNcent];
            int Numbertoosmall = 0;
            Double ProbabilitySum = 0.0;
            string documentit = "";
            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
            {
                RemoveCluster[ClusterIndex] = false;
                if(Dist.RunningPWC.C_k_[ClusterIndex] < (double) Program.minimumclustercount - 0.5)
                    RemoveCluster[ClusterIndex] = true;
                PWCUtility.SynchronizeMPIvariable(ref RemoveCluster[ClusterIndex]);
                if (RemoveCluster[ClusterIndex])
                {
                    ++Numbertoosmall;
                    documentit += " " + ClusterIndex.ToString() + " # " + Dist.RunningPWC.C_k_[ClusterIndex].ToString("F1");
                    continue;
                }
                ProbabilitySum += Dist.RunningPWC.P_k_[ClusterIndex];
            }
            if (Numbertoosmall == 0)
                return true;

            //  Need to change solution as one or more clusters too small
            PWCUtility.SALSAPrint(1, Numbertoosmall.ToString() + " Clusters Too Small " + documentit);
            if (!Dist.BestPWC.SolutionSet)
            {
                PWCUtility.SALSAPrint(1, "Solution Not changed as no best solution");
                return true;
            }
            if (Dist.BestPWC.Ncent != Dist.RunningPWC.Ncent)
            {
                PWCUtility.SALSAPrint(1, " Best Solution Not Used to restart as Number of Centers " + Dist.BestPWC.Ncent + " Different from Running Solution with " + Dist.RunningPWC.Ncent);
                return true;
            }

            ClusteringSolution.CopySolution(Dist.BestPWC, Dist.RunningPWC);
            Dist.RunningPWC.OldHammy = 0.0;
            Dist.RunningPWC.PairwiseHammy = 0.0;
            Dist.RunningPWC.Axset = false;
            Dist.RunningPWC.ClustertoSplit = -1;

            for (int clusterindex = 0; clusterindex < Dist.RunningPWC.Ncent; clusterindex++)
            {
                if (!RemoveCluster[clusterindex])
                {
                    if (Program.ContinuousClustering)
                        Dist.RunningPWC.P_k_[clusterindex] = Dist.RunningPWC.P_k_[clusterindex]/ProbabilitySum;
                    continue;
                }
            }
            int oldNcent = Dist.RunningPWC.Ncent;
            for (int clusterindex = oldNcent-1; clusterindex >= 0; clusterindex--)
            {
                if (RemoveCluster[clusterindex])
                    ClusteringSolution.RemoveCluster(RunningPWC, clusterindex); // This both changes cluster count and shifts up clusters
            }
            Dist.ActualMaxNcent = Dist.RunningPWC.Ncent;
            return false;

        }   // End CheckValidSolution

    }	// End class dist

}   // End namespace cluster
