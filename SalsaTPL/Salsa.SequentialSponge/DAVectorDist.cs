using System;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using MPI;

using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    //  Generic Parallel Clustering
    class ParallelClustering
    {
        public static ClusteringSolution RunningSolution;    // The main Solution
        public static ClusteringSolution SavedSolution;      // Saved Solution
        public static ClusteringSolution BestSolution;       // Current Best Solution

    }   //  End Parallel Clustering

    // Control Kmeans
    class ControlKmeans
    {
        public static bool initialized = false;

        public ControlKmeans()
        {

            //allocate memory on first and indeed only call
            if (!initialized)
            {
                if (!Program.DoKmeans)
                {
                    Exception e = DAVectorUtility.SALSAError(" Invalid Kmeans Request");
                    throw (e);
                }

                initialized = true;

                Program.UseSponge = false;
                Program.ContinuousClustering = false;
                Program.MaxNumberSplitClusters = Math.Max(Program.MaxNumberSplitClusters, 1);

                ParallelClustering.RunningSolution = new ClusteringSolution(Program.UseSponge);
                ParallelClustering.SavedSolution = ParallelClustering.RunningSolution;
                ParallelClustering.BestSolution = ParallelClustering.RunningSolution;
                DAVectorUtility.SALSAPrint(0, "Clustering Solutions Created");

            }	//end Initialization

            InitializeSolution(ParallelClustering.RunningSolution);

            ClusteringSolution.ClustersDeleted = 0;
            ClusteringSolution.ClustersMoved = 0;
            ClusteringSolution.ClustersSplit = 0;

            //  Run Kmeans
            Kmeans.RunKmeans(ParallelClustering.RunningSolution.Y_k_i_, ParallelClustering.RunningSolution.OccupationCounts_k_,
                ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_, out ParallelClustering.RunningSolution.Ncent_Global,
                out ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth);
            return;

        }   // End ControlKmeans   

        public static void SetSolution(ClusteringSolution StartSolution, int Ncent_GlobalINPUT)
        {
            
            for (int RealClusterIndex = 0; RealClusterIndex < Ncent_GlobalINPUT; RealClusterIndex++)
            {
                StartSolution.P_k_[0] = 1.0;
                StartSolution.SplitPriority_k_[RealClusterIndex] = -1;
                StartSolution.Splittable_k_[RealClusterIndex] = 0;
                StartSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
                StartSolution.LocalStatus[RealClusterIndex] = 1;
                int CreatedIndex = ClusteringSolution.SetCreatedIndex(RealClusterIndex);
            }

            StartSolution.Ncent_Global = Ncent_GlobalINPUT;
            StartSolution.Ncent_ThisNode = Ncent_GlobalINPUT;
            StartSolution.SetActiveClusters();
        }

        public static void InitializeSolution(ClusteringSolution StartSolution)
        {
            // Temperature always 1 for K means
            StartSolution.SpongeCluster = -1;

            StartSolution.Temperature = 1.0;
            Program.ActualStartTemperature = StartSolution.Temperature; // For Output
            Program.TargetEndTemperature = StartSolution.Temperature;   // For Output

            DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global.ToString() + " Kmeans");
            StartSolution.ActualCoolingFactor = Program.InitialCoolingFactor;

            StartSolution.PairwiseHammy = 0.0;
            ControlKmeans.SetSolution(StartSolution, Program.InitialNcent);

            return;

        }   // End InitializeSolution

        public static void CaptureKmeans(int[] FullAssignment, out double[][] KmeanCenters)
        {   // Save Kmeans assignment and centers

            Parallel.For(0, Kmeans._parallelOptions.MaxDegreeOfParallelism, Kmeans._parallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    //                    FullAssignment[alpha + DAVectorUtility.PointStart_Process] = Kmeans.InitialPointAssignment[alpha];
                    FullAssignment[alpha + DAVectorUtility.PointStart_Process] = KmeansTriangleInequality.NearestCentertoPoint[alpha];

                }   // End loop over Points
            }); // End Sum over Threads
            if (DAVectorUtility.MPI_Size > 1)
            {
                int MPItag = 100;
                for (int mpiloop = 1; mpiloop < DAVectorUtility.MPI_Size; mpiloop++)
                {
                    if (DAVectorUtility.MPI_Rank == 0)
                    {
                        MPIPacket<int> fromsource;
                        fromsource = DAVectorUtility.MPI_communicator.Receive<MPIPacket<int>>(mpiloop, MPItag);
                        for (int index = 0; index < fromsource.NumberofPoints; index++)
                            FullAssignment[index + fromsource.FirstPoint] = fromsource.Marray[index];
                    }
                    else
                    {
                        if (DAVectorUtility.MPI_Rank == mpiloop)
                        {
                            MPIPacket<int> tosend = new MPIPacket<int>(DAVectorUtility.PointCount_Process);
                            tosend.FirstPoint = DAVectorUtility.PointStart_Process;
                            tosend.NumberofPoints = DAVectorUtility.PointCount_Process;
                            for (int index = 0; index < DAVectorUtility.PointCount_Process; index++)
                                tosend.Marray[index] = FullAssignment[index + DAVectorUtility.PointStart_Process];
                            DAVectorUtility.MPI_communicator.Send<MPIPacket<int>>(tosend, 0, MPItag);
                        }
                    }
                    DAVectorUtility.MPI_communicator.Barrier();
                }
            }

            KmeanCenters = Kmeans.ClusterCenter;

        }   // End CaptureKmeans


    }   // End ControlKmeans

    // Calculate Vector Deterministic Annealing Algorithm
    class VectorAnnealIterate
    {

        public static int[] ListofClusterstoSplit;  // List of Clusters to Split
        public static double[] EigsofClusterstoSplit;  // List of Eigenvalues of Clusters to Split
        public static int[] PrioritiesofClusterstoSplit;    // Priorities
        public static int Numberthatcanbesplit = 0; // Length of entries in ListofClusterstoSplit which is at most Program.MaxNumberSplitClusters

        public static double Tmin;  //minimal temperature 
        public static int ActualMaxNcent;   // Maximum reduced if small clusters
        public static int CountValidityFailures = 0;    // Count failures in validity checks of final solution
        public static bool OnLastleg = false;   // If True we are on final convergence step with full number of clusters
        public static bool ArithmeticError = false; // Set True if overflow
        public static double MeanClusterCount = 0.0;    // Set to Mean Cluster Count 
        public static double LocalUselessCalcs = 0.0;   // Useless Calcs this iteration
        public static double LocalUsefulCalcs = 0.0;   // Useful Calcs this iteration
        public static double LocalIgnoredCalcs = 0.0;   // Ignored Calcs this iteration
        public static double PointswithClusterCount1 = 0.0; // Number of Clusters with Cluster Count 1
        public static double C_k_Sum = 0.0; // Set to Sum of C(k)'s for debugging
        public static double MalphaDiffAvg = 0.0;   // M value change for test
        public static double AverageY_k_SquaredChange = 0.0;   // Average over Clusters of Squared Changes in Y
        public static double NumberCountChanges = 0.0;  // Average over points of # Cluster Changes
        public static double TemperatureatLastCloseTest = -1.0;

        public static int EMIterationCount = 0; // iteration number in EM method
        public static int EMIterationStepCount = 0;
        public static int EMIterationStepCount1 = -1;    // M convergence
        public static int EMIterationStepCount2 = -1;    // Y convergence
        public static int Extra_EMIterationCount = 0; // Extra iterations in EM method
        public static int ActualWaititerations = 0;
        public static int countAfterFixingClusterCount = 0;   // Count loops after end on number of centers
        public static int IterationNumberPrinted = -1;             // Iteration Number Printed
        public static int convergence = 0;  // If zero NOT converged; = 1 Converged but more to go; =2 Converged at end
        public static int CountBetweenSplits = 0; // Count actual iterations between splits. Limit WaitIterations

        public static bool HitConvergenceLoopLimit = false;
        public static bool FinalLoop = false;  // If true just converging -- No splits
        public static int SplitFailures = 0;    // Counts splitting failures
        public static bool DistributeNextTime = false;
        public static bool initialized = false;
        public static bool HammyNotSet = true;  // If True Hamiltonian Set
        public static double ChangeinHammy = 0.0;   // Change in Hamiltonian

        public static bool StopReason1 = false;
        public static bool StopReason2 = false;
        public static bool StopReason3 = false;
        public static bool StopReason4 = false;
        public static bool StopReason5 = false;
        public static bool StopReason6 = false;

        /*
         * The basic function in the program. It implements all steps of the Deterministic Annealing Vector Clustering algorithm.
         *  
         */
        public void ControlVectorSpongeDA()
        {

            //allocate memory on first and indeed only call
            if (!initialized)
            {
                if(!Program.ContinuousClustering)
                {
                    Exception e = DAVectorUtility.SALSAError(" Invalid Continuous Clustering");
                    throw (e);
                }
                
                initialized = true;
                int cachelinesize = Program.cachelinesize;

                ParallelClustering.RunningSolution = new ClusteringSolution(Program.UseSponge);
                ParallelClustering.SavedSolution = new ClusteringSolution(Program.UseSponge);
                ParallelClustering.BestSolution = new ClusteringSolution(Program.UseSponge);
                DAVectorUtility.SALSAPrint(0, "Clustering Solutions Created");

                Program.InitialNcent = 1;                           //  We only support starting with one center (plus sponge if needed)
                if (Program.UseSponge)
                    ++Program.InitialNcent;

                Program.MaxNumberSplitClusters = Math.Max(Program.MaxNumberSplitClusters, 1);
                VectorAnnealIterate.ListofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
                VectorAnnealIterate.PrioritiesofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
                VectorAnnealIterate.EigsofClusterstoSplit = new double[Program.MaxNumberSplitClusters];

            }	//end Initialization

            //  Do EM calculation
            VectorAnnealIterate.ActualMaxNcent = Program.maxNcentperNode;
            VectorAnnealIterate.OnLastleg = false;

            //  Set up Triangle Inequality
            if (Program.UseTriangleInequality_DA > 0)
            {
                DAVectorUtility.SALSAPrint(0, "Triangle Inequality Initialized Option " + Program.UseTriangleInequality_DA.ToString());
                DATriangleInequality.SetTriangleInequalityParameters(Program.UseTriangleInequality_DA, Program.MaxClusterLBsperPoint_DA, Program.MaxCentersperCenter_DA,
                     Program.TriangleInequality_Delta1_old_DA, Program.TriangleInequality_Delta1_current_DA, Program.OldCenterOption_DA);

                DATriangleInequality.InitializeTriangleInequality(Program.PointPosition, Program._parallelOptions, ParallelClustering.RunningSolution.Y_k_i_, ParallelClustering.RunningSolution.Sigma_k_i_,
                    ParallelClustering.RunningSolution.LocalStatus, ClusteringSolution.RealClusterIndices, ParallelClustering.RunningSolution,
                    Program.maxNcentTOTAL, Program.maxNcentTOTALforParallelism_DA, Program.ParameterVectorDimension);

                DATriangleInequality.SetExternalFunctions(ParallelClustering.RunningSolution.SetClusterRadius);
            }

            // Initialize Clusters
            if (Program.RestartTemperature > 0.0)
                // Restart from Previous Runs
            {
                VectorAnnealIterate.Restart();
                ClusteringSolution.CopySolution(ParallelClustering.RunningSolution, ParallelClustering.BestSolution);

                //  Set up initial cluster and Sponge if it exists
                if (Program.UseTriangleInequality_DA > 0)
                {
                    for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                        DATriangleInequality.AddCenter(ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]);
                }
                DistributedClusteringSolution.ManageMajorSynchronization(true);
                VectorAnnealIterate.ActualWaititerations = 10 * Program.Waititerations;
                ParallelClustering.RunningSolution.ActualCoolingFactor = Program.FineCoolingFactor1;
                if(ParallelClustering.RunningSolution.Temperature < Program.CoolingTemperatureSwitch)
                    ParallelClustering.RunningSolution.ActualCoolingFactor = Program.FineCoolingFactor2;
            }
            else
            {
                // Set Initial Temperature Cluster Center P_k_ and occupations Malpha_k_
                VectorAnnealIterate.InitializeSolution(ParallelClustering.RunningSolution);

                //  Set up initial cluster and Sponge if it exists
                if (Program.UseTriangleInequality_DA > 0)
                {
                    for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                        DATriangleInequality.AddCenter(ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex]);
                }

                DistributedClusteringSolution.ManageMajorSynchronization(true);
                VectorAnnealIterate.ActualWaititerations = Program.Waititerations; // Wait this number of Temperature iterations before splitting
                // This is changed in special circumstances where more convergence needed
            }


            VectorAnnealIterate.EMIterationCount = 0; // Initialize iteration count -- there is no limit
            VectorAnnealIterate.EMIterationStepCount = 0;   // This is number of counts in current step converging EM at given temperature
            VectorAnnealIterate.EMIterationStepCount1 = -1;
            VectorAnnealIterate.EMIterationStepCount2 = -1;
            VectorAnnealIterate.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached (this decreases temperature)
            int HammyViolations = 0;    // Counts number of (illegal) consecutive INCREASES in Hamiltonian
            VectorAnnealIterate.CountBetweenSplits = 0; // Count iterations between splits. Limit WaitIterations
            bool CurrentJobFinished = false;

            //  Variables tracking convergence
            VectorAnnealIterate.HitConvergenceLoopLimit = false;// If true, EMIterationStepCount has hit limit
            VectorAnnealIterate.SplitFailures = 0;
            VectorAnnealIterate.FinalLoop = false; // If true we are in a special loop freezing current case

            VectorAnnealIterate.HammyNotSet = true;    // To note that OldHammy not set

            ClusteringSolution.ClustersDeleted = 0;
            ClusteringSolution.ClustersMoved = 0;
            ClusteringSolution.ClustersSplit = 0;

            //  Proper Deterministic Annealing
            //	Loop over EM calculations
            while (Program.Refinement)
            {
                VectorAnnealIterate.EMIterationCount++; // Increment EM Loop Count
                VectorAnnealIterate.EMIterationStepCount++;    // Iterate within a fixed temperature converging M and p

                //  Set Cooling Factors
                Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
                Program.FineCoolingFactor = Program.FineCoolingFactor1;
                if (ParallelClustering.RunningSolution.Temperature < Program.CoolingTemperatureSwitch)
                {
                    Program.InitialCoolingFactor = Program.InitialCoolingFactor2;
                    Program.FineCoolingFactor = Program.FineCoolingFactor2;
                }

                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                ParallelClustering.RunningSolution.DiffMalpha_k_Set = DAVectorUtility.MPI_communicator.Allreduce<int>(ParallelClustering.RunningSolution.DiffMalpha_k_Set, Operation<int>.Min);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

                DAVectorEMIterate();   // Basic Update loop
                ParallelClustering.RunningSolution.SolutionSet = true;
                ParallelClustering.RunningSolution.IterationSetAt = VectorAnnealIterate.EMIterationCount;

                //	Now see if we are done
                string ReasonforConvergence;
                double MvalueChange = Program.Malpha_MaxChange;
                if (VectorAnnealIterate.FinalLoop)
                    MvalueChange = Program.Malpha_MaxChange1;
                convergence = VectorAnnealIterate.convergenceTest(MvalueChange, out ReasonforConvergence);

                DAVectorUtility.SynchronizeMPIvariable(ref convergence);
                //  If not converged at this temperature and cluster number just proceed with while(true) iteration
                if (convergence == 0)
                    continue;

                if (VectorAnnealIterate.EMIterationStepCount1 == -2)
                    ++Program.NumberMfailures;
                if (VectorAnnealIterate.EMIterationStepCount1 >= 0)
                {
                    Program.AccumulateMvalues += VectorAnnealIterate.EMIterationStepCount1;
                    ++Program.NumberMsuccesses;
                }
                VectorAnnealIterate.EMIterationStepCount1 = -1;
                if (VectorAnnealIterate.EMIterationStepCount2 == -2)
                    ++Program.NumberYfailures;
                if (VectorAnnealIterate.EMIterationStepCount2 >= 0)
                {
                    Program.AccumulateYvalues += VectorAnnealIterate.EMIterationStepCount2;
                    ++Program.NumberYsuccesses;
                }
                VectorAnnealIterate.EMIterationStepCount2 = -1;

                // Case we are converged = 1 in middle or = 2 at end of refinement
                Program.NumberIterationSteps++; // Accumulate diagnostic statistics
                Program.IterationsperStepSum += VectorAnnealIterate.EMIterationStepCount;

                //  Update Solution
                DistributedClusteringSolution.ManageMajorSynchronization(true);

                //  Calculate ClusterSquaredWidth_k_
                //  Used in Print and Calculation of Objective Function and Split decisions
                ParallelClustering.RunningSolution.SetClusterWidths();

                VectorAnnealIterate.EMIterationStepCount = 0;
                VectorAnnealIterate.EMIterationStepCount1 = -1;
                VectorAnnealIterate.EMIterationStepCount2 = -1;

                //  Update Hamiltonian and see if decreasing
                bool decreasing = UpdateHamiltonian();

                //  Case when at end of stage (e.g. given number of clusters or of whole job). This still needs to be iterated to low Temperatures
                if (convergence == 2 || VectorAnnealIterate.FinalLoop)
                {
                    if (VectorAnnealIterate.countAfterFixingClusterCount < Program.Iterationatend)         //    do Iterationatend iterations after reaching the maximum cluster#
                    {
                        if (VectorAnnealIterate.countAfterFixingClusterCount == 0)
                        {   // First step of "justconverging stage"
                            HammyViolations = 0;
                            Program.ActualEndTemperature = ParallelClustering.RunningSolution.Temperature;
                            DAVectorUtility.InterimTiming();
                            Program.TimeatSplittingStop = DAVectorUtility.HPDuration;
                            Program.InitialCoolingFactor = Program.InitialCoolingFactor1;
                            ParallelClustering.RunningSolution.ActualCoolingFactor = Program.InitialCoolingFactor;
                        } 

                        //  Check Freezing measure -- only use to stop if Y_k_ and P(k) converged
                        int toobigfreezing01 = 0;
                        int toobigfreezing2 = 0;
                        double freezemax01 = 0.0;
                        double freezemax2 = 0.0;
                        for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                        {
                            if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0 || ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 2)
                                continue;
                            if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                            {
                                freezemax2 = Math.Max(freezemax2, ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex]);
                                if (ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] > Program.FreezingLimit)
                                    ++toobigfreezing2;
                            }
                            else
                            {
                                freezemax01 = Math.Max(freezemax01, ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex]);
                                if (ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] > Program.FreezingLimit)
                                    ++toobigfreezing01;
                            }
                        }

                        if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                        {
                            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                            freezemax2 = DAVectorUtility.MPI_communicator.Allreduce<double>(toobigfreezing2, Operation<double>.Max);
                            toobigfreezing2 = DAVectorUtility.MPI_communicator.Allreduce<int>(toobigfreezing2, Operation<int>.Add);
                            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
                        }
                        double freezemax = Math.Max(freezemax01, freezemax2);
                        DAVectorUtility.SynchronizeMPIvariable(ref toobigfreezing01);
                        int toobigfreezing = toobigfreezing01 + toobigfreezing2;
                        if (toobigfreezing == 0)
                            DAVectorUtility.SALSAPrint(1, " ****** Stop as Freezing Measures small " + ParallelClustering.RunningSolution.FreezingMeasure_k_[0]);

                        bool TemperatureLimit = (ParallelClustering.RunningSolution.Temperature < VectorAnnealIterate.Tmin);
                        DAVectorUtility.SynchronizeMPIvariable(ref TemperatureLimit);
                        VectorAnnealIterate.countAfterFixingClusterCount++;
                        if(decreasing&&(!TemperatureLimit))
                            VectorAnnealIterate.countAfterFixingClusterCount = Math.Min(Program.Iterationatend - 10, VectorAnnealIterate.countAfterFixingClusterCount);
                        string looptype = "Final Clean Up";
                        int printtype = -1;
                        if ((VectorAnnealIterate.countAfterFixingClusterCount == (Program.Iterationatend - 1))
                            || (toobigfreezing == 0) || (!decreasing))
                            printtype = ParallelClustering.RunningSolution.IterationSetAt;
                        looptype += " Frz Viol# " + toobigfreezing.ToString() + " Max " + freezemax.ToString("E3") + " " + ReasonforConvergence;
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
                        DAVectorUtility.SynchronizeMPIvariable(ref decreasing);
                        if (DAVectorUtility.MPI_Rank == 0)
                        {
                            CurrentJobFinished = !((toobigfreezing > 0) && decreasing);
                            if (VectorAnnealIterate.ArithmeticError)
                                CurrentJobFinished = true;
                        }
                        DAVectorUtility.SynchronizeMPIvariable(ref CurrentJobFinished);

                        // Iteration at end still going Reduce Temperature and proceed
                        if (!CurrentJobFinished)
                        {
                            ReduceTemperature();
                            continue;   // Continue while over EM Loop
                        }

                        // We have finished this job
                        if (VectorAnnealIterate.ArithmeticError)
                            StopReason1 = true;
                        if (!decreasing)
                            StopReason2 = true;
                        if (toobigfreezing > 0)
                            StopReason3 = true;

                    }   // End Convergence=2  or justconverging doing Final Program.Iterationatend iterations counted by countAfterFixingClusterCount

                    else
                    {   // Case when Program.Iterationatend iterations counted by countAfterFixingClusterCount exceeded
                        StopReason6 = true;
                        CurrentJobFinished = true;
                    }   // End Convergence=2 or justconverging case where extra iteration count completed
                    

                    // This completes processing for convergence=2 and justconverging=true case
                    if (CurrentJobFinished)
                    {
                        if (ParallelClustering.RunningSolution.Temperature < VectorAnnealIterate.Tmin)
                            StopReason4 = true;
                        if (VectorAnnealIterate.SplitFailures > 1)
                            StopReason5 = true;


                        // Real end of job except check validity
                        DAVectorUtility.StartSubTimer(1);
                        bool validity = VectorAnnealIterate.CheckValidSolution(true);
                        DAVectorUtility.StopSubTimer(1);
                        if ((!validity) || (!VectorAnnealIterate.FinalLoop))
                        {   // Restart with small clusters removed or do a final improved precision run
                            if (!validity)
                                DAVectorUtility.SALSAPrint(1, "Restart after small clusters removed");
                            if (!VectorAnnealIterate.FinalLoop)
                                DAVectorUtility.SALSAPrint(1, "Add in Y convergence Test");
                            VectorAnnealIterate.ActualMaxNcent = ParallelClustering.RunningSolution.Ncent_ThisNode;
                            ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;
                            ParallelClustering.RunningSolution.YPreviousSet = -1;
                            VectorAnnealIterate.OnLastleg = true;
                            VectorAnnealIterate.countAfterFixingClusterCount = Math.Max(1, Program.Iterationatend-100);
                            VectorAnnealIterate.CountBetweenSplits = 0;
                            VectorAnnealIterate.HammyNotSet = true;
                            HammyViolations = 0;
                            VectorAnnealIterate.FinalLoop = true;
                            CurrentJobFinished = false;
                            continue;
                        }

                        if (StopReason1)
                            Program.FinalReason += "Arithmetic Error ";
                        if (StopReason2)
                            Program.FinalReason += "Stop as Hamiltonian Increasing with Change " + ChangeinHammy.ToString("E4") + " ";
                        if (StopReason3)
                            Program.FinalReason += "Stop as Freezing Measures smaller than " + Program.FreezingLimit.ToString() + " ";
                        if (StopReason4)
                            Program.FinalReason += "Tmin Reached " + VectorAnnealIterate.Tmin.ToString("E4") + " ";
                        if (StopReason5)
                            Program.FinalReason += "Consecutive Split Failures ";
                        if (StopReason6)
                            Program.FinalReason += " Stop as Final Iteration Count larger than " + Program.Iterationatend.ToString() + " ";
                        DAVectorUtility.SALSAPrint(1, Program.FinalReason);
                        break;
                    }   // end case when task finished

                }   // End justconverging = true or convergence =2 case

                DAVectorUtility.StartSubTimer(1);
                bool convergedvalidity = VectorAnnealIterate.CheckValidSolution(false);
                DAVectorUtility.StopSubTimer(1);
                if (!convergedvalidity)
                {   // Restart with small clusters removed
                    PrintIteration(" Restart After Cluster Removal ", -1);
                    CountBetweenSplits = 0;
                    VectorAnnealIterate.ActualWaititerations = Program.Waititerations_Converge;
                    HammyNotSet = true;
                    continue;
                }
                //  This section results either in EM Loop Continue (for ongoing refinement for a nonsplittable case) or EM Loop Break (to end) 

                //  Convergence = 1 Case
                //	Converged so test for split -  -- take results from Rank 0 but rest do arithmetic
                //  For restarted task that decision has already been made

                // Need to decide if split to occur
                bool ResultofSplittingTest = false;
                ++CountBetweenSplits;
                if (CountBetweenSplits > VectorAnnealIterate.ActualWaititerations)
                {
                    CountBetweenSplits = 0;
                    VectorAnnealIterate.ActualWaititerations = Program.Waititerations;
                }
                if ((CountBetweenSplits > 0) && (VectorAnnealIterate.ActualWaititerations > 0) && (ParallelClustering.RunningSolution.Ncent_Global > 1))
                {
                    PrintIteration("Ongoing Annealing", -1);
                    ReduceTemperature();
                    continue;   // EM Loop for compulsory iterations between splits
                }

                bool ChangedClusters = VectorAnnealIterate.CheckNumberofClustersTooBig();
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                ChangedClusters = DAVectorUtility.MPI_communicator.Allreduce<bool>(ChangedClusters, Operation<bool>.LogicalOr);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
                if (ChangedClusters)
                    DistributedClusteringSolution.ManageMajorSynchronization(true);

                //  Decide if to split
                DAVectorUtility.StartSubTimer(0);
                ResultofSplittingTest = shouldweSplit();
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                {
                    DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                    DAVectorUtility.MPI_communicator.Allreduce<bool>(ResultofSplittingTest, Operation<bool>.LogicalOr);
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
                }
                else
                {
                    DAVectorUtility.SynchronizeMPIvariable(ref ResultofSplittingTest);
                    DAVectorUtility.SynchronizeMPIvariable(ref ParallelClustering.RunningSolution.ClustertoSplit);
                }
                DAVectorUtility.StopSubTimer(0);

                //  Diagnostic Output for splitting

                if( ResultofSplittingTest && (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0) )
                    this.diagnosticsplitprint(ResultofSplittingTest);


                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                ResultofSplittingTest = DAVectorUtility.MPI_communicator.Allreduce<bool>(ResultofSplittingTest, Operation<bool>.LogicalOr);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

                //	If split indicated perform this                          
                if (ResultofSplittingTest == true)
                {
                    if (ParallelClustering.RunningSolution.Ncent_Global > 1)
                    {
                        if (Program.ClusterCountOutput > 0)
                            VectorAnnealIterate.OutputClusteringResults("Inter2");
                    }

                    //  Need to perform already determined split (maybe zero in one node if distributed)
                    DAVectorUtility.StartSubTimer(5);

                    ClusteringSolution.ClustersSplit = 0;
                    string SplitString = "";
                    if (VectorAnnealIterate.Numberthatcanbesplit > 0)
                    {
                        for (int splitlistloop = 0; splitlistloop < VectorAnnealIterate.Numberthatcanbesplit; splitlistloop++)
                        {
                            ParallelClustering.RunningSolution.ClustertoSplit = VectorAnnealIterate.ListofClusterstoSplit[splitlistloop];
                            SplitString += ParallelClustering.RunningSolution.ClustertoSplit.ToString() + "(" + ParallelClustering.RunningSolution.C_k_[ParallelClustering.RunningSolution.ClustertoSplit].ToString()
                                + ")" + " " + VectorAnnealIterate.EigsofClusterstoSplit[splitlistloop].ToString("E4") + " * ";
                            this.dothesplit();

                        }
                        VectorAnnealIterate.Numberthatcanbesplit = 0;
                    }
                    if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
                    {
                        int TotalNumberSplit = ClusteringSolution.ClustersSplit;
                        if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                        {
                            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                            TotalNumberSplit = DAVectorUtility.MPI_communicator.Allreduce<int>(TotalNumberSplit, Operation<int>.Add);
                            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
                            if (TotalNumberSplit > 0)
                                DAVectorUtility.SALSASyncPrint(1, "Clusters Split ", SplitString);
                        }
                        else
                        {
                            if (TotalNumberSplit > 0) 
                                DAVectorUtility.SALSAPrint(1, "Clusters Split " + SplitString);
                        }
                    }
                    DAVectorUtility.StopSubTimer(5);

                    DistributedClusteringSolution.ManageMajorSynchronization(true);
                    VectorAnnealIterate.SplitFailures = 0;
                    VectorAnnealIterate.ActualWaititerations = Program.Waititerations;
                }   // End ResultofSplittingTest == true (either changed number of clusters or spun off a convergence task)

                //  Final portion of loop changing Temperature if needed
                if (ResultofSplittingTest == false)
                {   //	Reduce T and continue iterating if converged and no split
                    ReduceTemperature();
                }
                else
                {   // Don't reduce T as clusters split
                    ParallelClustering.RunningSolution.ActualCoolingFactor = Program.FineCoolingFactor;
                    CountBetweenSplits = 0;
                }

            }   // End while EM Loop

            //  Check if solution best
            if (!ParallelClustering.BestSolution.SolutionSet)
                return;
            if (ParallelClustering.BestSolution.Ncent_Global != ParallelClustering.RunningSolution.Ncent_Global)
            {
                DAVectorUtility.SALSAPrint(1, " Best Solution Not Used as Number of Centers " + ParallelClustering.BestSolution.Ncent_Global + " Different from Running Solution with " + ParallelClustering.RunningSolution.Ncent_Global);
                return;
            }
            bool changesolution = ParallelClustering.BestSolution.PairwiseHammy < ParallelClustering.RunningSolution.PairwiseHammy;
            DAVectorUtility.SynchronizeMPIvariable(ref changesolution);
            if (changesolution)
            {
                DAVectorUtility.SALSAPrint(1, " Solution at Iteration " + ParallelClustering.BestSolution.IterationSetAt.ToString() + " Chisq " + ParallelClustering.BestSolution.PairwiseHammy.ToString("E4")
                    + " Taken rather than Iteration " + ParallelClustering.RunningSolution.IterationSetAt.ToString() + " Chisq " + ParallelClustering.RunningSolution.PairwiseHammy.ToString("E4"));
                ClusteringSolution.CopySolution(ParallelClustering.BestSolution, ParallelClustering.RunningSolution);
            }
            else if (VectorAnnealIterate.ArithmeticError)
            {
                DAVectorUtility.SALSAPrint(1, " Solution at Iteration " + ParallelClustering.BestSolution.IterationSetAt.ToString() 
                    + " Taken rather than Iteration " + ParallelClustering.RunningSolution.IterationSetAt.ToString() + " Due to Arithmetic Error ");
                ClusteringSolution.CopySolution(ParallelClustering.BestSolution, ParallelClustering.RunningSolution);
            }
            int save = Program.ClusterPrintNumber;
            Program.ClusterPrintNumber = ParallelClustering.RunningSolution.Ncent_ThisNode;
            PrintIteration(" Final Solution ", ParallelClustering.RunningSolution.IterationSetAt);
            Program.ClusterPrintNumber = save;
            return;

        }	// End of ControlVectorSpongeDA()

        public void ReduceTemperature()
        {
            // Add in a Sponge if needed
            if (VectorAnnealIterate.AddSpongeCluster())
            {
                ParallelClustering.BestSolution.SolutionSet = false;
                CountBetweenSplits = 0;
                HammyNotSet = true;
                return;
            }
            else 
            {
                if (Program.ChangeSpongeFactor(ParallelClustering.RunningSolution.ActualCoolingFactor, ParallelClustering.RunningSolution))
                    HammyNotSet = true;
            }
            //  Reduce Cluster Sigmas for those being Annealed
            if (Program.ChangeClusterSigmas(ParallelClustering.RunningSolution.ActualCoolingFactor, ParallelClustering.RunningSolution))
            {
                ParallelClustering.RunningSolution.SetClusterSizes();
                ParallelClustering.RunningSolution.SetClusterWidths();
                HammyNotSet = true;
            }

            //  Switch into Distributed Mode if necessary
            if(!ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                bool godistributed = (ParallelClustering.RunningSolution.Temperature <= Program.TemperatureLimitforDistribution) 
                    || ( (Program.ClusterLimitforDistribution > 0 ) && (ParallelClustering.RunningSolution.Ncent_Global >= Program.ClusterLimitforDistribution));
                DAVectorUtility.SynchronizeMPIvariable(ref godistributed);
                if(godistributed)
                {
                    if (VectorAnnealIterate.DistributeNextTime)
                    {
                        ParallelClustering.RunningSolution.DistributedExecutionMode = true;

                        int newsplitnumber = Program.MaxNumberSplitClusters / DAVectorUtility.MPI_Size;
                        if ((newsplitnumber * Program.MaxNumberSplitClusters) < Program.MaxNumberSplitClusters)
                            ++newsplitnumber;
                        Program.MaxNumberSplitClusters = newsplitnumber;
                        ++DAVectorUtility.MPIREDUCETiming1;
                        Program.ActualTemperatureforDistribution = ParallelClustering.RunningSolution.Temperature;
                        DAVectorUtility.InterimTiming();
                        Program.TimeatDistribution = DAVectorUtility.HPDuration;
                        Program.ActualClusterNumberforDistribution = ParallelClustering.RunningSolution.Ncent_Global;
                        DAVectorUtility.SALSAPrint(1, "Start Distributed Execution Mode " + ParallelClustering.RunningSolution.Temperature.ToString("F2") + " Iteration Count "
                            + ParallelClustering.RunningSolution.IterationSetAt.ToString() + " " + EMIterationCount.ToString() + " Clusters " + ParallelClustering.RunningSolution.Ncent_Global.ToString()
                            + " Time " + DAVectorUtility.HPDuration.ToString("F0"));
                        DistributedClusteringSolution.ManageMajorSynchronization(false);
                        ParallelClustering.BestSolution.SolutionSet = false;
                        CountBetweenSplits = 0;
                        HammyNotSet = true;
                        return;
                    }
                    else
                    {
                        VectorAnnealIterate.DistributeNextTime = true;
                        VectorAnnealIterate.CompleteCleanUp();
                        return;
                    }
                }
            }

            //  Generate Magic Temperature Actions
            if (Program.magicindex < Program.MagicTemperatures.Length)
            {
                bool abracadabra = ParallelClustering.RunningSolution.Temperature <= Program.MagicTemperatures[Program.magicindex];
                DAVectorUtility.SynchronizeMPIvariable(ref abracadabra);
                if (abracadabra)
                {
//                    Dist.OutputClusterLabels("Temp" + ParallelClustering.RunningDAVectorSolt.Temperature.ToString("F2"));
                    VectorAnnealIterate.CompleteCleanUp();
                    if(ParallelClustering.RunningSolution.DistributedExecutionMode)
                        DistributedClusteringSolution.ManageMajorSynchronization(true);
                    CountBetweenSplits = 0;
                    ++Program.magicindex;
                    return;
                }
            }
            ParallelClustering.RunningSolution.Temperature = ParallelClustering.RunningSolution.ActualCoolingFactor * ParallelClustering.RunningSolution.Temperature;
            ++Program.NumberTemperatureSteps;
            DAVectorUtility.TemperatureValues.Add(ParallelClustering.RunningSolution.Temperature);
            DAVectorUtility.ClusterCountValues.Add(ParallelClustering.RunningSolution.Ncent_Global);

        }   // End ReduceTemperature

        public void PrintIteration(string looptype, int linecheck)
        {
            if (linecheck < 0)
            {
                if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval != 0)
                    return;
            }
            else if (linecheck == VectorAnnealIterate.IterationNumberPrinted)
                return;
            VectorAnnealIterate.IterationNumberPrinted = ParallelClustering.RunningSolution.IterationSetAt;

            DAVectorUtility.InterimTiming();
            string endinfo = "";
            if (!Program.CalculateIndividualWidths)
                endinfo = " Average Width " + ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth.ToString("E3");
            else
            {
                endinfo = " Average Widths ";
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    endinfo += ParallelClustering.RunningSolution.AverageWidth[VectorIndex].ToString("E3") + " ";
            }
            if (Program.SigmaMethod == 3)
                endinfo += " Sigma[0] Coeff " + Program.SigmaVectorParameters_i_[0].ToString("E4");

            DAVectorUtility.SALSAPrint(1, "B) Clusters " + ParallelClustering.RunningSolution.Ncent_Global.ToString() + " " + looptype + " Cluster # Chg " 
                + VectorAnnealIterate.NumberCountChanges.ToString("E3") + " M-Change " + VectorAnnealIterate.MalphaDiffAvg.ToString("E3")
                + " Y Chg " + (VectorAnnealIterate.AverageY_k_SquaredChange / ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth).ToString("E3")
                + " Cnvg " + VectorAnnealIterate.convergence.ToString() + " Iter " + VectorAnnealIterate.EMIterationCount.ToString() + " Major " + Program.NumberMajorSynchs1.ToString()
                + " T " + ParallelClustering.RunningSolution.Temperature.ToString("E4") + " PWHammy " + ParallelClustering.RunningSolution.PairwiseHammy.ToString("E4")
                + " Useful Calcs " + VectorAnnealIterate.LocalUsefulCalcs.ToString("E4") + " Useless Calcs " + VectorAnnealIterate.LocalUselessCalcs.ToString("E4")
                + " Ignored Calcs " + VectorAnnealIterate.LocalIgnoredCalcs.ToString("E4")
                + " Arithmetic Error " + VectorAnnealIterate.ArithmeticError.ToString() + " Mean Cluster Count per point " + VectorAnnealIterate.MeanClusterCount.ToString("F2")
                + " Pts with Just 1 Cluster " + VectorAnnealIterate.PointswithClusterCount1.ToString("F0") + " Sum of C(k) " + VectorAnnealIterate.C_k_Sum.ToString("F2") + endinfo
                +  " Time " + DAVectorUtility.HPDuration.ToString("F0"));

            clusterprint();

        }   // End Print Iteration from getDist

        public void diagnosticsplitprint(bool ResultofSplittingTest)
        {
            string endinfo = "";
            if (!Program.CalculateIndividualWidths)
                endinfo = " Average Width " + ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth.ToString("E3");
            else
            {
                endinfo = " Average Widths ";
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    endinfo += ParallelClustering.RunningSolution.AverageWidth[VectorIndex].ToString("E3") + " ";
            }
            if (Program.SigmaMethod == 3)
                endinfo += " Sigma[0] Coeff " + Program.SigmaVectorParameters_i_[0].ToString("E4");
            DAVectorUtility.InterimTiming();
            string nextline1 = "A) Clusters " + ParallelClustering.RunningSolution.Ncent_Global.ToString() + " Iter " + VectorAnnealIterate.EMIterationCount.ToString()
                + " Major " + Program.NumberMajorSynchs1.ToString()
                + " T " + ParallelClustering.RunningSolution.Temperature.ToString("E4") + " PWHammy " + ParallelClustering.RunningSolution.PairwiseHammy.ToString("E4")
                + " Useful Calcs " + VectorAnnealIterate.LocalUsefulCalcs.ToString("E4") + " Useless Calcs " + VectorAnnealIterate.LocalUselessCalcs.ToString("E4")
                + " Ignored Calcs " + VectorAnnealIterate.LocalIgnoredCalcs.ToString("E4")
                + " Arithmetic Error " + VectorAnnealIterate.ArithmeticError.ToString() + " Mean Cluster Count per point " + VectorAnnealIterate.MeanClusterCount.ToString("F2")
                + " Pts with Just 1 Cluster " + VectorAnnealIterate.PointswithClusterCount1.ToString("F0") + " Sum of C(k) " + VectorAnnealIterate.C_k_Sum.ToString("F2") + endinfo
                + " Time " + DAVectorUtility.HPDuration.ToString("F0");
            if (ParallelClustering.RunningSolution.ClustertoSplit >= 0)
            {
                if (ResultofSplittingTest)
                    nextline1 += " C# To Split " + ParallelClustering.RunningSolution.ClustertoSplit.ToString();
                else
                    nextline1 += " No Split";
                nextline1 += " " + ResultofSplittingTest.ToString() + " Status " + ParallelClustering.RunningSolution.Splittable_k_[ParallelClustering.RunningSolution.ClustertoSplit].ToString()
                    + " Eig " + ParallelClustering.RunningSolution.Eigenvalue_k[ParallelClustering.RunningSolution.ClustertoSplit].ToString("E4");
            }
            DAVectorUtility.SALSAPrint(1, nextline1);
            clusterprint();
            DAVectorUtility.SALSAPrint(1, " ");
            return;

        }   // End diagnosticsplitprint

        public void clusterprint()
        {
            string nextline = "";
            int ClusterIndex = ParallelClustering.RunningSolution.SpongeCluster;
            if (ClusterIndex < 0)
                ClusterIndex = 0;
            int count = 0;
            while (true)
            {
                string spongelabel = "";
                string Pformat = "F4";
                if (ParallelClustering.RunningSolution.LocalStatus[ClusterIndex] >= 0)
                {
                    if (ClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                    {
                        if (count > 0)
                        {
                            ClusterIndex++;
                            continue;
                        }
                        spongelabel = "Sponge " + Program.SpongeFactor.ToString("F2") + " ";
                        Pformat = "E3";
                    }
                    double tmp = ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[ClusterIndex];
                    if (count != 0)
                        nextline += "* ";
                    nextline += spongelabel + ClusterIndex.ToString() + "(" + ParallelClustering.RunningSolution.LocalCreatedIndex[ClusterIndex].ToString() + ") C " + ParallelClustering.RunningSolution.C_k_[ClusterIndex].ToString("F1") + " (Frz "
                        + ParallelClustering.RunningSolution.FreezingMeasure_k_[ClusterIndex].ToString("F6") + ") " + "[Wdth " + tmp.ToString("F4") + "] ";
                    if (Program.ContinuousClustering)
                        nextline += "P " + ParallelClustering.RunningSolution.P_k_[ClusterIndex].ToString(Pformat) + " ";
                }
                if ( (count >= Program.ClusterPrintNumber) || (ClusterIndex >= ParallelClustering.RunningSolution.Ncent_ThisNode-1) )
                    break;
                count++;
                if ((count == 1) && (ParallelClustering.RunningSolution.SpongeCluster >= 0))
                    ClusterIndex = 0;
                else
                    ClusterIndex++;
            }
            int ClusterLimit = Math.Max(ClusterIndex + 1, ParallelClustering.RunningSolution.Ncent_ThisNode - Program.ClusterPrintNumber);
            for (int ClusterEndIndex = ClusterLimit; ClusterEndIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; ClusterEndIndex++)
            {
                if (ClusterEndIndex == ParallelClustering.RunningSolution.SpongeCluster)
                    continue;
                if (ParallelClustering.RunningSolution.LocalStatus[ClusterEndIndex] < 0)
                    continue;
                double tmp = ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[ClusterEndIndex];
                if (ClusterEndIndex != 0)
                    nextline += "* ";
                nextline += ClusterEndIndex.ToString() + "(" + ParallelClustering.RunningSolution.LocalCreatedIndex[ClusterEndIndex].ToString() + ") C " + ParallelClustering.RunningSolution.C_k_[ClusterEndIndex].ToString("F1") + " (Frz "
                    + ParallelClustering.RunningSolution.FreezingMeasure_k_[ClusterEndIndex].ToString("F6") + ") " + "[Wdth " + tmp.ToString("F4") + "] ";
                if (Program.ContinuousClustering)
                    nextline += "P " + ParallelClustering.RunningSolution.P_k_[ClusterEndIndex].ToString("F4") + " ";
            }
            DAVectorUtility.SALSAPrint(1, nextline);

        }   // End clusterprint()

        //  Update Hamiltonian and return boolean decreasing to indicate if decreasing
        public bool UpdateHamiltonian()
        {
            double hammy2 = 0.0;
            double hammy01 = 0.0;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0 || ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 2)
                    continue;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                    hammy2 += 0.5 * ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                else
                    hammy01 += 0.5 * ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
            }
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                hammy2 = DAVectorUtility.MPI_communicator.Allreduce<double>(hammy2, Operation<double>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            ParallelClustering.RunningSolution.PairwiseHammy = hammy2 + hammy01;

            //  Set Best Solution and progress flag decreasing
            bool decreasing;
            if (HammyNotSet)
            {
                HammyNotSet = false;
                decreasing = true;
            }
            else
                decreasing = ParallelClustering.RunningSolution.PairwiseHammy <= ParallelClustering.RunningSolution.OldHammy;
            DAVectorUtility.SynchronizeMPIvariable(ref decreasing);
            if (decreasing)
                ClusteringSolution.CopySolution(ParallelClustering.RunningSolution, ParallelClustering.BestSolution);
            VectorAnnealIterate.ChangeinHammy = ParallelClustering.RunningSolution.PairwiseHammy - ParallelClustering.RunningSolution.OldHammy;
            ParallelClustering.RunningSolution.OldHammy = ParallelClustering.RunningSolution.PairwiseHammy;
            return decreasing;

        }   // End UpdateHamiltonian()

        //  Process CreatedIndex for a Point Cluster pointer
        // Note ActiveCluster is defined even for Remote Clusters
        public static void ClusterPointersforaPoint( int alpha, int IndirectClusterIndex, ref int RealClusterIndex, ref int ActiveClusterIndex, ref int RemoteIndex)
        {
            int CreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
            int MappedClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
            int ClusterIterationNo = ClusteringSolution.UniversalMapping[CreatedIndex].IterationSet;
            if (ClusterIterationNo < ClusteringSolution.CurrentIteration)
            {
                int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                string errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString() 
                    + " Created Index " + CreatedIndex.ToString() + " Actual Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Cluster Iteration " + ClusterIterationNo.ToString() + " Full Set Created Indices ";
                for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
                    errormessage += ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop].ToString() + " ";
                Exception e = DAVectorUtility.SALSAError(errormessage);
                throw (e);
            }
            if (MappedClusterIndex > 0)
            {
                RealClusterIndex = MappedClusterIndex - 1;
                if (RealClusterIndex >= ParallelClustering.RunningSolution.Ncent_ThisNode)
                {
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    string errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString() + " Bad Cluster Number " + RealClusterIndex.ToString()
                        + " Created Index " + CreatedIndex.ToString() + " Actual Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Cluster Iteration " + ClusterIterationNo.ToString() + " Full Set Created Indices ";
                    for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
                        errormessage += ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop].ToString() + " ";
                    Exception e = DAVectorUtility.SALSAError(errormessage);
                    throw (e);
                }
                ActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[RealClusterIndex];
                return;
            }
            else if (MappedClusterIndex == 0)
            {
                int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                string errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString() + " Zero Mapped Index " 
                    + " Created Index " + CreatedIndex.ToString() + " Actual Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Cluster Iteration "
                    + ClusterIterationNo.ToString() + " Full Set Created Indices ";
                for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
                    errormessage += ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop].ToString() + " ";
                Exception e = DAVectorUtility.SALSAError(errormessage);
                throw (e);
            }
            else
            {
                if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
                {
                    Exception e = DAVectorUtility.SALSAError(" DAVectorEMIterate() Illegal Created Index " + CreatedIndex.ToString() + " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString());
                    throw (e);
                }
                RemoteIndex = -MappedClusterIndex - 1;
                if (RemoteIndex >= DistributedClusteringSolution.StorageforTransportedClusters.SizeOfTransportedArray)
                {
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    string errormessage = " Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Cluster Index " + IndirectClusterIndex.ToString() + " Remote Index " + RemoteIndex.ToString()
                        + " Created Index " + CreatedIndex.ToString() + " Actual Iteration " + ClusteringSolution.CurrentIteration.ToString() + " Cluster Iteration " + ClusterIterationNo.ToString() + " Full Set Created Indices ";
                    for (int errorloop = 0; errorloop < IndirectSize; errorloop++)
                        errormessage += ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][errorloop].ToString() + " ";
                    Exception e = DAVectorUtility.SALSAError(errormessage);
                    throw (e);
                }
                ActiveClusterIndex = RemoteIndex + ClusteringSolution.NumberLocalActiveClusters;
                return;
            }

        }   // End ClusterPointersforaPoint( int alpha, int IndirectClusterIndex, int IndirectSize, ref int RealClusterIndex, ref int ActiveClusterIndex, ref int RemoteIndex)

        //  Update in EM order Malpha_k_ setting Previous_Malpha_k_, P_k_ Y_k_ and Differences in Malpha for convergence test
        //  Set new values of C_k_, P_k_ Y_k_ and Freezing Factors
        public void DAVectorEMIterate()
        {
            DAVectorUtility.StartSubTimer(2);
            double SpongeTerm = 0.0;
            if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                SpongeTerm = Program.SpongeFactor * Program.SpongeFactor;
            double MinP_k_ = 0.0000001;

            GlobalReductions.FindIndirectVectorDoubleSum FindDiffMalpha_k_ = null;
            GlobalReductions.FindIndirectVectorDoubleSum FindC_k_ = null;
            GlobalReductions.FindIndirectVectorDoubleSum FindFreezingMeasure_k_ = null;
            DistributedReductions.FindIndirectMultiVectorDoubleSum Find3Components = null;
            int BeginFindDiffMalpha_k_ = -1;
            int BeginFindC_k_ = -1;
            int BeginFindFreezingMeasure_k_ = -1;

            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                Find3Components = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginFindDiffMalpha_k_ = Find3Components.AddComponents(1);
                BeginFindC_k_ = Find3Components.AddComponents(1);
                BeginFindFreezingMeasure_k_ = Find3Components.AddComponents(1);
                Find3Components.NodeInitialize();
            }
            else
            {
                FindDiffMalpha_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount,
                    ClusteringSolution.NumberLocalActiveClusters);
                FindC_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount,
                    ClusteringSolution.NumberLocalActiveClusters);
                FindFreezingMeasure_k_ = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount,
                    ClusteringSolution.NumberLocalActiveClusters);
            }

            GlobalReductions.FindBoolOr FindArithmeticError = new GlobalReductions.FindBoolOr(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleMean FindMeanClusterCount = new GlobalReductions.FindDoubleMean(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleSum FindSinglePoints = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

            VectorAnnealIterate.LocalUsefulCalcs = 0.0;
            VectorAnnealIterate.LocalUselessCalcs = 0.0;
            VectorAnnealIterate.LocalIgnoredCalcs = 0.0;
            GlobalReductions.FindDoubleSum FindUsefulCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleSum FindUselessCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);
            GlobalReductions.FindDoubleSum FindIgnoredCalcs = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (Threadno) =>
            {

                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    Find3Components.ThreadInitialize((int)Threadno);
                else
                {
                    FindDiffMalpha_k_.startthread((int)Threadno);
                    FindC_k_.startthread((int)Threadno);
                    FindFreezingMeasure_k_.startthread((int)Threadno);
                }
                
                //  Arrays used in each thread and re-used by points in thread
                int ArraySize = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, ClusteringSolution.MaximumClustersperPoint);
                double[] AbsChangeinMalpha_k_ = new double[ArraySize];  //    Store change in Malpha for each point
                double[] Malphax1minusMalpha = new double[ArraySize];   // For Freezing Calculation
                double[] Save_CurrentP_k_TimesExponential = new double[ArraySize];  // Save P values for clusters
                int[] ActiveClustersperPoint = new int[ArraySize];  // For storing Cluster pointers used in this point

                double[] Save_Term_NegativeExponential = new double[ArraySize]; // Save terms in  exponential
                int[] ThreadStorePosition = new int[ArraySize];  // Used to store location of Thread storage in distributed Execution Case

                //  Loop over Points
                int indexlen = DAVectorUtility.PointsperThread[Threadno];
                int beginpoint = DAVectorUtility.StartPointperThread[Threadno] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int MinimumTerm_IndirectClusterIndex = -1;
                    double MinimumTerm_P_k_ = 0.0;
                    double MinimumTerm_NegativeExponential = 0.0;

                    //  Number of Clusters used for this point
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];

                    //  Accumulate mean number of clusters per point
                    FindMeanClusterCount.addapoint((int) Threadno, (double)IndirectSize);

                    //  Set count of points with just one cluster
                    double value = 0.0;
                    if (IndirectSize == 1)
                        value = 1.0;
                    FindSinglePoints.addapoint((int) Threadno, value);
                    double NumSponge = 0.0;

                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {   // First Pass loop sets up a safer piece of arithmetic by subtracting smallest term in exponential
                        // It will cancel in top and bottom when calculating M

                        int RealClusterIndex = -1;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        ClusterPointersforaPoint( alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                        if (RemoteIndex < 0)
                        {
                            ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                                ThreadStorePosition[IndirectClusterIndex] = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][Threadno];
                        }
                        else
                            ThreadStorePosition[IndirectClusterIndex] = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][Threadno];

                        double Term_NegativeExponential = 0.0;
                        if ((ParallelClustering.RunningSolution.SpongeCluster >= 0) && (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster))
                        {
                            Term_NegativeExponential = SpongeTerm;
                            NumSponge = 1.0;
                        }
                        else
                            Term_NegativeExponential = DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(alpha, ActiveClusterIndex, ParallelClustering.RunningSolution);
                        Term_NegativeExponential = Term_NegativeExponential / ParallelClustering.RunningSolution.Temperature;
                        Save_Term_NegativeExponential[IndirectClusterIndex] = Term_NegativeExponential;

                        //  Set P_k_ for this cluster
                        double CurrentP_k_;
                        if (RemoteIndex < 0)
                            CurrentP_k_ = ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                        else
                            CurrentP_k_ = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedP_t[RemoteIndex];
                        Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = CurrentP_k_;

                        //  Set Selected (aka Best) Index to use for subtraction
                        if((MinimumTerm_IndirectClusterIndex == -1)||(Term_NegativeExponential < MinimumTerm_NegativeExponential))
                        {
                            MinimumTerm_NegativeExponential = Term_NegativeExponential;
                            MinimumTerm_IndirectClusterIndex = IndirectClusterIndex;
                            MinimumTerm_P_k_ = CurrentP_k_;
                        }

                    }   // End First Pass Loop over Clusters for this point

                    bool LocalArithmeticError = false;
                    double MalphaDenominator = 0.0;
                    bool ThereisanOKTerm_NegativeExponential = false;

                    //  Second Pass: Now calculate quantities subtracting smallest term at top and bottom
                    double NumUseful = 0.0;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        double ExponentiatedTerm;
                        double SubtractedTerm_NegativeExponential = Save_Term_NegativeExponential[IndirectClusterIndex] - MinimumTerm_NegativeExponential;   // Subtract minimum term -- cancels between numerator and denominator


                        if (SubtractedTerm_NegativeExponential < (2.0 * Program.ExpArgumentCut1))
                        {
                            NumUseful += 1.0;
                            ExponentiatedTerm = Math.Exp(-0.5 * SubtractedTerm_NegativeExponential);
                            if (Double.IsNaN(ExponentiatedTerm) || Double.IsInfinity(ExponentiatedTerm))
                            {
                                ExponentiatedTerm = 0.0;
                                LocalArithmeticError = true;
                            }
                            else
                                ThereisanOKTerm_NegativeExponential = true;
                        }
                        else
                            ExponentiatedTerm = 0.0;

                        Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = Math.Max(MinP_k_, Save_CurrentP_k_TimesExponential[IndirectClusterIndex]) * ExponentiatedTerm;
                        MalphaDenominator += Save_CurrentP_k_TimesExponential[IndirectClusterIndex];

                    }   // End Second Pass Loop over Clusters for this point

                    //  Sponge Term is never rejected
                    FindUsefulCalcs.addapoint( (int) Threadno, NumUseful - NumSponge);
                    FindUselessCalcs.addapoint((int)Threadno, (double)IndirectSize - NumUseful);
                    FindIgnoredCalcs.addapoint((int)Threadno, (double)ParallelClustering.RunningSolution.Ncent_ThisNode - (double)IndirectSize); 

                    //  Cope with ill defined arithmetic
                    if (ThereisanOKTerm_NegativeExponential)
                    {
                        MalphaDenominator = 1.0 / MalphaDenominator;
                        if (Double.IsNaN(MalphaDenominator) || Double.IsInfinity(MalphaDenominator))
                        {
                            ThereisanOKTerm_NegativeExponential = false;
                            LocalArithmeticError = true;
                        }
                    }
                    if(!ThereisanOKTerm_NegativeExponential)
                    {
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                            Save_CurrentP_k_TimesExponential[IndirectClusterIndex] = 0.0;
                        if (MinimumTerm_IndirectClusterIndex < 0)
                        {
                            Exception e = DAVectorUtility.SALSAError(" Error in Choosing Default Index Point " +  (alpha + DAVectorUtility.PointStart_Process).ToString()
                                + " Number of Clusters " + IndirectSize.ToString());
                            throw e;
                        }
                        Save_CurrentP_k_TimesExponential[MinimumTerm_IndirectClusterIndex] = Math.Max( MinP_k_, MinimumTerm_P_k_);
                        MalphaDenominator = 1.0 / Save_CurrentP_k_TimesExponential[MinimumTerm_IndirectClusterIndex];
                    }
                    //  End Cope with ill defined arithmetic

                    //  Third Pass -- Finally set M and its change
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        double CandidateMvalue = Save_CurrentP_k_TimesExponential[IndirectClusterIndex] * MalphaDenominator;
                        if (Double.IsNaN(CandidateMvalue) || Double.IsInfinity(CandidateMvalue) )
                        {
                            string message = "";
                            for(int loop = 0; loop < IndirectSize; loop++)
                            {
                                int CreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][loop];
                                double Cvalue = 0.0;
                                double Pvalue = 0.0;
                                int ClusterIndex = ClusteringSolution.UniversalMapping[CreatedIndex].Availability;
                                if (ClusterIndex > 0)
                                {
                                    Cvalue = ParallelClustering.RunningSolution.C_k_[ClusterIndex - 1];
                                    Pvalue = ParallelClustering.RunningSolution.P_k_[ClusterIndex - 1];
                                }
                                else
                                    Pvalue = DistributedClusteringSolution.StorageforTransportedClusters.TotalTransportedP_t[-ClusterIndex - 1];
                                message += CreatedIndex.ToString() + " P " + Pvalue.ToString("F5") + " C " + Cvalue.ToString("F5")
                                    + " Temp " + Save_CurrentP_k_TimesExponential[loop].ToString("E4") + " Term " + Save_Term_NegativeExponential[loop].ToString("E6") + " Old Malpha "
                                    + ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][loop].ToString("F5") + " * ";
                            }
                            Exception e = DAVectorUtility.SALSAError("Arithmetic Error Point " + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Number of Clusters " + IndirectSize.ToString()
                                + " Indirect Index " + IndirectClusterIndex.ToString() + " Created Index " + ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex].ToString()
                                + " Selected Indirect Index " + MinimumTerm_IndirectClusterIndex.ToString() + " ZeeX " + MalphaDenominator.ToString("E4") + "\n" + message);
                            throw (e);
                        }

                        //  Finally set basic accumulation values
                        if (ParallelClustering.RunningSolution.DiffMalpha_k_Set >= 0)
                            AbsChangeinMalpha_k_[IndirectClusterIndex] = Math.Abs(CandidateMvalue - ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex]);

                        ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = CandidateMvalue;
                        Malphax1minusMalpha[IndirectClusterIndex] = CandidateMvalue * (1.0 - CandidateMvalue);
                    }
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                        ParallelClustering.RunningSolution.LegalCluster(alpha, IndirectClusterIndex);

                    //  Finally Load Contributions -- Distributed Cluster Execution
                    if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    {
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                        {
                            int LocationInThreadArray = ThreadStorePosition[IndirectClusterIndex];
                            if ((LocationInThreadArray < 0) || (LocationInThreadArray >= DistributedClusteringSolution.NodeAccMetaData.NumberofPointsperThread[Threadno]))
                            {
                                Exception e = DAVectorUtility.SALSAError("Bad Thread Location " + LocationInThreadArray.ToString() + " Thread " + Threadno.ToString()
                                    + " Number Clusters " + IndirectSize.ToString() + " Cluster Index " + ActiveClustersperPoint[IndirectClusterIndex].ToString()
                                    + " Number Local " + ClusteringSolution.NumberLocalActiveClusters.ToString());
                                throw(e);
                            }
                            if (ParallelClustering.RunningSolution.DiffMalpha_k_Set >= 0) 
                                Find3Components.addapoint((int)Threadno, LocationInThreadArray, BeginFindDiffMalpha_k_, AbsChangeinMalpha_k_[IndirectClusterIndex]);
                            Find3Components.addapoint((int)Threadno, LocationInThreadArray, BeginFindC_k_,
                                ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex]);
                            Find3Components.addapoint((int)Threadno, LocationInThreadArray, BeginFindFreezingMeasure_k_, Malphax1minusMalpha[IndirectClusterIndex]);
                        }
                    }
                    else
                    {   //  Finally Load Contributions -- Global Cluster Execution
                        if (ParallelClustering.RunningSolution.DiffMalpha_k_Set >= 0)
                            FindDiffMalpha_k_.addapoint(Threadno, IndirectSize, ActiveClustersperPoint, AbsChangeinMalpha_k_);
                        FindC_k_.addapoint(Threadno, IndirectSize, ActiveClustersperPoint, ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha]);
                        FindFreezingMeasure_k_.addapoint(Threadno, IndirectSize, ActiveClustersperPoint, Malphax1minusMalpha);
                    }
                    //  Error is always Global over nodes
                    FindArithmeticError.addapoint(Threadno, LocalArithmeticError);

                }   // End loop over points 
            }); // End parallel Thread loop
            DAVectorUtility.StopSubTimer(2);

            //  Parallel Accumulations. First those using all nodes
            DAVectorUtility.StartSubTimer(3);
            FindArithmeticError.sumoverthreadsandmpi();
            VectorAnnealIterate.ArithmeticError = FindArithmeticError.TotalOr;
            FindMeanClusterCount.sumoverthreadsandmpi();
            VectorAnnealIterate.MeanClusterCount = FindMeanClusterCount.Totalmean;
            FindSinglePoints.sumoverthreadsandmpi();
            VectorAnnealIterate.PointswithClusterCount1 = FindSinglePoints.Total;

            FindUsefulCalcs.sumoverthreadsandmpi();
            FindUselessCalcs.sumoverthreadsandmpi();
            FindIgnoredCalcs.sumoverthreadsandmpi();
            VectorAnnealIterate.LocalUsefulCalcs = FindUsefulCalcs.Total;
            VectorAnnealIterate.LocalUselessCalcs = FindUselessCalcs.Total;
            VectorAnnealIterate.LocalIgnoredCalcs = FindIgnoredCalcs.Total;
            Program.SumUsefulCalcs += VectorAnnealIterate.LocalUsefulCalcs;
            Program.SumUselessCalcs += VectorAnnealIterate.LocalUselessCalcs;
            Program.SumIgnoredCalcs += VectorAnnealIterate.LocalIgnoredCalcs;

            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {   // Distributed Mode Cluster Accumulations
                Find3Components.sumoverthreadsandmpi();

                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int AccumulationPosition = ClusteringSolution.LocalNodeAccPosition[LocalActiveClusterIndex];
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    ParallelClustering.RunningSolution.C_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindC_k_];
                    ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindFreezingMeasure_k_];
                    if (ParallelClustering.RunningSolution.DiffMalpha_k_Set >= 0) // Set DiffMalpha_k_
                        ParallelClustering.RunningSolution.DiffMsummed_k_[RealClusterIndex] = Find3Components.TotalVectorSum[AccumulationPosition][BeginFindDiffMalpha_k_];
                }
            }
            else
            {   // Non distributed Mode
                if (ParallelClustering.RunningSolution.DiffMalpha_k_Set >= 0)
                {   // Set DiffMalpha_k_
                    FindDiffMalpha_k_.sumoverthreadsandmpi();
                    for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                    {
                        int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                        ParallelClustering.RunningSolution.DiffMsummed_k_[RealClusterIndex] = FindDiffMalpha_k_.TotalVectorSum[LocalActiveClusterIndex];
                    }
                }

                FindC_k_.sumoverthreadsandmpi();
                FindFreezingMeasure_k_.sumoverthreadsandmpi();
                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    ParallelClustering.RunningSolution.C_k_[RealClusterIndex] = FindC_k_.TotalVectorSum[LocalActiveClusterIndex];
                    ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] = FindFreezingMeasure_k_.TotalVectorSum[LocalActiveClusterIndex];

                }   // End Code setting  C_k_ FreezingMeasure_k_
            }

            // Broadcast zero size information for global clusters
            bool[] ZeroSizeClusters = new bool[ClusteringSolution.NumberGlobalClusters];
            int CountGlobalClusters = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if ((ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0) || (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 1))
                    continue;
                bool zerosizecluster = false;
                if (RealClusterIndex != ParallelClustering.RunningSolution.SpongeCluster)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;
                ZeroSizeClusters[CountGlobalClusters] = zerosizecluster;
                ++CountGlobalClusters;
            }
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            DAVectorUtility.MPI_communicator.Broadcast<bool>(ref ZeroSizeClusters, 0);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);

            CountGlobalClusters = 0;
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {   // Distributed Mode Cluster Accumulations

                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    bool largeenoughcluster = (ParallelClustering.RunningSolution.C_k_[RealClusterIndex] > Program.CountforCluster_C_ktobezero);
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                    {
                        largeenoughcluster = !ZeroSizeClusters[CountGlobalClusters];
                        ++CountGlobalClusters;
                    }
                    if (largeenoughcluster)
                        ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] = ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] / ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                }
            }
            else
            {   // Non distributed Mode
                for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
                {
                    int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                    bool zerosizecluster = false;
                    if ((ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 0) || (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 1))
                    {
                        zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                        ++CountGlobalClusters;
                    }
                        if(!zerosizecluster)
                            ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] = ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] / ParallelClustering.RunningSolution.C_k_[RealClusterIndex];

                }   // End Code normalizing FreezingMeasure_k_
            }

            DAVectorUtility.StopSubTimer(3);

            DAVectorUtility.StartSubTimer(9);
            //  Indicate status of setting Malpha differences
            ++ParallelClustering.RunningSolution.DiffMalpha_k_Set;

            //  Set new value of P_k_ including case where value for Sponge Cluster is special
            double wgt = 1.0 / DAVectorUtility.PointCount_Global;
            double C_k_Sum01 = 0.0;
            double C_k_Sum2 = 0.0;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0)
                    continue;
                ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                    C_k_Sum2 += ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                else
                    C_k_Sum01 += ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
            }
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                C_k_Sum2 = DAVectorUtility.MPI_communicator.Allreduce<double>(C_k_Sum2, Operation<double>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            VectorAnnealIterate.C_k_Sum = C_k_Sum01 + C_k_Sum2;

            // Correctness check
            double CDiff = C_k_Sum - (double)DAVectorUtility.PointCount_Global;
            if (Math.Abs(CDiff) > 0.5)
            {
                Exception e = DAVectorUtility.SALSAError("C_k_Sum bad " + VectorAnnealIterate.C_k_Sum.ToString("F2") + " Should be " + DAVectorUtility.PointCount_Global.ToString());
                throw (e);
            }

            if ((ParallelClustering.RunningSolution.SpongeCluster != -1) && (Program.SpongePoption == 1))
            {
                wgt = 1.0 - ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster];
                ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster] = Program.SpongePWeight / ParallelClustering.RunningSolution.Ncent_Global;
                wgt = (1.0 - ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster]) / wgt;
                for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                {
                    if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                        continue;
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 0)
                        ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                }
            }
            DAVectorUtility.StopSubTimer(9);
            DAVectorUtility.StartSubTimer(10);

            //  Set Y_k_
            GlobalReductions.FindIndirectVectorDoubleSum[] FindY_k_;
            FindY_k_ = new GlobalReductions.FindIndirectVectorDoubleSum[Program.ParameterVectorDimension];
            DistributedReductions.FindIndirectMultiVectorDoubleSum FindY_k_Component = null;
            int BeginFindY_k_ = -1;
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                FindY_k_Component = new DistributedReductions.FindIndirectMultiVectorDoubleSum();
                BeginFindY_k_ = FindY_k_Component.AddComponents(Program.ParameterVectorDimension);
                FindY_k_Component.NodeInitialize();
            }
            else
            {
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    FindY_k_[VectorIndex] = new GlobalReductions.FindIndirectVectorDoubleSum(DAVectorUtility.ThreadCount,
                        ClusteringSolution.NumberLocalActiveClusters);
                }
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadIndex) =>
            {
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    FindY_k_Component.ThreadInitialize((int)ThreadIndex);
                else
                {
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        FindY_k_[VectorIndex].startthread((int)ThreadIndex);
                }
                int ArraySize = Math.Min(ClusteringSolution.NumberAvailableActiveClusters, ClusteringSolution.MaximumClustersperPoint);
                int[] ActiveClustersperPoint = new int[ArraySize];
                double[][] CenterPositionsperCluster = new double[Program.ParameterVectorDimension][];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    CenterPositionsperCluster[VectorIndex] = new double[ArraySize];
                double[] TerminY_k_ = new double[Program.ParameterVectorDimension];

                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    int ThreadStorePosition = -1;
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        int RealClusterIndex = -2;
                        int RemoteIndex = -1;
                        int ActiveClusterIndex = -1;
                        ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                        ActiveClustersperPoint[IndirectClusterIndex] = ActiveClusterIndex;
                        if (RemoteIndex < 0)
                        {
                            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                                ThreadStorePosition = ClusteringSolution.LocalThreadAccPosition[ActiveClusterIndex][ThreadIndex];
                        }
                        else
                            ThreadStorePosition = DistributedClusteringSolution.StorageforTransportedClusters.TransportedThreadAccPosition[RemoteIndex][ThreadIndex];

                        if ( (ParallelClustering.RunningSolution.SpongeCluster >= 0) && (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster))
                            continue;
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        {
                            TerminY_k_[VectorIndex] = Program.PointPosition[alpha][VectorIndex] * ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            CenterPositionsperCluster[VectorIndex][IndirectClusterIndex] = TerminY_k_[VectorIndex];
                        }
                        if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                        {
                            FindY_k_Component.addapoint((int)ThreadIndex, ThreadStorePosition, BeginFindY_k_, Program.ParameterVectorDimension, TerminY_k_);
                            int CreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        }
                    }
                    if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
                    {
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            FindY_k_[VectorIndex].addapoint((int)ThreadIndex, IndirectSize, ActiveClustersperPoint, CenterPositionsperCluster[VectorIndex]);
                    }
                }
            });  // End Parallel Section calculating center positions

            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                FindY_k_Component.sumoverthreadsandmpi();
            else
            {
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    FindY_k_[VectorIndex].sumoverthreadsandmpi();
            }
                // End section calculating centers

            //  Normalize Centers correctly
            CountGlobalClusters = 0;
            for (int ActiveClusterIndex = 0; ActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; ActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[ActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 0)
                    continue;
                bool zerosizecluster = false;
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)
                {
                    zerosizecluster = ZeroSizeClusters[CountGlobalClusters];
                    ++CountGlobalClusters;
                }
                else if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    zerosizecluster = ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.CountforCluster_C_ktobezero;

                if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                    continue;

                if (!zerosizecluster)
                {   // If zero size cluster, leave Y unchanged
                    if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    {
                        int NodeAccumulationIndex = ClusteringSolution.LocalNodeAccPosition[ActiveClusterIndex];
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        {
                            if (ParallelClustering.RunningSolution.YPreviousSet > -1)
                            {
                                ParallelClustering.RunningSolution.YPrevious_k_i_[RealClusterIndex][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                                ParallelClustering.RunningSolution.YPreviousActuallySet[RealClusterIndex] = true;
                            }
                            ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex] =
                                FindY_k_Component.TotalVectorSum[NodeAccumulationIndex][VectorIndex] / ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                        }
                    }
                    else
                    {
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        {
                            if (ParallelClustering.RunningSolution.YPreviousSet > -1)
                            {
                                ParallelClustering.RunningSolution.YPrevious_k_i_[RealClusterIndex][VectorIndex] = ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex];
                                ParallelClustering.RunningSolution.YPreviousActuallySet[RealClusterIndex] = true;
                            }
                            ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][VectorIndex] = FindY_k_[VectorIndex].TotalVectorSum[ActiveClusterIndex] / ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
                        }
                    }

                    if (Program.SigmaMethod > 1)
                        Program.CalculateSigma(ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex], ref ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex]);
                }
            }   // End section setting Y_k_

            //  Section calculating change in Y
            VectorAnnealIterate.AverageY_k_SquaredChange = 0.0;
            if (ParallelClustering.RunningSolution.YPreviousSet > -1)
            {
                double SumY_k_SquaredChange01 = 0.0;
                double SumY_k_SquaredChange2 = 0.0;
                double NumGlobalSet = 0.0;
                double NumDistributedSet = 0.0;
                for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                {
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] <= 0)
                        continue;
                    if (!ParallelClustering.RunningSolution.YPreviousActuallySet[RealClusterIndex])
                        continue;
                    double SquaredDifference_k_ = DAVectorParallelism.getSquaredScaledDistancebetweenVectors(ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex], ParallelClustering.RunningSolution.YPrevious_k_i_[RealClusterIndex],  
                        ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex]);
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                    {
                        SumY_k_SquaredChange2 += SquaredDifference_k_;
                        NumDistributedSet += 1.0;
                    }
                    else
                    {
                        SumY_k_SquaredChange01 += SquaredDifference_k_;
                        NumGlobalSet += 1.0;
                    }
                }
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                {
                    DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                    SumY_k_SquaredChange2 = DAVectorUtility.MPI_communicator.Allreduce<double>(SumY_k_SquaredChange2, Operation<double>.Add);
                    NumDistributedSet = DAVectorUtility.MPI_communicator.Allreduce<double>(NumDistributedSet, Operation<double>.Add);
                    DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
                }
                if( (NumDistributedSet + NumGlobalSet) > 0.5 )
                    VectorAnnealIterate.AverageY_k_SquaredChange = (SumY_k_SquaredChange01 + SumY_k_SquaredChange2) / (NumDistributedSet + NumGlobalSet);
            }   // End computation of change of vector difference squared 
            ++ParallelClustering.RunningSolution.YPreviousSet;

            ParallelClustering.RunningSolution.CorrelationsSet = false;
            DAVectorUtility.StopSubTimer(10);
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                DistributedClusteringSolution.MinorSynchronizationTransportDistributedClusterCenters();
            return;

        }	// End DAVectorEMIterate()

        // The change of sum(points) Delta(Malpha)/# Points should be less than argument ChangeLimit summed over points/clusters
        //  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
        //  Note divided by number of points here not in calculation
        // Note hitting Program.ConvergenceLoopLimit limit is NOT fatal. Just continues
        public static int convergenceTest(double ChangeLimit, out string ReasonforConvergence)
        {
            ReasonforConvergence = "Not Converged";
            if (VectorAnnealIterate.ArithmeticError)
            {
                VectorAnnealIterate.OnLastleg = true;
                ReasonforConvergence = "Arithmetic Error";
                return 2;
            }
            VectorAnnealIterate.HitConvergenceLoopLimit = false;
            if (ParallelClustering.RunningSolution.DiffMalpha_k_Set < 1)
                return 0;
            double MalphaDiffAvg01 = 0.0;
            double MalphaDiffAvg2 = 0.0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                    MalphaDiffAvg2 += ParallelClustering.RunningSolution.DiffMsummed_k_[RealClusterIndex];
                else
                    MalphaDiffAvg01 += ParallelClustering.RunningSolution.DiffMsummed_k_[RealClusterIndex];
            }
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                MalphaDiffAvg2 = DAVectorUtility.MPI_communicator.Allreduce<double>(MalphaDiffAvg2, Operation<double>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            VectorAnnealIterate.MalphaDiffAvg = (MalphaDiffAvg01 + MalphaDiffAvg2) / DAVectorUtility.PointCount_Global;

            ++Program.NumberMdiffSums;  // Accumulate diagnostic statistics
            Program.MdiffSum += VectorAnnealIterate.MalphaDiffAvg;

            bool notconverged1 = VectorAnnealIterate.MalphaDiffAvg > ChangeLimit;
            if (VectorAnnealIterate.EMIterationStepCount1 < 0)
                VectorAnnealIterate.EMIterationStepCount1 = -2;
            if ((VectorAnnealIterate.EMIterationStepCount1 < 0) && (!notconverged1))
                VectorAnnealIterate.EMIterationStepCount1 = VectorAnnealIterate.EMIterationStepCount;
 
            bool notconverged2 = false;
            double YchangeTest = 0.0;
            ReasonforConvergence = "M cnvgd";
            if ( VectorAnnealIterate.FinalLoop )
            {
                YchangeTest = ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth;
                notconverged2 = VectorAnnealIterate.AverageY_k_SquaredChange > Program.YChangeSquared * YchangeTest;
                if (VectorAnnealIterate.EMIterationStepCount2 < 0)
                    VectorAnnealIterate.EMIterationStepCount2 = -2;
                if ((VectorAnnealIterate.EMIterationStepCount2 < 0) && (!notconverged2))
                    VectorAnnealIterate.EMIterationStepCount2 = VectorAnnealIterate.EMIterationStepCount;
                ReasonforConvergence = "Y and M cnvgd";
            }
            bool notconverged = notconverged1 || notconverged2;
            DAVectorUtility.SynchronizeMPIvariable(ref notconverged);
            if (notconverged) 
            {
                if (notconverged1)
                    ReasonforConvergence = "Converging Malpha " + VectorAnnealIterate.MalphaDiffAvg.ToString("E4") + " " + ChangeLimit.ToString("E4") + " ";
                if (notconverged2)
                    ReasonforConvergence += "Converging Center Change " + VectorAnnealIterate.AverageY_k_SquaredChange.ToString("E4") + " " +  YchangeTest.ToString("E4");
                if (VectorAnnealIterate.EMIterationStepCount > Program.ConvergenceLoopLimit)
                {
                    if (notconverged)
                        DAVectorUtility.SALSAPrint(1, " Too many Iterations " + VectorAnnealIterate.EMIterationStepCount.ToString() + " Warning Message " + ReasonforConvergence);
                    VectorAnnealIterate.HitConvergenceLoopLimit = true;
                }
                else
                    return 0;
            }

            //  "Converged" either due StepCount Limit or Malpha (plus Y) change
            if (VectorAnnealIterate.OnLastleg || VectorAnnealIterate.FinalLoop)
                return 2;

            bool toomanyclusters = (ParallelClustering.RunningSolution.Ncent_ThisNode >= VectorAnnealIterate.ActualMaxNcent)|| (ParallelClustering.RunningSolution.Ncent_Global >= Program.maxNcentTOTAL);
            DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
            toomanyclusters = DAVectorUtility.MPI_communicator.Allreduce<bool>(toomanyclusters, Operation<bool>.LogicalOr);
            DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);

            bool toolowtemperature = ParallelClustering.RunningSolution.Temperature < VectorAnnealIterate.Tmin;
            DAVectorUtility.SynchronizeMPIvariable(ref toolowtemperature);
            if ( toomanyclusters || toolowtemperature || (VectorAnnealIterate.SplitFailures > 1))
            {
                string reason = "";
                if (toomanyclusters)
                    reason = "Too many Clusters";
                if (toolowtemperature)
                    reason += " Temperature Limit Reached";
                if (VectorAnnealIterate.SplitFailures > 1)
                    reason += " Too Many Split Failures";
                DAVectorUtility.SALSAPrint(1, " On last Leg due to " + reason);
                VectorAnnealIterate.OnLastleg = true;
                ReasonforConvergence = " On last Leg due to " + reason;
                return 2;
            }
            return 1;

        }   // End convergenceTest

        public void SaveCurrentTask()
        {   // Save current task that should be split but we need to converge first

            ClusteringSolution.CopySolution(ParallelClustering.RunningSolution, ParallelClustering.SavedSolution);
            return;

        }   // End saving current task that should be split but we need to converge first

        public void RestorePreviousTask()
        {   // Restore previous task that should be split but we needed to converge current cluster configuration first

            ClusteringSolution.CopySolution(ParallelClustering.SavedSolution, ParallelClustering.RunningSolution);
            DistributedClusteringSolution.ManageMajorSynchronization(true);
            return;

        }   // End Restore previous task that should be split but we needed to converge current cluster configuration first

        // Decide if to split based on negative eigenvalue of second derivative matrix
        // MinimumEigenvalue is MINIMUM eigenvalue
        // ClustertoSplit is cluster number of this
        // In distributed mode ONLY distributed clusters split and this is done asynchronously  in each node 
        //  In non distributed mode, clusters are split synchronously
        public bool shouldweSplit()
        {
            //	Calculate Cluster with minimum eigenvalue -- find cluster number and eigenvalue (which could be positive)
            VectorAnnealIterate.Numberthatcanbesplit = 0;
            int LimitonSplits = Math.Min(Program.MaxNumberSplitClusters, VectorAnnealIterate.ActualMaxNcent - ParallelClustering.RunningSolution.Ncent_ThisNode);
            ParallelClustering.RunningSolution.ClustertoSplit = -1;
            bool eigenvaluesexist = false;

            vectorclass vc = new vectorclass();
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                {
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] != 2)
                    {
                        ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 0;
                        continue;
                    }
                }
                else
                {
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] != 1)
                    {
                        ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 0;
                        continue;
                    }
                }
                ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 1;
                if(ParallelClustering.RunningSolution.SplitPriority_k_[RealClusterIndex] == 0)
                    ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 0;

                if ((ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] <= Program.MinimumScaledWidthsquaredtosplit) || (ParallelClustering.RunningSolution.C_k_[RealClusterIndex] <= Program.ToosmalltoSplit) || (RealClusterIndex==ParallelClustering.RunningSolution.SpongeCluster))
                    ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 0;
            }
            if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref ParallelClustering.RunningSolution.Splittable_k_, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            }

            //  Set up eigenvalues
            if (!Program.CalculateEigenvaluesfromMatrix)
                vc.SetAllEigenvaluesIteratively(ParallelClustering.RunningSolution);

            // Only implement Continuous Clustering 
            //  Unlike pairwise call geteigenvalues separately for each cluster
            //   Set Correlation and Change Direction
            ParallelClustering.RunningSolution.SetClusterCorrelations();
            int[] eigenvaluenegative = new int[ParallelClustering.RunningSolution.Ncent_ThisNode];
            int[] eigenvaluestatus = new int[ParallelClustering.RunningSolution.Ncent_ThisNode];
            for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.RunningSolution.Ncent_ThisNode; ClusterToRefine++)
            {
                eigenvaluenegative[ClusterToRefine] = -1;
                eigenvaluestatus[ClusterToRefine] = -1;
                if (ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] == 0)
                {
                    ParallelClustering.RunningSolution.Eigenvalue_k[ClusterToRefine] = 0.0;
                    continue;
                }
                if (Program.CalculateEigenvaluesfromMatrix)
                {
                    double[,] secondderivmatrix = new double[Program.ParameterVectorDimension, Program.ParameterVectorDimension];
                    for (int VectorIndex1 = 0; VectorIndex1 < Program.ParameterVectorDimension; VectorIndex1++)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Program.ParameterVectorDimension; VectorIndex2++)
                        {
                            secondderivmatrix[VectorIndex1, VectorIndex2] = -ParallelClustering.RunningSolution.Correlation_k_i_j[ClusterToRefine][VectorIndex1, VectorIndex2];
                            if (VectorIndex1 == VectorIndex2)
                                secondderivmatrix[VectorIndex1, VectorIndex2] += ParallelClustering.RunningSolution.C_k_[ClusterToRefine] * 0.5;
                        }
                    }
                    vc.getEigenvaluefromMatrix(secondderivmatrix);
                }
                else
                    vc.getEigenvaluefromIteration(ClusterToRefine);
                eigenvaluestatus[ClusterToRefine] =vc.EigenStatus;
                if ( vc.EigenStatus > 0)
                {
                    eigenvaluesexist = true;
                    if (Program.CalculateEigenvaluesfromMatrix)
                    {   // In iterative case, data already stored
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            ParallelClustering.RunningSolution.Eigenvector_k_i[ClusterToRefine][VectorIndex] = vc.Eigenvector[VectorIndex];
                        ParallelClustering.RunningSolution.Eigenvalue_k[ClusterToRefine] = vc.Eigenvalue;
                    }
                    
                    bool EigenvalueNegative = vc.Eigenvalue < 0.0;
                    if (EigenvalueNegative)
                        eigenvaluenegative[ClusterToRefine] = 0;
                    else
                        eigenvaluenegative[ClusterToRefine] = 1;
                }
                if (ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine] >= 1)
                {
                    ++ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine];
                    if (ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine] >= 9)
                        ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine] = -1;

                    if (ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine] < 6)
                    {
                        ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] = 0;
                        continue;
                    }

                }
            }
            if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref eigenvaluenegative, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            }

            for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.RunningSolution.Ncent_ThisNode; ClusterToRefine++)
            {
                if (ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] == 0)
                    continue;
                if (eigenvaluenegative[ClusterToRefine] == -1)
                    continue;
                double CurrentClusterMinimumEigenvalue = ParallelClustering.RunningSolution.Eigenvalue_k[ClusterToRefine];
                bool EigenvalueNegative = false;
                if (eigenvaluenegative[ClusterToRefine] == 0)
                    EigenvalueNegative = true;
                ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] = 2;
                if (EigenvalueNegative)
                    ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] = 3;
            }
            if((DAVectorUtility.MPI_Rank == 0) || ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                for (int ClusterToRefine = 0; ClusterToRefine < ParallelClustering.RunningSolution.Ncent_ThisNode; ClusterToRefine++)
                {
                    if (ParallelClustering.RunningSolution.Splittable_k_[ClusterToRefine] == 0)
                        continue;
                    if (eigenvaluenegative[ClusterToRefine] == -1)
                        continue; 
                    int CurrentClusterPriority = ParallelClustering.RunningSolution.SplitPriority_k_[ClusterToRefine];
                    if (CurrentClusterPriority > 0) CurrentClusterPriority += -10;
                    double CurrentClusterMinimumEigenvalue = ParallelClustering.RunningSolution.Eigenvalue_k[ClusterToRefine];

                    if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
                    {
                        string correlationmessage = "";
                        if( Program.CalculateCorrelationMatrix && Program.ParameterVectorDimension == 2)
                            correlationmessage = " Correls " + ParallelClustering.RunningSolution.Correlation_k_i_j[ClusterToRefine][0, 0].ToString("E4")
                                + " " + ParallelClustering.RunningSolution.Correlation_k_i_j[ClusterToRefine][1, 1].ToString("E4") + " " +
                                ParallelClustering.RunningSolution.Correlation_k_i_j[ClusterToRefine][0, 1].ToString("E4");
                        string eigenvectormessageY = "";
                        string eigenvectormessageDeltaY = "";
                        if( Program.Printeigenvectors)
                        {
                            eigenvectormessageY = " Y ";
                            eigenvectormessageDeltaY = "Delta Y ";
                            for(int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            {
                                eigenvectormessageY += ParallelClustering.RunningSolution.Y_k_i_[ClusterToRefine][VectorIndex].ToString("E4") + " ";
                                eigenvectormessageDeltaY +=  ParallelClustering.RunningSolution.DC_k_DY_k_i_[ClusterToRefine][VectorIndex].ToString("E4") + " ";
                            }
                        }

                        if ((ClusterToRefine < Program.ClusterPrintNumber) || (ClusterToRefine >= (ParallelClustering.RunningSolution.Ncent_ThisNode - Program.ClusterPrintNumber)))
                            DAVectorUtility.SALSAPrint(1, "       Cluster " + ClusterToRefine.ToString() + "(" + ParallelClustering.RunningSolution.LocalCreatedIndex[ClusterToRefine].ToString() + ")"
                                + " Status " + eigenvaluestatus[ClusterToRefine] .ToString() + " Priority " + CurrentClusterPriority.ToString() +
                                " Eigen " + CurrentClusterMinimumEigenvalue.ToString("E4") + " Size " + ParallelClustering.RunningSolution.C_k_[ClusterToRefine].ToString("F1")
                                + eigenvectormessageY + "  " + eigenvectormessageDeltaY + correlationmessage);
                    }

                    bool EigenvalueNegative = false;
                    if (eigenvaluenegative[ClusterToRefine] == 0)
                        EigenvalueNegative = true;
                    if (EigenvalueNegative)
                    {   // Candidate for Split List

                        double test1 = 2.0 * ( (double) DAVectorUtility.PointCount_Global) / ((double) ParallelClustering.RunningSolution.Ncent_ThisNode);
                        if (ParallelClustering.RunningSolution.C_k_[ClusterToRefine] >= test1)
                            CurrentClusterPriority += 20;
                        if (ParallelClustering.RunningSolution.C_k_[ClusterToRefine] >= (3.0 * test1) )
                            CurrentClusterPriority += 20;

                        if (VectorAnnealIterate.Numberthatcanbesplit == 0)     // Initialize Split List
                        {
                            VectorAnnealIterate.Numberthatcanbesplit = 1;
                            VectorAnnealIterate.EigsofClusterstoSplit[0] = CurrentClusterMinimumEigenvalue;
                            VectorAnnealIterate.ListofClusterstoSplit[0] = ClusterToRefine;
                            VectorAnnealIterate.PrioritiesofClusterstoSplit[0] = CurrentClusterPriority;
                        }
                        else     // Add to Split List
                        {
                            int position = VectorAnnealIterate.Numberthatcanbesplit;
                            for (int positionloop = 0; positionloop < VectorAnnealIterate.Numberthatcanbesplit; positionloop++)
                            {
                                if (CurrentClusterPriority < VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop])
                                    continue;
                                if (CurrentClusterPriority == VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop])
                                {
                                    bool eigtest = (CurrentClusterMinimumEigenvalue >= VectorAnnealIterate.EigsofClusterstoSplit[positionloop]);
                                    if (eigtest)
                                        continue;
                                }
                                position = positionloop;
                                break;
                            }
                            if (position >= LimitonSplits)
                                continue;
                            for (int positionloop = VectorAnnealIterate.Numberthatcanbesplit - 1; positionloop >= position; positionloop--)
                            {
                                if (positionloop == (LimitonSplits - 1))
                                    continue;
                                VectorAnnealIterate.EigsofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.EigsofClusterstoSplit[positionloop];
                                VectorAnnealIterate.ListofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.ListofClusterstoSplit[positionloop];
                                VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop + 1] = VectorAnnealIterate.PrioritiesofClusterstoSplit[positionloop];
                            }
                            VectorAnnealIterate.Numberthatcanbesplit = Math.Min(VectorAnnealIterate.Numberthatcanbesplit + 1, LimitonSplits);
                            VectorAnnealIterate.EigsofClusterstoSplit[position] = CurrentClusterMinimumEigenvalue;
                            VectorAnnealIterate.ListofClusterstoSplit[position] = ClusterToRefine;
                            VectorAnnealIterate.PrioritiesofClusterstoSplit[position] = CurrentClusterPriority;
                        }
                    }
                }
            }
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                eigenvaluesexist = DAVectorUtility.MPI_communicator.Allreduce<bool>(eigenvaluesexist, Operation<bool>.LogicalOr);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            else
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref VectorAnnealIterate.Numberthatcanbesplit, 0);
                DAVectorUtility.MPI_communicator.Broadcast<double>(ref VectorAnnealIterate.EigsofClusterstoSplit, 0);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref VectorAnnealIterate.ListofClusterstoSplit, 0);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref VectorAnnealIterate.PrioritiesofClusterstoSplit, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            }
            if (!eigenvaluesexist)
            {
                VectorAnnealIterate.SplitFailures++;
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref VectorAnnealIterate.SplitFailures, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                return false;
            }
            if (VectorAnnealIterate.Numberthatcanbesplit == 0)
                return false;
            ParallelClustering.RunningSolution.ClustertoSplit = VectorAnnealIterate.ListofClusterstoSplit[0];
            return true;

        }   // End shouldweSplit

        //  Do a split on the identified cluster
        //  New clusters stored in last cluster and old position
        //  Cluster count is incremented by one
        //  ClustertoSplit is cluster number to split
        public void dothesplit()
        {   // This just gets Malpha_k_ P_k_ and Y_k_

            if (ParallelClustering.RunningSolution.ClustertoSplit < 0)
                return;
            int NewCenterIndex = ParallelClustering.RunningSolution.Ncent_ThisNode;
            int OldCenterIndex = ParallelClustering.RunningSolution.ClustertoSplit;
            ++ClusteringSolution.ClustersSplit;

            //  Calculate Perturbed center positions
            // First Estimate Scale as minumum of one that produces .05 change in population of new centers and one that shifts 0.1 times cluster width
            double ambitiousscale = 0.0;
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                ambitiousscale += ParallelClustering.RunningSolution.Eigenvector_k_i[OldCenterIndex][VectorIndex] * ParallelClustering.RunningSolution.DC_k_DY_k_i_[OldCenterIndex][VectorIndex];
            ambitiousscale = Math.Abs(0.025 * ParallelClustering.RunningSolution.C_k_[OldCenterIndex] / ambitiousscale);
            double actualscale = Math.Min(ambitiousscale, 0.2 * Math.Sqrt(ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[OldCenterIndex]));
            actualscale = Math.Max(actualscale, 0.05 * Math.Sqrt(ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[OldCenterIndex]));

            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                double tmp1 = ParallelClustering.RunningSolution.Y_k_i_[OldCenterIndex][VectorIndex];
                double tmp2 = ParallelClustering.RunningSolution.Eigenvector_k_i[OldCenterIndex][VectorIndex] * actualscale
                    * Math.Sqrt(ParallelClustering.RunningSolution.Sigma_k_i_[OldCenterIndex][VectorIndex]);
                ParallelClustering.RunningSolution.Y_k_i_[OldCenterIndex][VectorIndex] = tmp1 + tmp2;
                ParallelClustering.RunningSolution.Y_k_i_[NewCenterIndex][VectorIndex] = tmp1 - tmp2;
            }
            ParallelClustering.RunningSolution.YPreviousActuallySet[OldCenterIndex] = false;
            ParallelClustering.RunningSolution.YPreviousActuallySet[NewCenterIndex] = false;
            Program.CalculateSigma(ParallelClustering.RunningSolution.Y_k_i_[NewCenterIndex], ref ParallelClustering.RunningSolution.Sigma_k_i_[NewCenterIndex]);
            ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;

            ParallelClustering.RunningSolution.Splittable_k_[OldCenterIndex] = -1;
            ParallelClustering.RunningSolution.Splittable_k_[NewCenterIndex] = -1;
            ParallelClustering.RunningSolution.SplitPriority_k_[OldCenterIndex] = 2;
            ParallelClustering.RunningSolution.SplitPriority_k_[NewCenterIndex] = 2;
            ParallelClustering.RunningSolution.LocalStatus[NewCenterIndex] = ParallelClustering.RunningSolution.LocalStatus[OldCenterIndex];
            ParallelClustering.RunningSolution.LocalSplitCreatedIndex[NewCenterIndex] = 0;
            int CreatedIndex_child = ClusteringSolution.SetCreatedIndex(NewCenterIndex);
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                ParallelClustering.RunningSolution.LocalSplitCreatedIndex[OldCenterIndex] = -1 - CreatedIndex_child;
                ParallelClustering.RunningSolution.LocalSplitCreatedIndex[NewCenterIndex] = +1 + ParallelClustering.RunningSolution.LocalCreatedIndex[OldCenterIndex];
            }
            ParallelClustering.RunningSolution.P_k_[OldCenterIndex] = 0.5 * ParallelClustering.RunningSolution.P_k_[OldCenterIndex];
            ParallelClustering.RunningSolution.P_k_[NewCenterIndex] = ParallelClustering.RunningSolution.P_k_[OldCenterIndex];
            ParallelClustering.RunningSolution.C_k_[OldCenterIndex] = 0.5 * ParallelClustering.RunningSolution.C_k_[OldCenterIndex];
            ParallelClustering.RunningSolution.C_k_[NewCenterIndex] = ParallelClustering.RunningSolution.C_k_[OldCenterIndex];

            // Increase number of Clusters
            ParallelClustering.RunningSolution.Ncent_ThisNode++;

            //  Exit if Distributed
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                return;

            ParallelClustering.RunningSolution.SetActiveClusters();
            if (Program.UseTriangleInequality_DA > 0)
                DATriangleInequality.SplitCenter(OldCenterIndex, ParallelClustering.RunningSolution.LocalCreatedIndex[OldCenterIndex], NewCenterIndex, CreatedIndex_child);

            int OldActiveClusterIndex = ClusteringSolution.ActiveClusterIndices[OldCenterIndex];

            //  Parallel Section Splitting Malpha_k_ for non distributed mode
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {   // Perturb Malpha_k_
                    int IndirectClusterIndex = ParallelClustering.RunningSolution.MapClusterToIndirect(alpha, OldActiveClusterIndex );
                    if (IndirectClusterIndex < 0)
                        continue;

                    double newvalueofM = ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] * 0.5;

                    ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] = newvalueofM;
                    int NewIndirectClusterIndex = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][NewIndirectClusterIndex] = CreatedIndex_child;
                    ++ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][NewIndirectClusterIndex] = newvalueofM;
                }
            }); //  End Parallel Section Splitting Cluster

        }   // End dothesplit

            public static bool CheckNumberofClustersTooBig()
        {
            bool changed = false;
            GlobalReductions.FindIntSum NumberChanged = new GlobalReductions.FindIntSum(DAVectorUtility.ThreadCount);

            //  Check Number of Clusters for each point. This could be stricter
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    if (ParallelClustering.RunningSolution.NumClusters_alpha_[alpha] <= ClusteringSolution.TargetClustersperPoint)
                        continue;
                    ParallelClustering.RunningSolution.SetClustersforaPoint(alpha);
                    NumberChanged.addapoint(ThreadNo, 1);
                }

            }); //  End Section Setting Cluster numbers for points
            NumberChanged.sumoverthreadsandmpi();
            
            if (NumberChanged.TotalInt != 0)
            {
                if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0) 
                    DAVectorUtility.SALSAPrint(1, " Points Changed " + NumberChanged.TotalInt.ToString());
                changed = true;
            }
            return changed;

        }   // End CheckNumberofClustersTooBig()

        public static void CompleteCleanUp()
        {
            DAVectorUtility.SALSAPrint(0, "Complete Clean Up T " + ParallelClustering.RunningSolution.Temperature);
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    ParallelClustering.RunningSolution.
                        SetClustersforaPoint(alpha);
                }

            }); //  End Section Setting Cluster numbers for points

            VectorAnnealIterate.HammyNotSet = true;
            ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;

        }   // End CompleteCleanUp()

        //	Find initial Temperature Cluster Center and other initial parameters
        public static void InitializeSolution(ClusteringSolution StartSolution)
        {
            
            //  Deterministic Annealing
            double[] Initial_Y = new double[Program.ParameterVectorDimension];
            int SpongePosition = StartSolution.SpongeCluster;
            double InitialM = 1.0;
            StartSolution.P_k_[0] = 1.0;
            int FirstRealCluster = 0;
            if (SpongePosition >= 0)
            {
                InitialM = 0.5;
                StartSolution.P_k_[0] = 0.5;
                StartSolution.P_k_[1] = 0.5;
                FirstRealCluster = 1;
            }
            StartSolution.SplitPriority_k_[FirstRealCluster] = -1;
            StartSolution.Splittable_k_[FirstRealCluster] = 0;
            StartSolution.LocalSplitCreatedIndex[FirstRealCluster] = 0;
            StartSolution.LocalStatus[FirstRealCluster] = 1;
            int CreatedIndex = ClusteringSolution.SetCreatedIndex(FirstRealCluster);

            GlobalReductions.FindArrayMean SystemCoG = new GlobalReductions.FindArrayMean(DAVectorUtility.ThreadCount, Program.ParameterVectorDimension);
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                if( SpongePosition >=0)
                    StartSolution.Y_k_i_[SpongePosition][VectorIndex] = 0.0;    // Set "center" of sponge to Origin
            }

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    SystemCoG.addapoint(ThreadIndex, Program.PointPosition[alpha] );

                    StartSolution.Map_alpha_PointertoCreatedIndex[alpha][FirstRealCluster] = CreatedIndex;
                    StartSolution.M_alpha_kpointer_[alpha][FirstRealCluster] = InitialM;
                    if (SpongePosition >= 0)
                    {
                        StartSolution.M_alpha_kpointer_[alpha][SpongePosition] = InitialM;
                        StartSolution.Map_alpha_PointertoCreatedIndex[alpha][SpongePosition] = StartSolution.LocalCreatedIndex[SpongePosition];
                    }
                }
            });  // End Parallel Section calculating initial center positions

            SystemCoG.sumoverthreadsandmpi();
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                Initial_Y[VectorIndex] = SystemCoG.Totalmean[VectorIndex];
                StartSolution.Y_k_i_[FirstRealCluster][VectorIndex] = Initial_Y[VectorIndex];
            }
            Program.CalculateSigma(Initial_Y, ref StartSolution.Sigma_k_i_[FirstRealCluster]);
            
            GlobalReductions.FindDoubleMean InitialAverages = new GlobalReductions.FindDoubleMean(DAVectorUtility.ThreadCount);

            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadIndex) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadIndex];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadIndex] - DAVectorUtility.PointStart_Process;

                for (int LocalPointIndex = beginpoint; LocalPointIndex < indexlen + beginpoint; LocalPointIndex++)
                    InitialAverages.addapoint(ThreadIndex, DAVectorParallelism.getSquaredScaledDistancePointActiveCluster(LocalPointIndex, FirstRealCluster, StartSolution));
            });  // End Parallel Section calculating average distance
            InitialAverages.sumoverthreadsandmpi();

            //  Estimate of Initial Temperature is average scaled squared distance
            //  Fudge factor of 1.5 over estimated critical temperature for first split
            StartSolution.Temperature = 1.5 * InitialAverages.Totalmean;
            if (Program.Tminimum > 0.0)
                VectorAnnealIterate.Tmin = Program.Tminimum;
            else
                VectorAnnealIterate.Tmin = StartSolution.Temperature / Math.Abs(Program.Tminimum);
            Program.ActualStartTemperature = StartSolution.Temperature; // For Output
            Program.TargetEndTemperature = VectorAnnealIterate.Tmin;   // For Output
            DAVectorUtility.SALSAPrint(1, "Points " + DAVectorUtility.PointCount_Global.ToString() + " Initial Temperature " + StartSolution.Temperature.ToString("F4") + " Minimum " + VectorAnnealIterate.Tmin.ToString("F4") );
            
            StartSolution.ActualCoolingFactor = Program.InitialCoolingFactor;
            
            StartSolution.PairwiseHammy = 0.0;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                StartSolution.PairwiseHammy += 0.5 * ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[RealClusterIndex] * ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
            return;

        }	// End InitializeSolution

        //  This assumes essentially that Sponge Used   and is cluster 0 in both files
        //  Assumes that second run processes sponge from first
        public static void Restart()
        {
            //  Reset Input File selection
            Program.SelectedInputLabel = Program.RestartSelectedInputLabel;
            Program.InputFileType = 1;

            // Read old Center Positions
            Program.Replicate = 1;
            ParallelClustering.RunningSolution.Ncent_Global = 0;
            if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                ParallelClustering.RunningSolution.Ncent_Global = 1;
            --ParallelClustering.RunningSolution.Ncent_ThisNode;
            DAVectorReadData.ReadDataFromFile(Program.RestartClusterFile, 2);
            int InitialClusterCount = ParallelClustering.RunningSolution.Ncent_Global;
            int ExtraClusterCount = 0;
            if (Program._configurationManager.DAVectorSpongeSection.LabelFile.Length > 0)
            {
                DAVectorReadData.ReadDataFromFile(Program._configurationManager.DAVectorSpongeSection.LabelFile, 1);
                DAVectorReadData.ReadDataFromFile(Program._configurationManager.DAVectorSpongeSection.LabelFile, 2);
                ExtraClusterCount = ParallelClustering.RunningSolution.Ncent_Global - InitialClusterCount;
            }
            ParallelClustering.RunningSolution.Ncent_ThisNode = ParallelClustering.RunningSolution.Ncent_Global;
            Program.InitialNcent = ParallelClustering.RunningSolution.Ncent_Global;

            // Set up Clusters including CreatedIndex. This is in Global not Distributed mode
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = 1.0/(  (double) ParallelClustering.RunningSolution.Ncent_ThisNode);
                ParallelClustering.RunningSolution.FreezingMeasure_k_[RealClusterIndex] = 0.0;
                ParallelClustering.RunningSolution.SplitPriority_k_[RealClusterIndex] = -1;
                ParallelClustering.RunningSolution.Splittable_k_[RealClusterIndex] = 0;
                ParallelClustering.RunningSolution.LocalSplitCreatedIndex[RealClusterIndex] = 0;
                ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] = 0;
                if (RealClusterIndex != ParallelClustering.RunningSolution.SpongeCluster)
                {
                    ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] = 1;
                    Program.CalculateSigma(ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex], ref ParallelClustering.RunningSolution.Sigma_k_i_[RealClusterIndex]);
                }
                int CreatedIndex = ClusteringSolution.SetCreatedIndex(RealClusterIndex); // Sets LocalCreatedIndex
                ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex] = 0;
            }

            ParallelClustering.RunningSolution.Temperature = Program.RestartTemperature;
            if (Program.Tminimum > 0.0)
                VectorAnnealIterate.Tmin = Program.Tminimum;
            else
                VectorAnnealIterate.Tmin = ParallelClustering.RunningSolution.Temperature / Math.Abs(Program.Tminimum);
            Program.ActualStartTemperature = ParallelClustering.RunningSolution.Temperature; // For Output
            Program.TargetEndTemperature = VectorAnnealIterate.Tmin;   // For Output

            ParallelClustering.RunningSolution.SetActiveClusters();

            //  Set initial clusters for every point!
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {

                    int NumClustersperPoint = 0;
                    int AssignedCluster = Program.PointLabel[alpha];
                    if (AssignedCluster >= ParallelClustering.RunningSolution.Ncent_Global)
                    {
                        Exception e = DAVectorUtility.SALSAError(" Illegal Cluster Number " + AssignedCluster.ToString() + " At Point "
                            + (alpha + DAVectorUtility.PointStart_Process).ToString() + " Total " + ParallelClustering.RunningSolution.Ncent_Global.ToString() );
                        throw (e);
                    }
                    if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                    {   // Always have a sponge entry for each point
                        NumClustersperPoint = 1;
                        ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][0] = ParallelClustering.RunningSolution.LocalCreatedIndex[ParallelClustering.RunningSolution.SpongeCluster];
                        ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][0] = 1.0;
                    }
                    if ( (AssignedCluster != 0) || (ParallelClustering.RunningSolution.SpongeCluster < 0) )
                    {
                        ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][NumClustersperPoint] = ParallelClustering.RunningSolution.LocalCreatedIndex[AssignedCluster];
                        ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][NumClustersperPoint] = 1.0;
                        if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                            ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][0] = 0.0;
                        ++NumClustersperPoint;
                    }
                    ParallelClustering.RunningSolution.NumClusters_alpha_[alpha] = NumClustersperPoint;
                }
            }); //  End Parallel Section

            ParallelClustering.RunningSolution.FindOccupationCounts();
            double wgt = 1.0 / DAVectorUtility.PointCount_Global;
            for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
            {
                if (ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex] == 0)
                {
                    // Exception e = DAVectorUtility.SALSAError(" No points for Cluster Number " + RealClusterIndex.ToString() + " Total " + ParallelClustering.RunningDAVectorSolt.Ncent_Global.ToString());
                    // throw (e);
                    DAVectorUtility.SALSAPrint( 0, " No points for Cluster Number " + RealClusterIndex.ToString() + " Total " + ParallelClustering.RunningSolution.Ncent_Global.ToString());
                }
                ParallelClustering.RunningSolution.C_k_[RealClusterIndex] = (double)ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex];
                ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.RunningSolution.C_k_[RealClusterIndex];
            }
            string message = "Restart Set up Initial Clusters " + InitialClusterCount.ToString() + " Extra Clusters " + ExtraClusterCount.ToString();
            if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                message += " Sponge Count " + ParallelClustering.RunningSolution.OccupationCounts_k_[ParallelClustering.RunningSolution.SpongeCluster].ToString();
            else
                message += " No Sponge";
            DAVectorUtility.SALSAPrint(1, message);

            if ((ParallelClustering.RunningSolution.SpongeCluster != -1) && (Program.SpongePoption == 1))
            {
                wgt = 1.0 - ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster];
                ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster] = Program.SpongePWeight / ParallelClustering.RunningSolution.Ncent_Global;
                wgt = (1.0 - ParallelClustering.RunningSolution.P_k_[ParallelClustering.RunningSolution.SpongeCluster]) / wgt;
                for (int RealClusterIndex = 0; RealClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; RealClusterIndex++)
                {
                    if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                        continue;
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] > 0)
                        ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = wgt * ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                }
            }
            ParallelClustering.RunningSolution.SetClusterWidths();
            double AverageClusterwidth = ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth;
            double[] oldwidth = new double[Program.ParameterVectorDimension];
            if (Program.CalculateIndividualWidths)
            {
                for(int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    oldwidth[VectorIndex] = ParallelClustering.RunningSolution.AverageWidth[VectorIndex];
            }
            Program.ActualSpongeTemperature = ParallelClustering.RunningSolution.Temperature;
            Program.ActualWidth = AverageClusterwidth;
            DAVectorUtility.InterimTiming();
            Program.TimeatSponge = DAVectorUtility.HPDuration;

            VectorAnnealIterate.DistributeNextTime = true;
            ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;

            VectorAnnealIterate.CompleteCleanUp();
            ParallelClustering.RunningSolution.SetClusterWidths();
            string widthmessage = "";
            if (Program.CalculateIndividualWidths)
            {
                for(int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    widthmessage += VectorIndex.ToString() + " Component Widths " + oldwidth[VectorIndex].ToString("E4") + " " 
                        + ParallelClustering.RunningSolution.AverageWidth[VectorIndex].ToString("E4") + " ";
            }
            else
            {
                widthmessage = "Total width " + ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth.ToString("E4");
            }
            DAVectorUtility.SALSAPrint(1,widthmessage);

        }   // End Restart


        public static bool CheckCloseClusters(int LocalActiveClusterIndex1)
        {
            int RealClusterIndex1 = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex1];
            for (int LocalActiveClusterIndex2 = 0; LocalActiveClusterIndex2 < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex2++)
            {
                if(LocalActiveClusterIndex2 == LocalActiveClusterIndex1)
                    continue;
                int RealClusterIndex2 = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex2];
                if (RealClusterIndex2 == ParallelClustering.RunningSolution.SpongeCluster)
                    continue;
                int test = ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex1] - ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex2];
                if (test > 0)
                    continue;
                if ( (test == 0) && (LocalActiveClusterIndex1 < LocalActiveClusterIndex2) )
                    continue;
                double TestDistce = DAVectorParallelism.getSquaredScaledDistanceTweenActiveClusters(LocalActiveClusterIndex1, LocalActiveClusterIndex2, ParallelClustering.RunningSolution);
                if (TestDistce < Program.ScaledSquaredDistanceatClosenessTest)
                    return true;
            }
            return false;

        }   // End CheckCloseClusters

        // Return False if current final solution has a problem
        //  Return True if there is no problem or no way of fixing problem
        public static bool CheckValidSolution(bool UseBestSolution)
        {
            if (UseBestSolution)
            {
                if (!ParallelClustering.BestSolution.SolutionSet)
                {
                    DAVectorUtility.SALSAPrint(1, "Solution Not changed as no best solution");
                }
                if ((ParallelClustering.BestSolution.Ncent_ThisNode != ParallelClustering.RunningSolution.Ncent_ThisNode) || (ParallelClustering.BestSolution.IterationSetAt != ParallelClustering.RunningSolution.IterationSetAt))
                    DAVectorUtility.SALSAPrint(1, " Best Solution Not Used to restart as Number of Centers " + ParallelClustering.BestSolution.Ncent_ThisNode + " Different from Running Solution with "
                        + ParallelClustering.RunningSolution.Ncent_ThisNode + " Best Iteration# " + ParallelClustering.BestSolution.IterationSetAt + " Current " + ParallelClustering.RunningSolution.IterationSetAt);
                else
                {
                    bool hammychange = ParallelClustering.RunningSolution.PairwiseHammy > ParallelClustering.BestSolution.PairwiseHammy;
                    DAVectorUtility.SynchronizeMPIvariable(ref hammychange);
                    if (hammychange)
                    {
                        ClusteringSolution.CopySolution(ParallelClustering.BestSolution, ParallelClustering.RunningSolution);
                        DistributedClusteringSolution.ManageMajorSynchronization(true);
                        DAVectorUtility.SALSAPrint(1, "Best Solution Used");
                        ParallelClustering.RunningSolution.OldHammy = 0.0;
                        ParallelClustering.RunningSolution.PairwiseHammy = 0.0;
                        ParallelClustering.RunningSolution.Eigenvectorset = false;
                        ParallelClustering.RunningSolution.ClustertoSplit = -1;
                    }
                }
            }

            bool CheckCloseness = false;
            if (ParallelClustering.RunningSolution.Temperature < Program.TemperatureforClosenessTest)
            {
                if ((VectorAnnealIterate.TemperatureatLastCloseTest < 0.0) || (ParallelClustering.RunningSolution.Temperature < VectorAnnealIterate.TemperatureatLastCloseTest - 0.25))
                    CheckCloseness = true;
            }
            DAVectorUtility.SynchronizeMPIvariable(ref CheckCloseness);
            if(CheckCloseness)
                VectorAnnealIterate.TemperatureatLastCloseTest = ParallelClustering.RunningSolution.Temperature;

            bool[] RemoveCluster = new bool[ParallelClustering.RunningSolution.Ncent_ThisNode];
            int Numbertoosmall01 = 0;
            int Numbertoosmall2 = 0;
            double ProbabilitySum01 = 0.0;
            double ProbabilitySum2 = 0.0;
            string documentit = "";
            ParallelClustering.RunningSolution.FindOccupationCounts();
            double C_k_Test = Program.MinimumCountforCluster_C_k;
            if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                C_k_Test = Program.MinimumCountforCluster_C_kwithSponge;
            int NumberSmall_C = 0;
            int NumberSmall_OccCount = 0;
            int NumberClose = 0;
            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                RemoveCluster[RealClusterIndex] = false;
                if (RealClusterIndex == ParallelClustering.RunningSolution.SpongeCluster)
                {
                    if(ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                        ProbabilitySum2 += ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                    else
                        ProbabilitySum01 += ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                    continue;
                }
                if (!(ParallelClustering.RunningSolution.DistributedExecutionMode && (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] < 2)) )
                {
                    if (ParallelClustering.RunningSolution.C_k_[RealClusterIndex] < C_k_Test)
                    {
                        RemoveCluster[RealClusterIndex] = true;
                        ++NumberSmall_C;
                    }
                    else if (ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex] < Program.MinimumCountforCluster_Points)
                    {
                        RemoveCluster[RealClusterIndex] = true;
                        ++NumberSmall_OccCount;
                    }
                    if (CheckCloseness && (!RemoveCluster[RealClusterIndex]))
                    {
                        RemoveCluster[RealClusterIndex] = CheckCloseClusters(LocalActiveClusterIndex);
                        if(RemoveCluster[RealClusterIndex])
                            ++NumberClose;
                    }
                }
            }

            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                NumberSmall_C = DAVectorUtility.MPI_communicator.Allreduce<int>(NumberSmall_C, Operation<int>.Add);
                NumberSmall_OccCount = DAVectorUtility.MPI_communicator.Allreduce<int>(NumberSmall_OccCount, Operation<int>.Add);
                NumberClose = DAVectorUtility.MPI_communicator.Allreduce<int>(NumberClose, Operation<int>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            else
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<bool>(ref RemoveCluster, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
            }
            Program.TotalClustersDeleted_CSmall += NumberSmall_C;
            Program.TotalClustersDeleted_OccCount += NumberSmall_OccCount;
            Program.TotalClustersDeleted_Close += NumberClose;

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex]; 
                if (RemoveCluster[RealClusterIndex])
                {
                    if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                        ++Numbertoosmall2;
                    else
                        ++Numbertoosmall01;
                    documentit += RealClusterIndex.ToString() +  "(" + ParallelClustering.RunningSolution.LocalCreatedIndex[RealClusterIndex].ToString() + ") C_k_ " + ParallelClustering.RunningSolution.C_k_[RealClusterIndex].ToString("F1")
                            + " Points " + ParallelClustering.RunningSolution.OccupationCounts_k_[RealClusterIndex].ToString();
                    if (Program.ParameterVectorDimension == 2)
                        documentit +=
                            " Y " + ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][0].ToString("F3") + " " + ParallelClustering.RunningSolution.Y_k_i_[RealClusterIndex][1].ToString("F3") + " * ";
                    else
                        documentit += " * ";
                    continue;
                }
                if (ParallelClustering.RunningSolution.LocalStatus[RealClusterIndex] == 2)
                    ProbabilitySum2 += ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                else
                    ProbabilitySum01 += ParallelClustering.RunningSolution.P_k_[RealClusterIndex];
                continue;
            }
            
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIREDUCETiming3);
                Numbertoosmall2 = DAVectorUtility.MPI_communicator.Allreduce<int>(Numbertoosmall2, Operation<int>.Add);
                ProbabilitySum2 = DAVectorUtility.MPI_communicator.Allreduce<double>(ProbabilitySum2, Operation<double>.Add);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIREDUCETiming3);
            }
            double ProbabilitySum = ProbabilitySum01 + ProbabilitySum2;
            int GlobalNumberTooSmall = Numbertoosmall01 + Numbertoosmall2;
            if (GlobalNumberTooSmall == 0)
                return true;

            //  Need to change solution as one or more clusters too small
            if (VectorAnnealIterate.EMIterationCount % Program.PrintInterval == 0)
            {
                if (ParallelClustering.RunningSolution.DistributedExecutionMode)
                    DAVectorUtility.SALSASyncPrint(1, GlobalNumberTooSmall.ToString() + " Clusters Too Small ", documentit);
                else
                    DAVectorUtility.SALSAPrint(1, GlobalNumberTooSmall.ToString() + " Clusters Too Small " + documentit);
            }

            for (int LocalActiveClusterIndex = 0; LocalActiveClusterIndex < ClusteringSolution.NumberLocalActiveClusters; LocalActiveClusterIndex++)
            {
                int RealClusterIndex = ClusteringSolution.RealClusterIndices[LocalActiveClusterIndex];
                if (!RemoveCluster[RealClusterIndex])
                {
                    if (Program.ContinuousClustering)
                        ParallelClustering.RunningSolution.P_k_[RealClusterIndex] = ParallelClustering.RunningSolution.P_k_[RealClusterIndex]/ProbabilitySum;
                    continue;
                }
            }
            if (!ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                int oldNcent = ParallelClustering.RunningSolution.Ncent_ThisNode;
                for (int clusterindex = oldNcent - 1; clusterindex >= 0; clusterindex--)
                {
                    if (RemoveCluster[clusterindex])
                        ParallelClustering.RunningSolution.RemoveCluster(clusterindex); // This both changes cluster count and shifts up clusters
                }
            }
            else
            {
                for (int clusterindex = 0; clusterindex < ParallelClustering.RunningSolution.Ncent_ThisNode; clusterindex++)
                {
                    if (RemoveCluster[clusterindex])
                        ParallelClustering.RunningSolution.LocalStatus[clusterindex] = -1;
                }
            }
            DistributedClusteringSolution.ManageMajorSynchronization(true);
            ClusteringSolution.CopySolution(ParallelClustering.RunningSolution, ParallelClustering.BestSolution);
            return false;

        }   // End CheckValidSolution

        public static bool AddSpongeCluster()
        {                
            if (ParallelClustering.RunningSolution.SpongeCluster >= 0) return false;
            int CurrentMax = ParallelClustering.RunningSolution.Ncent_ThisNode;
            if (CurrentMax >= Program.maxNcentperNode)
                return false;

            double AverageClusterwidth = ParallelClustering.RunningSolution.TotaloverVectorIndicesAverageWidth;

            bool Addasponge = (Program.CreateSpongeScaledSquaredWidth > 0.0) && (AverageClusterwidth <= Program.CreateSpongeScaledSquaredWidth);
            if (!Addasponge)
                Addasponge = (Program.SpongeTemperature1 > 0.0) && (ParallelClustering.RunningSolution.Temperature < Program.SpongeTemperature1);
            DAVectorUtility.SynchronizeMPIvariable(ref Addasponge);
            if (!Addasponge)
                return false;
            VectorAnnealIterate.OutputClusteringResults("StartSponge");

            //  Note Cluster numbers different in each node if distributed execution. Created Index identical across nodes
            ParallelClustering.RunningSolution.SpongeCluster = CurrentMax;
            int SpongeCreatedIndex = 0;
            if (ParallelClustering.RunningSolution.DistributedExecutionMode)
            {
                if(DAVectorUtility.MPI_Rank == 0)
                    SpongeCreatedIndex = ClusteringSolution.SetCreatedIndex(CurrentMax);
                DAVectorUtility.StartSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                DAVectorUtility.MPI_communicator.Broadcast<int>(ref SpongeCreatedIndex, 0);
                DAVectorUtility.StopSubTimer(DAVectorUtility.MPIBROADCASTTiming);
                if (DAVectorUtility.MPI_Rank != 0)
                {
                    ClusteringSolution.UniversalMapping[SpongeCreatedIndex] = new ClusterIndirection(ClusteringSolution.CurrentIteration, 1 + CurrentMax);
                    ParallelClustering.RunningSolution.LocalCreatedIndex[CurrentMax] = SpongeCreatedIndex;
                }
            }
            else
                SpongeCreatedIndex = ClusteringSolution.SetCreatedIndex(CurrentMax);
            
            Program.ActualSpongeTemperature = ParallelClustering.RunningSolution.Temperature;
            Program.ActualWidth = AverageClusterwidth;
            DAVectorUtility.InterimTiming();
            Program.TimeatSponge = DAVectorUtility.HPDuration;
            DAVectorUtility.SALSAPrint(1, "Sponge added with Created Index " + SpongeCreatedIndex.ToString() + " Time " + Program.TimeatSponge.ToString("F0"));

            //  Add Sponge to every point!
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                double[] workingvector = new double[Program.ParameterVectorDimension];
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {   
                    int NewIndirectClusterIndex = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][NewIndirectClusterIndex] = SpongeCreatedIndex;
                    ++ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][NewIndirectClusterIndex] = 0.0;
                }
            }); //  End Parallel Section

            //  Set up Sponge Cluster
            ParallelClustering.RunningSolution.Splittable_k_[CurrentMax] = 0;
            ParallelClustering.RunningSolution.SplitPriority_k_[CurrentMax] = 0;
            ParallelClustering.RunningSolution.C_k_[CurrentMax] = 0.0;
            ParallelClustering.RunningSolution.ClusterScaledSquaredWidth_k_[CurrentMax] = 0.0;
            ParallelClustering.RunningSolution.FreezingMeasure_k_[CurrentMax] = 0.0;
            ParallelClustering.RunningSolution.LocalStatus[CurrentMax] = 0;
            ParallelClustering.RunningSolution.LocalSplitCreatedIndex[CurrentMax] = 0;

            ParallelClustering.RunningSolution.P_k_[CurrentMax] = Program.SpongePWeight / ParallelClustering.RunningSolution.Ncent_Global;
            double wgt = (1.0 - ParallelClustering.RunningSolution.P_k_[CurrentMax]);
            for (int ClusterIndex = 0; ClusterIndex < ParallelClustering.RunningSolution.Ncent_ThisNode; ClusterIndex++)
            {
                if (ParallelClustering.RunningSolution.LocalStatus[ClusterIndex] < 0)
                    continue;
                ParallelClustering.RunningSolution.P_k_[ClusterIndex] = wgt * ParallelClustering.RunningSolution.P_k_[ClusterIndex];
            }

            ParallelClustering.RunningSolution.DiffMalpha_k_Set = -1;
            ParallelClustering.RunningSolution.Ncent_ThisNode = CurrentMax + 1;
            ++ParallelClustering.RunningSolution.Ncent_Global;
            DistributedClusteringSolution.ManageMajorSynchronization(true);
            return true;

        }   // End AddSpongeCluster()

        public static void CalculateClusterStatus()
        {
            if(Program.DoLCMS)
                LCMSAnalyze.LCMSCalculateClusterStatus();
        }

        // Unlike Earlier versions, cluster labels are 0-based
        public static void OutputClusteringResults(string FileLabel)
        {
            bool DoFileOutput = true;
            if (FileLabel.Length == 0)
                DoFileOutput = false;

            Program.ClusterNumberOutput = ParallelClustering.RunningSolution.Ncent_Global;

            //  Set Global Counts and Positions
            double[][] GlobalCenters = new double[ParallelClustering.RunningSolution.Ncent_Global][];
            for (int GlobalClusterIndex = 0; GlobalClusterIndex < ParallelClustering.RunningSolution.Ncent_Global; GlobalClusterIndex++)
                GlobalCenters[GlobalClusterIndex] = new double[Program.ParameterVectorDimension];
            ClusteringSolution.SetGlobalClusterNumbers();

            // Generate Cluster Labels (cluster number) from current solution
            int[] labels = new int[DAVectorUtility.PointCount_Process];
            string[] Extralabels = new string[DAVectorUtility.PointCount_Process];

            //  Parallel Section setting cluster labels
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;

                for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                {
                    double distmax = -1.0;
                    int NearestClusterCreatedIndex = 0;
                    int IndirectSize = ParallelClustering.RunningSolution.NumClusters_alpha_[alpha];
                    for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                    {
                        if (ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex] > distmax)
                        {
                            distmax = ParallelClustering.RunningSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            NearestClusterCreatedIndex = ParallelClustering.RunningSolution.Map_alpha_PointertoCreatedIndex[alpha][IndirectClusterIndex];
                        }
                    }
                    int position = ClusteringSolution.UniversalMapping[NearestClusterCreatedIndex].GlobalClusterNumber;
                    if (position < 0)
                    {
                        Exception e = DAVectorUtility.SALSAError(" CreatedIndex has no Global Position " + NearestClusterCreatedIndex.ToString());
                        throw(e);
                    }
                    labels[alpha] = position;
                    Extralabels[alpha] = "";
                    if (Program.CompareSolution > 0)
                    {   // Add older clustering labels
                        int MedeaIndex = Program.MedeaClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int MclustIndex = Program.MclustClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int GoldenIndex = Program.GoldenPeaks.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        int OurIndex = Program.OurClusters.PointstoClusterIDs[alpha + DAVectorUtility.PointStart_Process];
                        Extralabels[alpha] = OurIndex.ToString() + " " + MedeaIndex.ToString() + " " + MclustIndex.ToString() + " " + GoldenIndex.ToString();
                    }
                }
            });  // End Parallel Section setting cluster label


            string directory = Path.GetDirectoryName(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string file = Path.GetFileNameWithoutExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile) + FileLabel + "-M" + Program.maxNcentperNode.ToString()
                + "-C" + ParallelClustering.RunningSolution.Ncent_Global.ToString() + Path.GetExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string ClusternumberFileName = Path.Combine(directory, file);

            int MPItag = 100;
            if (DAVectorUtility.MPI_Rank == 0)
            {
                if(DoFileOutput)
                    VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, ref labels, ref Program.PointPosition, ref Extralabels, DAVectorUtility.PointCount_Process, DAVectorUtility.PointStart_Process, false);

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                    int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;

                    for (int index = beginpoint; index < indexlen + beginpoint; index++)
                    {
                        Program.ClusterAssignments[index + DAVectorUtility.PointStart_Process] = labels[index];
                    }
                });  // End Parallel Section setting cluster assignments in process 0

                for (int MPISource = 1; MPISource < DAVectorUtility.MPI_Size; MPISource++)
                {
                    MPIPacket<int> fromsource;
                    fromsource = DAVectorUtility.MPI_communicator.Receive<MPIPacket<int>>(MPISource, MPItag);
                    int AwayArraySize = fromsource.NumberofPoints;
                    for (int index = 0; index < AwayArraySize; index++)
                        Program.ClusterAssignments[index + fromsource.FirstPoint] = fromsource.Marray[index];

                    double[][] AwayPointPositions = new double[AwayArraySize][];
                    for (int index = 0; index < AwayArraySize; index++)
                        AwayPointPositions[index] = new double[Program.ParameterVectorDimension];
                    MPIPacket<double> fromsourcedouble;
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        fromsourcedouble = DAVectorUtility.MPI_communicator.Receive<MPIPacket<double>>(MPISource, MPItag);
                        for (int LocalPointIndex = 0; LocalPointIndex < AwayArraySize; LocalPointIndex++)
                            AwayPointPositions[LocalPointIndex][VectorIndex] = fromsourcedouble.Marray[LocalPointIndex];
                    }
                    if(DoFileOutput)
                        VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, ref fromsource.Marray, ref AwayPointPositions, ref Extralabels, AwayArraySize, fromsource.FirstPoint, true);
                }

                int CenterLabelSize = ParallelClustering.RunningSolution.Ncent_Global;
                if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                    --CenterLabelSize;
                int[] CenterLabel = new int[CenterLabelSize];
                double[][] CenterPositions = new double[CenterLabelSize][];
                string[] ExtraCenterlabels = new string[CenterLabelSize];

                int decrement = 0;
                for (int StrippedClusterIndex = 0; StrippedClusterIndex < CenterLabelSize; StrippedClusterIndex++)  
                {
                    CenterPositions[StrippedClusterIndex] = new double[Program.ParameterVectorDimension];
                    if (ClusteringSolution.TotalClusterSummary.SpongeCluster == StrippedClusterIndex)
                        decrement = 1;
                    int ClusterIndex = StrippedClusterIndex + decrement;
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        CenterPositions[StrippedClusterIndex][VectorIndex] = ClusteringSolution.TotalClusterSummary.CenterPosition[ClusterIndex][VectorIndex];
                    CenterLabel[StrippedClusterIndex] = Math.Max(100000000, CenterLabelSize);
                    int OtherIndices = -1;
                    ExtraCenterlabels[StrippedClusterIndex] = "";
                    if (Program.CompareSolution > 0)
                    {
                        int OurIndex = ClusterIndex;
                        ExtraCenterlabels[StrippedClusterIndex] = OurIndex.ToString() + " " + OtherIndices.ToString() + " " + OtherIndices.ToString() + " " + OtherIndices.ToString();
                    }

                }

                // In first column for centers, we output cluster number + Total Point Count
                if(DoFileOutput)
                    VectorAnnealIterate.WriteClusterFile(ClusternumberFileName, ref CenterLabel, ref CenterPositions, ref ExtraCenterlabels, CenterLabelSize, DAVectorUtility.PointCount_Global, true);

            }   // End Root Process that receives cluster assignments from afar
            else
            {
                MPIPacket<int> tosend = new MPIPacket<int>(DAVectorUtility.PointCount_Process);
                tosend.FirstPoint = DAVectorUtility.PointStart_Process;
                tosend.NumberofPoints = DAVectorUtility.PointCount_Process;
                labels.CopyTo(tosend.Marray, 0);
                DAVectorUtility.MPI_communicator.Send<MPIPacket<int>>(tosend, 0, MPItag);
                MPIPacket<double> tosenddouble = new MPIPacket<double>(DAVectorUtility.PointCount_Process);
                for(int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    tosenddouble.FirstPoint = DAVectorUtility.PointStart_Process;
                    tosenddouble.NumberofPoints = DAVectorUtility.PointCount_Process;
                    for (int LocalPointIndex = 0; LocalPointIndex < DAVectorUtility.PointCount_Process; LocalPointIndex++)
                        tosenddouble.Marray[LocalPointIndex] = Program.PointPosition[LocalPointIndex][VectorIndex];
                    DAVectorUtility.MPI_communicator.Send<MPIPacket<double>>(tosenddouble, 0, MPItag);
                }
            }
            DAVectorUtility.MPI_communicator.Broadcast<int>(ref Program.ClusterAssignments, 0);
            DAVectorUtility.MPI_communicator.Barrier();

            DAVectorUtility.PreciseTimer.Start();   // Restart Timer
            return;

        }   // End OutputClusteringResults

        public static void Output3DClusterLabels(string extraname)
        {
            if (DAVectorUtility.MPI_Rank != 0)
                return;
            if (Program.RW3DData <= 0)
                return;

            string directory = Path.GetDirectoryName(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string file = Path.GetFileNameWithoutExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile) + "-Plot3D" + "-M" + Program.maxNcentperNode.ToString()
                + "-C" + ParallelClustering.RunningSolution.Ncent_Global.ToString() + Path.GetExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string ClusternumberFileName = Path.Combine(directory, extraname + "-" + file);

            Write3DClusterFile(ClusternumberFileName, ref Program.ClusterAssignments, 0, false);

        }   // End Output3DClusterLabels()

        public static void WriteClusterFile(string fname, ref int[] labels, ref double [][] PointPositions, ref string[] ExtraLabels, int dataPoints, int startposition, bool append)
        {
            FileMode mode = append ? FileMode.Append : FileMode.Create;

            using (FileStream stream = new FileStream(fname, mode, FileAccess.Write, FileShare.None))
            {
                using (StreamWriter writer = new StreamWriter(stream))
                {
                    for (int i = 0; i < dataPoints; i++)
                    {
                        string line = (i + startposition).ToString();
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            line += " " + PointPositions[i][VectorIndex].ToString("F6");
                        if (Program.ParameterVectorDimension == 2)
                            line += " 0.0";
                        line += " " + labels[i].ToString() + " " + ExtraLabels[i];
                        writer.WriteLine(line);
                    }
                }
            }
            if (Program.SigmaMethod > 0)
            {
                string fname1 = fname.Replace(".txt", "Scaled_.txt");
                using (FileStream stream1 = new FileStream(fname1, mode, FileAccess.Write, FileShare.None))
                {
                    using (StreamWriter writer1 = new StreamWriter(stream1))
                    {
                        double tmp;
                        for (int i = 0; i < dataPoints; i++)
                        {
                            string line = (i + startposition).ToString();
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            {
                                tmp = PointPositions[i][VectorIndex];
                                if ((Program.SigmaMethod >= 2) && (VectorIndex == 0) && (tmp != 0.0) )
                                    tmp = Math.Log(tmp);
                                tmp = tmp / Program.SigmaVectorParameters_i_[VectorIndex];
                                line += " " + tmp.ToString("F6");
                            }
                            if (Program.ParameterVectorDimension == 2)
                                line += " 0.0";
                            line += " " + labels[i].ToString();
                            writer1.WriteLine(line);
                        }
                    }
                }
            }
        }   // End WriteClusterFile

        public static void Write3DClusterFile(string fname, ref int[] labels, int startposition, bool append)
        {
            FileMode mode = append ? FileMode.Append : FileMode.Create;

            using (FileStream stream = new FileStream(fname, mode, FileAccess.Write, FileShare.None))
            {
                using (StreamWriter writer = new StreamWriter(stream))
                {
                    for (int i = 0; i < DAVectorUtility.PointCount_Global; i++)
                    {
                        string line = (i + startposition).ToString();
                        for (int VectorIndex = 0; VectorIndex < Program.RW3DData; VectorIndex++)
                            line += " " + Program.FullPoint3DPosition[i][VectorIndex].ToString("F6");
                        if (Program.RW3DData == 2)
                            line += " 0.0";
                        line += " " + labels[i].ToString();
                        writer.WriteLine(line);
                    }
                }
            }

        }   // End Write3DClusterFile

        // Unlike Earlier versions, cluster labels are 0-based
        public static void SimpleOutputClusteringResults(string FileLabel, double[][] ClusterCenters)
        {   // Output using precalculated Cluster Assignment Array
            string directory = Path.GetDirectoryName(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string file = Path.GetFileNameWithoutExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile) + FileLabel + "-M" + Program.maxNcentperNode.ToString()
                + "-C" + ParallelClustering.RunningSolution.Ncent_Global.ToString() + Path.GetExtension(Program._configurationManager.DAVectorSpongeSection.ClusterFile);
            string ClusternumberFileName = Path.Combine(directory, file);

            int MPItag = 100;
            if (DAVectorUtility.MPI_Rank == 0)
            {
                VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, ref Program.ClusterAssignments, ref Program.PointPosition, DAVectorUtility.PointCount_Process,
                    DAVectorUtility.PointStart_Process, DAVectorUtility.PointStart_Process, false);

                for (int MPISource = 1; MPISource < DAVectorUtility.MPI_Size; MPISource++)
                {
                    MPIPacket<int> fromsource;
                    fromsource = DAVectorUtility.MPI_communicator.Receive<MPIPacket<int>>(MPISource, MPItag);
                    int AwayArraySize = fromsource.NumberofPoints;

                    double[][] AwayPointPositions = new double[AwayArraySize][];
                    for (int index = 0; index < AwayArraySize; index++)
                        AwayPointPositions[index] = new double[Program.ParameterVectorDimension];
                    MPIPacket<double> fromsourcedouble;
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        fromsourcedouble = DAVectorUtility.MPI_communicator.Receive<MPIPacket<double>>(MPISource, MPItag);
                        for (int LocalPointIndex = 0; LocalPointIndex < AwayArraySize; LocalPointIndex++)
                            AwayPointPositions[LocalPointIndex][VectorIndex] = fromsourcedouble.Marray[LocalPointIndex];
                    }
                    VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, ref Program.ClusterAssignments, ref AwayPointPositions, AwayArraySize, fromsource.FirstPoint, fromsource.FirstPoint, true);
                }

                int CenterLabelSize = ParallelClustering.RunningSolution.Ncent_Global;
                if (ParallelClustering.RunningSolution.SpongeCluster >= 0)
                    --CenterLabelSize;
                int[] CenterLabel = new int[CenterLabelSize];

                int decrement = 0;
                for (int StrippedClusterIndex = 0; StrippedClusterIndex < CenterLabelSize; StrippedClusterIndex++)
                {
                    if (ClusteringSolution.TotalClusterSummary.SpongeCluster == StrippedClusterIndex)
                        decrement = 1;
                    int ClusterIndex = StrippedClusterIndex + decrement;
                    CenterLabel[StrippedClusterIndex] = Math.Max(100000000, CenterLabelSize);
                }

                // In first column for centers, we output cluster number + DAVectorUtility.PointCount_Global
                VectorAnnealIterate.SimpleWriteClusterFile(ClusternumberFileName, ref CenterLabel, ref ClusterCenters, CenterLabelSize, DAVectorUtility.PointCount_Global, 0, true);

            }   // End Root Process that receives cluster assignments from afar
            else
            {
                MPIPacket<int> tosend = new MPIPacket<int>(DAVectorUtility.PointCount_Process);
                tosend.FirstPoint = DAVectorUtility.PointStart_Process;
                tosend.NumberofPoints = DAVectorUtility.PointCount_Process;
                for (int index = 0; index < tosend.NumberofPoints; index++)
                    tosend.Marray[index] = Program.ClusterAssignments[index + tosend.FirstPoint];
                DAVectorUtility.MPI_communicator.Send<MPIPacket<int>>(tosend, 0, MPItag);
                MPIPacket<double> tosenddouble = new MPIPacket<double>(DAVectorUtility.PointCount_Process);
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    tosenddouble.FirstPoint = DAVectorUtility.PointStart_Process;
                    tosenddouble.NumberofPoints = DAVectorUtility.PointCount_Process;
                    for (int LocalPointIndex = 0; LocalPointIndex < DAVectorUtility.PointCount_Process; LocalPointIndex++)
                        tosenddouble.Marray[LocalPointIndex] = Program.PointPosition[LocalPointIndex][VectorIndex];
                    DAVectorUtility.MPI_communicator.Send<MPIPacket<double>>(tosenddouble, 0, MPItag);
                }
            }

            DAVectorUtility.PreciseTimer.Start();   // Restart Timer
            return;

        }   // End SimpleOutputClusterLabels(filename)


        public static void SimpleWriteClusterFile(string fname, ref int[] labels, ref double[][] PointPositions, int dataPoints, int startposition1, int startposition2, bool append)
        {
            FileMode mode = append ? FileMode.Append : FileMode.Create;

            using (FileStream stream = new FileStream(fname, mode, FileAccess.Write, FileShare.None))
            {
                using (StreamWriter writer = new StreamWriter(stream))
                {
                    for (int i = 0; i < dataPoints; i++)
                    {
                        string line = (i + startposition1).ToString();
                        for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            line += " " + PointPositions[i][VectorIndex].ToString("F6");
                        if (Program.ParameterVectorDimension == 2)
                            line += " 0.0";
                        line += " " + labels[i + startposition2].ToString();
                        writer.WriteLine(line);
                    }
                }
            }
        }   // End SimpleWriteClusterFile

    }	// End class dist

}   // End namespace cluster
