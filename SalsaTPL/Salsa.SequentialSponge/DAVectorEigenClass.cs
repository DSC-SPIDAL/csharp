using System;
using System.Threading;
using System.Threading.Tasks;
using MPI;

using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    //	Class vectorclass **************************************************

    class vectorclass
    {

        public double Eigenvalue; // Current eigenvalue
        public int EigenStatus; // Indicator of status of eigenvalue 
        public double[] Eigenvector;    // Eigenvector

        public double[][] CenterEigenvector;
        public double[] CenterEigenvalue;
        public double[] InitVector;
        public double[] FirstTerm;
        public int[] CenterEigenstatus;
        public int[] CenterEigenconvergence;

        public ClusteringSolution CurrentSolution;

        // Find the minimum eigenvalue of second derivative matrix -- called from shouldSplit
        public void getEigenvaluefromMatrix(double[,] SecondDerivMatrix)
        {
            
            Eigenvector = new double[Program.ParameterVectorDimension];
            if (!Program.CalculateEigenvaluesfromMatrix)
            {   // Re-use Earlier Calculation of results for all clusters
            }

            //  Calculate Eigenvalues from Matrix
            if (Program.ParameterVectorDimension != 2)
            {
                Exception e = DAVectorUtility.SALSAError(" Illegal Vector Dimension " + Program.ParameterVectorDimension.ToString());
                throw(e);
            }

            //  Case of Two Dimensions
            EigenStatus = 1;
            double tmp = SecondDerivMatrix[0, 0] - SecondDerivMatrix[1, 1];
            tmp = tmp * tmp + 4.0 * SecondDerivMatrix[0, 1] * SecondDerivMatrix[0, 1];
            Eigenvalue = 0.5 * (SecondDerivMatrix[0, 0] + SecondDerivMatrix[1, 1] - Math.Sqrt(tmp) );
            Eigenvector[0] = -SecondDerivMatrix[0, 1];
            Eigenvector[1] = SecondDerivMatrix[1, 1] - Eigenvalue;

            // Normalize
            tmp = 0.0;
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                tmp += Eigenvector[VectorIndex] * Eigenvector[VectorIndex];
            tmp = 1.0 / Math.Sqrt(tmp);
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                Eigenvector[VectorIndex] *= tmp;

        }   // End getEigenvalue(double[,] SecondDerivMatrix)

        public void SetAllEigenvaluesIteratively(ClusteringSolution Solution)
        {
            if (Solution.DistributedExecutionMode)
            {
                Exception e = DAVectorUtility.SALSAError(" Illegal Eigenvalue and Parallelization Combination ");
                throw (e);
            }
            if(Program.SigmaMethod > 0 )
            {
                Exception e = DAVectorUtility.SALSAError(" Illegal Eigenvalue and Sigma Method Combination " + Program.SigmaMethod.ToString());
                throw (e);
            }
            this.CurrentSolution = Solution;
            this.CenterEigenvector = this.CurrentSolution.Eigenvector_k_i;
            this.CenterEigenvalue = this.CurrentSolution.Eigenvalue_k;
            this.InitVector = new double[Program.ParameterVectorDimension];
            this.FirstTerm = new double[this.CurrentSolution.Ncent_Global];
            this.CenterEigenstatus = new int[this.CurrentSolution.Ncent_Global];
            this.CenterEigenconvergence = new int[this.CurrentSolution.Ncent_Global];

            Random random = new Random();
            double InitNorm = 0.0; 
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
            {
                InitVector[VectorIndex] = -0.5 + random.NextDouble();
                InitNorm += InitVector[VectorIndex] * InitVector[VectorIndex];
            }
            InitNorm = 1.0/Math.Sqrt(InitNorm);
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                InitVector[VectorIndex] *= InitNorm;

            //  Initialization Loop over Clusters
            int somethingtodo = 0;
            for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
            {
                this.CenterEigenconvergence[ClusterIndex] = 0;
                this.CenterEigenstatus[ClusterIndex] = 0;
                this.FirstTerm[ClusterIndex] = 0;
                if(this.CurrentSolution.Splittable_k_[ClusterIndex] != 1)
                    continue;
                ++somethingtodo;
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    this.CenterEigenvector[ClusterIndex][VectorIndex] = InitVector[VectorIndex];
            }   // End Loop over Clusters
            if( somethingtodo == 0 )
                return;

            GlobalReductions.FindVectorDoubleSum FindClusterFirstTerm = new GlobalReductions.FindVectorDoubleSum(DAVectorUtility.ThreadCount, this.CurrentSolution.Ncent_Global);
            GlobalReductions.FindDoubleSum FindNumberScalarProducts = new GlobalReductions.FindDoubleSum(DAVectorUtility.ThreadCount);

            for (int NumPowerIterations = 0; NumPowerIterations < Program.PowerIterationLimit; NumPowerIterations++)
            {
                somethingtodo = 0;
                for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
                {
                    if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
                        continue;
                    if( this.CurrentSolution.Splittable_k_[ClusterIndex] != 1)
                        continue;
                    if( this.CenterEigenconvergence[ClusterIndex] == 0 )
                        ++somethingtodo;
                }
                if(somethingtodo == 0 )
                    break;

                GlobalReductions.FindVectorDoubleSum3 FindNewPowerVectors = new GlobalReductions.FindVectorDoubleSum3(DAVectorUtility.ThreadCount,
                    Program.ParameterVectorDimension, this.CurrentSolution.Ncent_Global);

                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    FindNewPowerVectors.startthread(ThreadNo);
                    double[] PartVector = new double[Program.ParameterVectorDimension];

                    int indexlen = DAVectorUtility.PointsperThread[ThreadNo];
                    int beginpoint = DAVectorUtility.StartPointperThread[ThreadNo] - DAVectorUtility.PointStart_Process;
                    for (int alpha = beginpoint; alpha < indexlen + beginpoint; alpha++)
                    {
                        int IndirectSize = this.CurrentSolution.NumClusters_alpha_[alpha];
                        
                        for (int IndirectClusterIndex = 0; IndirectClusterIndex < IndirectSize; IndirectClusterIndex++)
                        {   // Loop over Clusters for this point
                            int RealClusterIndex = -1;
                            int RemoteIndex = -1;
                            int ActiveClusterIndex = -1;
                            VectorAnnealIterate.ClusterPointersforaPoint(alpha, IndirectClusterIndex, ref RealClusterIndex, ref ActiveClusterIndex, ref RemoteIndex);
                            if (this.CurrentSolution.Splittable_k_[RealClusterIndex] != 1)
                                continue;

                            double Mvalue = this.CurrentSolution.M_alpha_kpointer_[alpha][IndirectClusterIndex];
                            if (NumPowerIterations == 0)
                                FindClusterFirstTerm.addapoint(ThreadNo, Mvalue, RealClusterIndex);
                            double multiplier = 0.0;
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                            {
                                PartVector[VectorIndex] = this.CurrentSolution.Y_k_i_[RealClusterIndex][VectorIndex] - Program.PointPosition[alpha][VectorIndex];
                                multiplier += PartVector[VectorIndex] * CenterEigenvector[RealClusterIndex][VectorIndex];
                            }
                            FindNumberScalarProducts.addapoint(ThreadNo, 1.0);

                            double wgt = Mvalue * multiplier;
                            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                                PartVector[VectorIndex] *= wgt;
                            FindNewPowerVectors.addapoint(ThreadNo, PartVector, RealClusterIndex);
                        }
                    }   // End Loop over points
                }); // End loop initialing Point dependent quantities

                FindNewPowerVectors.sumoverthreadsandmpi();
                for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
                {
                    if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
                        continue;
                    if( (this.CurrentSolution.Splittable_k_[ClusterIndex] != 1) || (this.CenterEigenconvergence[ClusterIndex] != 0) )
                        continue;
                    double[] sums = new double[3];  // Old.New Old.Old New.New
                    for(int loop = 0; loop < 3; loop++)
                        sums[loop] = 0.0;
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        int TotalIndex = VectorIndex + ClusterIndex * Program.ParameterVectorDimension;
                        double newvalue = FindNewPowerVectors.TotalVectorSum[TotalIndex];
                        double oldvalue = CenterEigenvector[ClusterIndex][VectorIndex];
                        sums[0] += oldvalue * newvalue;
                        sums[1] += oldvalue * oldvalue;
                        sums[2] += newvalue * newvalue;
                        CenterEigenvector[ClusterIndex][VectorIndex] = newvalue;
                    }

                    //  Decide if finished and set eigenvalue
                    double CandidateEigenvalue = sums[0] / sums[1];
                    bool LegalEigenvalue = (CandidateEigenvalue > 0.0);
                    DAVectorUtility.SynchronizeMPIvariable(ref LegalEigenvalue);

                    //	Check if converged
                    //	Do this in one process ONLY 
                    if ((NumPowerIterations > 5) && LegalEigenvalue)
                    { // Arbitrary choice for Number of Power Iterations Cut

                        int EigenvalueDone = 0;
                        if (DAVectorUtility.MPI_Rank == 0)
                        { // Decisions can only be made in one process			
                            if (Math.Abs(CandidateEigenvalue - this.CenterEigenvalue[ClusterIndex]) > CandidateEigenvalue * Program.eigenvaluechange) ++EigenvalueDone;
                            double delta = sums[2] - 2.0 * sums[0] * CandidateEigenvalue + sums[1] * CandidateEigenvalue * CandidateEigenvalue;   // (Ax- Eigenvalue*Axold)**2
                            if (Math.Abs(delta) > CandidateEigenvalue * CandidateEigenvalue * Program.eigenvectorchange) ++EigenvalueDone;
                        }   // End Test on Convergence
                        DAVectorUtility.SynchronizeMPIvariable(ref EigenvalueDone);

                        if (EigenvalueDone == 0)
                            this.CenterEigenconvergence[ClusterIndex] = 1 + NumPowerIterations;
                    }
                    this.CenterEigenvalue[ClusterIndex] = CandidateEigenvalue;

                    //  Normalize current Power Vector to 1
                    double wgt = 1.0 / Math.Sqrt(sums[2]);
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        CenterEigenvector[ClusterIndex][VectorIndex] *= wgt;
                }   // End Loop over Clusters

            }   //  End Loop over NumPowerIterations

            FindClusterFirstTerm.sumoverthreadsandmpi();
            FindNumberScalarProducts.sumoverthreadsandmpi();
            Program.SumEigenSPCalcs += FindNumberScalarProducts.Total;

            for (int ClusterIndex = 0; ClusterIndex < this.CurrentSolution.Ncent_Global; ClusterIndex++)
            {
                this.CenterEigenstatus[ClusterIndex] = 0;
                if (this.CurrentSolution.LocalStatus[ClusterIndex] != 1)
                    continue;
                if ((this.CurrentSolution.Splittable_k_[ClusterIndex] != 1) || (this.CenterEigenconvergence[ClusterIndex] <= 0))
                    continue;
                this.CenterEigenstatus[ClusterIndex] = 1;
                this.FirstTerm[ClusterIndex] = FindClusterFirstTerm.TotalVectorSum[ClusterIndex];
                double tmp = this.CenterEigenvalue[ClusterIndex] / this.CurrentSolution.Temperature;
                this.CenterEigenvalue[ClusterIndex] = this.FirstTerm[ClusterIndex] - tmp;
            }


        }   // End SetEigenvaluesIteratively(ClusteringSolution Solution)

        public void getEigenvaluefromIteration(int ClusterIndex)
        {   // Eigenvector stored in Solution

            if(this.CenterEigenconvergence[ClusterIndex] <= 0 )
                this.CenterEigenstatus[ClusterIndex] = 0;
            this.EigenStatus = this.CenterEigenstatus[ClusterIndex];
            if (this.EigenStatus <= 0)
                this.Eigenvalue = 1.0;
            else
                this.Eigenvalue = this.CenterEigenvalue[ClusterIndex];

        }   // End getEigenvaluefromIteration(int ClusterIndex)

    }	// End vectorclass

}   // End namespace Salsa.DAVectorSponge
