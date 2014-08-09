using System;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using MPI;
using Salsa.Core;
using SALSALibrary;

namespace Salsa.DAVectorSponge
{
    class LCMSAnalyze
    {
        public static void InitializeLCMS()
        {
            Program.CalculateCorrelationMatrix = true;
            Program.CalculateIndividualWidths = true;
            Program.Printeigenvectors = true;
            Program.CalculateEigenvaluesfromMatrix = true;

        }   // End InitializeLCMS()

        public static void JustAnalyze()
        {
            SpongePointAnalysis();
            
            //  Set up to read full file
            Program.SelectedInputLabel = -100000000;
            Program.InputFileType = 1;
            Program.RestartTemperature = 0.025;

            Program.CompareSolution = 1;
            Program.ComparisonSelectedInputLabel = 2;
            Program.ComparisonInputFileType = 0;
            Program.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

            Program.Refinement = false;
            Program.UseSponge = true;

            //  Set Final Sigma[0]
            Program.SigmaVectorParameters_i_[0] = 0.00000598;
            Program.SigmaMethod = 2;


            //  Please set these
            Program.SelectedInputLabel = -1000000;    // Old Data Style
            Program.RestartSelectedInputLabel = -1000000; // Old Data Style

            // Modern Style
            Program.SelectedInputLabel = -100000000;    // New Data Style
            Program.RestartSelectedInputLabel = -100000000; // New Data Style

            //  Set Sponge
            Program.SpongeFactor1 = 2.0;
            Program.SpongeFactor2 = 2.0;
            Program.SpongeFactor1 = 1.0;
            Program.SpongeFactor2 = 1.0; 

            //  No Sponge
           /*
            Program.SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
            Program.SpongeTemperature2 = -1.0;
            Program.NearbySpongePointLimit = 4.0;
            Program.UseSponge = false;
            */

        }

        public static void JoinEquilibriate()
        {   // Join 2 Files Together (typically original and analysis of its sponge) and equilibriate
            //  Set Distance Matrix File to be output (ClusterFinal) of Previous job i.e. basic file
            //  Set Label File to be second (sponge analysis) or update

            SpongePointAnalysis();

            //  Set up to read full files
            Program.SelectedInputLabel = -100000000;
            Program.InputFileType = 1;
            Program.RestartSelectedInputLabel = -100000000;
            Program.RestartInputFileType = 1;
            Program.RestartClusterFile = Program._configurationManager.DAVectorSpongeSection.DistanceMatrixFile;

            Program.Refinement = true;
            Program.RestartTemperature = 0.025; // was 0.2
            Program.Tminimum = 0.005;


            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);
            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);

            //  Set Final Sigma[0] for low temperature restart
            Program.SigmaVectorParameters_i_[0] = 0.00000598;
            Program.SigmaMethod = 2;

            Program.CompareSolution = 1;
            Program.ComparisonSelectedInputLabel = 2;
            Program.ComparisonInputFileType = 0;
            Program.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";

            //  Please set these
            Program.UseSponge = true;
            Program.SpongeFactor1 = 2.0;
            Program.SpongeFactor2 = 2.0;

            // Old Sponge 1
            Program.SpongeFactor1 = 1.0;
            Program.SpongeFactor2 = 1.0;
            Program.RestartTemperature = 1.0; // was 0.025
            Program.SelectedInputLabel = -1000000;
            Program.RestartSelectedInputLabel = -1000000;

            // Redo DA2D

            //  No Sponge
            Program.SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
            Program.SpongeTemperature2 = -1.0;
            Program.NearbySpongePointLimit = 4.0;
            Program.UseSponge = false;
            Program.RestartTemperature = 0.4;
            Program.SelectedInputLabel = -100000000;
            Program.RestartSelectedInputLabel = -100000000;

        }   // End JoinEquilibriate()

        public static void ChangeSpongeFactor()
        {

            JoinEquilibriate();

            //  Reduce Sponge
            Program.SpongeFactor1 = 2.0;
            Program.SpongeFactor2 = 1.0;
            Program.SpongeTemperature1 = 5.0;
            Program.SpongeTemperature2 = 2.0;
            Program.RestartTemperature = 5.0;

            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);

        }   // End ChangeSpongeFactor()

        public static void SpongePointAnalysis()
        {   // Analysis of Sponge Points from a Previous run
            //  Set Distance Matrix File to be output (ClusterFinal) of Previous job

            FreshFullAnalysis();
            Program.Refinement = true;
            Program.MaxNumberSplitClusters = 32;

            // Cooling
            Program.InitialCoolingFactor1 = 0.9875;
            Program.FineCoolingFactor1 = 0.9975;
            Program.InitialCoolingFactor2 = 0.99375;
            Program.FineCoolingFactor2 = 0.999375;

            Program.InputFileType = 1;  // Output File of Previous Job
            Program.SelectedInputLabel = 0; // Select Sponge Points
            Program.RestartTemperature = -1.0;
            Program.CompareSolution = -1;

            //  Analyze Sponge 1
            Program.SpongeFactor2 = 1.0;

        }   // End SpongePointAnalysis()

        public static void FreshFullAnalysis()
        {   // Complete Fresh Analysis

            //  Center Parameters
            Program.maxNcentperNode = 50000;
            Program.maxNcentTOTAL = 100000;
            Program.targetNcentperPoint = 8;
            Program.maxNcentperPoint = 17;
            Program.targetMinimumNcentperPoint = 2;
            Program.maxNcentCreated = 5 * Program.maxNcentperNode;

            //  Ending Parameters
            Program.Iterationatend = 2000;
            Program.Malpha_MaxChange = 0.005;
            Program.Tminimum = 0.025;

            //  Split Parameters
            Program.MaxNumberSplitClusters = 256;
            Program.ClusterLimitforDistribution = 256;

            Program.Waititerations = 1;
            Program.Waititerations_Converge = 4;

            //  Set Print Options
            Program.PrintInterval = 50;
            Program.ClusterPrintNumber = 5;
            DAVectorUtility.DebugPrintOption = 1;

            //  Magic Temperatures
            Program.MagicTemperatures[0] = -1.0;
            Program.MagicTemperatures[1] = -1.0;
            Program.MagicTemperatures[2] = -1.0;
            Program.MagicTemperatures[3] = -1.0;
            Program.MagicTemperatures[4] = -1.0;

            // Cooling
            Program.InitialCoolingFactor1 = 0.9875;
            Program.FineCoolingFactor1 = 0.9975;
            Program.InitialCoolingFactor2 = 0.99375;
            Program.FineCoolingFactor2 = 0.999375;

            //  Slow Cooling
            Program.CoolingTemperatureSwitch = 30.0;
            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);
            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);

            //  Exponential Cuts
            Program.ExpArgumentCut1 = 20.0;
            Program.ExpArgumentCut2 = 40.0;

            // Correct center sizing parameters
            Program.MinimumScaledWidthsquaredtosplit = 0.5;
            Program.ScaledSquaredDistanceatClosenessTest = 0.35;
            Program.MinimumCountforCluster_C_kwithSponge = 1.5;
            Program.MinimumCountforCluster_C_k = 0.5;
            Program.MinimumCountforCluster_Points = -1;

            Program.TemperatureforClosenessTest = 4.0;
            Program.TemperatureforClosenessTest = 9.0;

            // Set Data Reading parameters including Charge
            //  DAVectorSpongeSection.DistanceMatrixFile is basic
            //  DAVectorSpongeSection.LabelFile is update
            Program.SelectedInputLabel = -60;
            Program.Replicate = 1;
            Program.ComparisonInputFileType = 0;
            Program.ComparisonClusterFile = "c:\\remote\\input\\DarTBFULLSorted.txt";
            Program.RestartSelectedInputLabel = -100000000;
            Program.RestartInputFileType = 1;
            Program.RestartClusterFile = Program._configurationManager.DAVectorSpongeSection.DistanceMatrixFile;
            Program.ComparisonSelectedInputLabel = 2;
            Program.SelectedInputLabel = 2; // Charge 2
            Program.InputFileType = 0; //  Original Harvard format
            Program.RestartTemperature = -1.0;
            Program.CompareSolution = -1;

            //  Sigma Annealing
            Program.SigmaVectorParameters_i_[0] = 0.000598;
            Program.SigmaVectorParameters_i_[1] = 2.35;
            //       Program.SigmaMethod = 2;
            Program.SigmaMethod = 3;
            Program.FinalTargetSigma0 = 0.00000598;
            Program.FinalTargetTemperature = 12.0;
            Program.FinalTargetTemperature = 30.0;

            //  Sponge Parameters
            Program.UseSponge = false;
            Program.SpongePWeight = 0.1;
            Program.CreateSpongeScaledSquaredWidth = 10.0;

            Program.SpongeFactor1 = 12.0;
            Program.SpongeFactor1 = 45.0;
            Program.SpongeFactor2 = 2.0;
            Program.SpongeTemperature1 = 9.0;
            Program.SpongeTemperature1 = 30.0;
            Program.SpongeTemperature2 = 2.5;
            Program.SpongeTemperature2 = 2.0;

            // No Sponge
            Program.CreateSpongeScaledSquaredWidth = -1.0; // Create a Sponge Cluster if doesn't exist already when average  squared scaled width reaches this
            Program.SpongeTemperature1 = -1.0;     // Minimum Temperature where Sponge Introduced
            Program.SpongeTemperature2 = -1.0;

            // Run F7 
            Program.SigmaVectorParameters_i_[0] = 0.000598;
            Program.SigmaMethod = 3;
            Program.FinalTargetSigma0 = 0.00000598;
            Program.FinalTargetTemperature = 12.0;

            Program.ScaledSquaredDistanceatClosenessTest = 0.5;
            Program.TemperatureforClosenessTest = 4.0;
            Program.MaxNumberSplitClusters = 256;
            Program.MinimumScaledWidthsquaredtosplit = 1.5;
            Program.ClusterLimitforDistribution = 256;
            Program.Tminimum = 0.1;
            Program.RestartTemperature = -1.0;
            Program.Iterationatend = 400;

            Program.UseSponge = false;
            Program.SpongeFactor1 = 12.0;
            Program.SpongeFactor2 = 3.0;
            Program.SpongeTemperature1 = 9.0;
            Program.SpongeTemperature2 = 3.0;
            Program.SpongePWeight = 0.1;
            Program.CreateSpongeScaledSquaredWidth = 10.0;

            Program.CoolingTemperatureSwitch = 12.0;
            Program.InitialCoolingFactor1 = 0.9875000;
            Program.FineCoolingFactor1 = 0.9975000;
            Program.FineCoolingFactor2 = 0.9993750;
            Program.InitialCoolingFactor2 = 0.9937500;

            Program.Malpha_MaxChange = 0.005;
            Program.Malpha_MaxChange1 = 0.005;
            Program.MinimumCountforCluster_C_kwithSponge = 1.5;
            Program.MinimumCountforCluster_C_k = 0.5;
            Program.MinimumCountforCluster_Points = -1;

            //  F7 PLUS
            /*
            Program.InitialCoolingFactor1 = Math.Sqrt(Program.InitialCoolingFactor1);
            Program.FineCoolingFactor1 = Math.Sqrt(Program.FineCoolingFactor1);
            Program.FineCoolingFactor2 = Math.Sqrt(Program.FineCoolingFactor2);
            Program.InitialCoolingFactor2 = Math.Sqrt(Program.InitialCoolingFactor2);
            */
 
    }   // End FreshFullAnalysis()

        public static void LCMSCalculateClusterStatus()
        {
            double cut = Program.NearbySpongePointLimit;
            if (cut < 0.0)
                cut = Program.SpongeFactor;
            if (cut < 0.0)
                cut = 1.0;
            ClusteringSolution.SetGlobalClusterNumbers();
            VectorAnnealIterate.OutputClusteringResults("");   // Set Point -- Cluster links

            Program.ClusterStatus = new ClusterQuality(ClusteringSolution.TotalClusterSummary.NumberofCenters, ClusteringSolution.TotalClusterSummary.SpongeCluster, Program.NumberNearbyClusters, cut);
            
            Program.ClusterStatus.SetClusterStatistics();
            Program.ClusterStatus.SetPointStatistics();
            Program.ClusterStatus.SetNearbyClusters();
            Program.ClusterStatus.OutputStatus();

            //  Set up Cluster Comparisons
            //  Convert Sponge to singleton clusters
            //  Note FullClusterNumber will end up as 1 more than total of clusters if there is a sponge
            if (Program.CompareSolution <= 0)
                return;
            int FullClusterNumber = ClusteringSolution.TotalClusterSummary.NumberofCenters;
            int SpongeCluster = ClusteringSolution.TotalClusterSummary.SpongeCluster;
            for (int GlobalPointIndex = 0; GlobalPointIndex < DAVectorUtility.PointCount_Global; GlobalPointIndex++)
            {
                int NewClusterIndex = Program.ClusterAssignments[GlobalPointIndex];
                if (NewClusterIndex == SpongeCluster)
                {
                    NewClusterIndex = FullClusterNumber;
                    ++FullClusterNumber;
                }
                Program.OurClusters.PointstoClusterIDs[GlobalPointIndex] = NewClusterIndex;
            }
            Program.OurClusters.setup();
            LCMSAnalyze.ClusterComparison();

        }   // End LCMSCalculateClusterStatus()

        public static void ClusterComparison()
        {
            if (Program.CompareSolution <= 0)
                return;

            
            //  GoldenExamination of Cuts and Differences
            GoldenExamination.NumberClusteringMethods = 4;

            GoldenExamination.OneDHistogramInterval = 0.1;
            GoldenExamination.TwoDHistogramInterval = 0.1;
            GoldenExamination.OneDHistogramMaxVal = 8.0;
            GoldenExamination.TwoDHistogramMaxVal = 8.0;
            GoldenExamination.OneDHistogramSize = (int)Math.Floor((GoldenExamination.OneDHistogramMaxVal + 0.0001) / GoldenExamination.OneDHistogramInterval) + 1;
            GoldenExamination.TwoDHistogramSize = (int)Math.Floor((GoldenExamination.TwoDHistogramMaxVal + 0.0001) / GoldenExamination.TwoDHistogramInterval) + 1;

            GoldenExamination.OneD_Center_HistogramInterval = 0.025;
            GoldenExamination.TwoD_Center_HistogramInterval = 0.025;
            GoldenExamination.OneD_Center_HistogramMaxVal = 8.0;
            GoldenExamination.TwoD_Center_HistogramMaxVal = 8.0;
            GoldenExamination.OneD_Center_HistogramSize = (int)Math.Floor((GoldenExamination.OneD_Center_HistogramMaxVal + 0.0001) / GoldenExamination.OneD_Center_HistogramInterval) + 1;
            GoldenExamination.TwoD_Center_HistogramSize = (int)Math.Floor((GoldenExamination.TwoD_Center_HistogramMaxVal + 0.0001) / GoldenExamination.TwoD_Center_HistogramInterval) + 1;

            if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
            {
                //  Handleloop = 0 Cut on 2D distance
                //  Handleloop = 1 Cut on each 1D distance (Removed)
                for (int handleloop = 0; handleloop < 1; handleloop++)
                {
                    for (int minloop = 0; minloop < 3; minloop++)
                    {
                        GoldenExamination GoldenExaminationContainer = new GoldenExamination(); // This sets Accumulation Record
                        GoldenExaminationContainer.HandleGoldenCluster = handleloop;
                        GoldenExaminationContainer.MinimumforAveraging = 3;
                        if (minloop == 1)
                            GoldenExaminationContainer.MinimumforAveraging = 15;
                        if (minloop == 2)
                            GoldenExaminationContainer.MinimumforAveraging = 40;
                        GoldenExaminationContainer.CutValue = 0.7;
                        GoldenExaminationContainer.CutMinimumforAveraging = GoldenExaminationContainer.MinimumforAveraging;
                        if (GoldenExaminationContainer.CutMinimumforAveraging > 10)
                            GoldenExaminationContainer.CutMinimumforAveraging = GoldenExaminationContainer.MinimumforAveraging - 4;


                        ArbitraryClustering[] ThreeMethods = new ArbitraryClustering[GoldenExamination.NumberClusteringMethods];
                        ThreeMethods[0] = Program.OurClusters;
                        ThreeMethods[1] = Program.MedeaClusters;
                        ThreeMethods[2] = Program.MclustClusters;
                        ThreeMethods[3] = Program.GoldenPeaks;
                        GoldenExaminationContainer.GoldenComparisons(Program.GoldenPeaks, ThreeMethods);
                        GoldenExaminationContainer.PrintAccumulation();
                    }
                }
            }

            //  Point Position Histograms
            GoldenExamination.OneDHistogramInterval = 0.05;
            GoldenExamination.TwoDHistogramInterval = 0.025;
            GoldenExamination.OneDHistogramMaxVal = 8.0;
            GoldenExamination.TwoDHistogramMaxVal = 8.0;
            GoldenExamination.OneDHistogramSize = (int)Math.Floor((GoldenExamination.OneDHistogramMaxVal + 0.0001) / GoldenExamination.OneDHistogramInterval) + 1;
            GoldenExamination.TwoDHistogramSize = (int)Math.Floor((GoldenExamination.TwoDHistogramMaxVal + 0.0001) / GoldenExamination.TwoDHistogramInterval) + 1;

            int ClusterCountcut = 5;
            for (int minloop = 0; minloop < 4; minloop++)
            {
                if (minloop == 2)
                    ClusterCountcut = 20;
                if (minloop == 3)
                    ClusterCountcut = 50;
                Program.OurClusters.HistogramPeaks(ClusterCountcut);
                Program.MedeaClusters.HistogramPeaks(ClusterCountcut);
                Program.MclustClusters.HistogramPeaks(ClusterCountcut);
                Program.GoldenPeaks.HistogramPeaks(ClusterCountcut);
                ClusterCountcut += 5;
            }

            //  Comparison of Basic Clusters
            DAVectorUtility.SALSAPrint(0, "\n**************** Statistics of Clustering in each method and selection to Golden Clusters\nWith means versus occupation count and 1D/2D Point-Center Histograms");
            Program.OurClusters.Statistics();
            Program.MedeaClusters.Statistics();
            Program.MclustClusters.Statistics();
            if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
                Program.GoldenPeaks.Statistics();
            
            if (Program.GoldenPeaks.MaxIndependentClusterIndex > 0)
            {
                DAVectorUtility.SALSAPrint(0, "\n Comparisons with Golden Clusters");
                Program.OurClusters.Difference(Program.GoldenPeaks);
                Program.MedeaClusters.Difference(Program.GoldenPeaks);
                Program.MclustClusters.Difference(Program.GoldenPeaks);
                Program.GoldenPeaks.Difference(Program.GoldenPeaks);
            }

            DAVectorUtility.SALSAPrint(0, "\n Comparisons with DAVector Clusters");
            Program.MedeaClusters.Difference(Program.OurClusters);
            Program.MclustClusters.Difference(Program.OurClusters);
            Program.OurClusters.Difference(Program.MedeaClusters);
            Program.OurClusters.Difference(Program.MclustClusters);
            
        }

    }   // End LCMSAnalyze

    public class ArbitraryClustering
    {
        public int NumberofPoints = 0;
        public int MaxClusterID = 0;
        public int MaxIndependentClusterIndex = 0;
        public string ClusterString = "";

        public int[] PointstoClusterIDs;
        public int[] ClusterIDtoIndex;
        public int[] ClusterIDtoGoldenIndex;
        public int[] ClusterIndextoID;
        public int[] ClusterCountsbyIndex;
        public int[] ClusterGoldenCountsbyIndex;

        public ArbitraryClustering(int NumberofPointsINPUT, string ClusterStringINPUT)
        {
            ClusterString = ClusterStringINPUT;
            NumberofPoints = NumberofPointsINPUT;
            PointstoClusterIDs = new int[NumberofPoints];
        }

        //  Rather Messy Indexing
        //  For each "type" (DAVS Medea Mcluster Golden) there is an array PointstoClusterIDs that maps each Point to a cluster. 
        //  This array is set in public static void ReadLabelsFromFile(string fname) for "foreign types"
        //  This array is set in LCMSCalculateClusterStatus() for "DAVS clusters" and ALL Sponger points are converted to individual clusters with one point each and stored at end of list
        //  PointstoClusterIDs is -1 if no cluster
        //  The Clusters are listed in order they they first appear in file
        //  this.MaxIndependentClusterIndex is number of clusters and index runs from 0 to this.MaxIndependentClusterIndex-1
        //  this.ClusterIndextoID[index] = is RawClusterIndex for this collection. RawClusterIndex will appear in PointstoClusterIDs for each point
        //  this.ClusterCountsbyIndex[index] is number of points in this cluster
        //  this.ClusterGoldenCountsbyIndex[index] is number of Golden points in this cluster
        //
        public void setup()
        {
            this.MaxClusterID = 0;
            for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
                this.MaxClusterID = Math.Max(this.MaxClusterID, this.PointstoClusterIDs[GlobalPointIndex]);
            ++this.MaxClusterID;
            this.ClusterIDtoIndex = new int[MaxClusterID];
            this.ClusterIDtoGoldenIndex = new int[MaxClusterID];
            for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
            {
                this.ClusterIDtoIndex[RawClusterIndex] = 0;
                this.ClusterIDtoGoldenIndex[RawClusterIndex] = 0;
            }
            for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
            {
                int RawClusterIndex = this.PointstoClusterIDs[GlobalPointIndex];
                if (RawClusterIndex < 0)
                    continue;
                this.ClusterIDtoIndex[RawClusterIndex]++;
                if (Program.GoldenPeaks.PointstoClusterIDs[GlobalPointIndex] >= 0)
                    this.ClusterIDtoGoldenIndex[RawClusterIndex]++;
            }
            this.MaxIndependentClusterIndex = 0;
            for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
            {
                int Count = this.ClusterIDtoIndex[RawClusterIndex];
                if (Count == 0)
                {
                    this.ClusterIDtoIndex[RawClusterIndex] = -1;
                    continue;
                }
                ++this.MaxIndependentClusterIndex;
            }

            this.ClusterIndextoID = new int[this.MaxIndependentClusterIndex];
            this.ClusterCountsbyIndex = new int[this.MaxIndependentClusterIndex];
            this.ClusterGoldenCountsbyIndex = new int[this.MaxIndependentClusterIndex];

            int index = 0;
            for (int RawClusterIndex = 0; RawClusterIndex < this.MaxClusterID; RawClusterIndex++)
            {
                int Count = this.ClusterIDtoIndex[RawClusterIndex];
                if (Count < 0)
                    continue;
                this.ClusterIndextoID[index] = RawClusterIndex;
                this.ClusterCountsbyIndex[index] = Count;
                this.ClusterGoldenCountsbyIndex[index] = this.ClusterIDtoGoldenIndex[RawClusterIndex];
                this.ClusterIDtoIndex[RawClusterIndex] = index;
                ++index;
            }
        }   // End setup

        public void Statistics()
        {   // Individual statistics for clustering

            for (int histogramloops = 0; histogramloops < 2; histogramloops++)
            {
                int[] ListofCounts;
                string HistogramType;
                int CountLimit;
                int[] Limits;
                if (histogramloops == 0)
                {
                    ListofCounts = this.ClusterCountsbyIndex;
                    CountLimit = 1;
                    HistogramType = "Full";
                    Limits = new int[4] {5, 10, 20 , 30};
                }
                else
                {
                    ListofCounts = this.ClusterGoldenCountsbyIndex;
                    CountLimit = 0;
                    HistogramType = "Select Golden Peaks";
                    Limits = new int[3] { 5, 15,40 };
                }
                int NumIntervals = Limits.Length + 1;
                int MaxCount = 0;
                double PWHammySingle = 0.0;
                double PWHammyReal = 0.0;
                double[] Width_xReal = new double[NumIntervals];
                double[] Width_yReal = new double[NumIntervals];
                int[] NumPeaksinCategory = new int[NumIntervals];
                int[] NumClustersinCategory = new int[NumIntervals];
                for (int Category = 0; Category < NumIntervals; Category++)
                {
                    Width_xReal[Category] = 0.0;
                    Width_yReal[Category] = 0.0;
                    NumPeaksinCategory[Category] = 0;
                    NumClustersinCategory[Category] = 0;
                }
                int NSingle = 0;
                int NReal = 0;

                for (int ClusterIndex = 0; ClusterIndex < this.MaxIndependentClusterIndex; ClusterIndex++)
                    MaxCount = Math.Max(MaxCount, ListofCounts[ClusterIndex]);
                int[] HistCounts = new int[MaxCount + 1];

                int GreaterThanOne = 0;
                double Average = 0.0;
                for (int ClusterIndex = 0; ClusterIndex < this.MaxIndependentClusterIndex; ClusterIndex++)
                {
                    int Count = ListofCounts[ClusterIndex];
                    if (Count < CountLimit)
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal Cluster Count " + this.ClusterString + " Index " + ClusterIndex.ToString()
                            + " Original Label " + this.ClusterIndextoID[ClusterIndex].ToString());
                        throw e;
                    }
                    HistCounts[Count]++;
                    if (Count <= 0)
                        continue;

                    int Category = NumIntervals - 1;
                    for (int FindCategory = 0; FindCategory < (NumIntervals - 1); FindCategory++)
                    {
                        if( Count < Limits[FindCategory])
                        {
                            Category = FindCategory;
                            break;
                        }
                    }
                    
                    int RawCluster = this.ClusterIndextoID[ClusterIndex];
                    if (Count > 1)
                    {   // Real Clusters
                        GreaterThanOne++;
                        Average += Count;

                        double LocalWidth_x = 0.0;
                        double LocalWidth_y = 0.0;
                        SetWidths(RawCluster, out LocalWidth_x, out LocalWidth_y);
                        NReal += Count;
                        Width_xReal[Category] += LocalWidth_x;
                        Width_yReal[Category] += LocalWidth_y;
                        NumPeaksinCategory[Category] += Count;
                        NumClustersinCategory[Category]++;
                        PWHammyReal += 0.5 * (LocalWidth_x + LocalWidth_y);
                    }
                    else
                    {   // Singletons
                        ++NSingle;
                        PWHammySingle += 0.5 * Program.SpongeFactor * Program.SpongeFactor;
                    }

                }
                if (GreaterThanOne > 0)
                    Average = Average / GreaterThanOne;

                double PWHammy = PWHammyReal + PWHammySingle;

                double Width_xReal_Total = 0.0;
                double Width_yReal_Total = 0.0;
                for (int Category = 0; Category < NumIntervals; Category++)
                {
                    Width_xReal_Total += Width_xReal[Category];
                    Width_yReal_Total += Width_yReal[Category];
                    if (NumPeaksinCategory[Category] == 0)
                        continue;
                    Width_xReal[Category] = Width_xReal[Category] / NumPeaksinCategory[Category];
                    Width_yReal[Category] = Width_yReal[Category] / NumPeaksinCategory[Category];
                }
                if (NReal > 0)
                {
                    Width_xReal_Total = Width_xReal_Total / NReal;
                    Width_yReal_Total = Width_yReal_Total / NReal;
                }

                string nextline = "\nBasic Statistics of Cluster " + ClusterString + " " + HistogramType + " Max Count " + MaxCount.ToString() + " Avg Count above 1 " + Average.ToString("F3")
                    + " Number of Clusters " + this.MaxIndependentClusterIndex.ToString() + " Above Count of One " + GreaterThanOne.ToString() + "\n"
                    + " Total Hamiltonian " + PWHammy.ToString("E4") + " Singletons " + NSingle.ToString() + " Single Hamiltonian " + PWHammySingle.ToString("E4")
                    + " > 1 Cluster Particles " + NReal.ToString() + " >1 Hamiltonian " + PWHammyReal.ToString("E4") + " Width-x " 
                    + Width_xReal_Total.ToString("E4") + " Width-y " + Width_yReal_Total.ToString("E4") + " Followed by Count Interval Averages preceded by #Clusters(#Peaks)\n";
                for (int Category = 0; Category < NumIntervals; Category++)
                {
                    string start = "Rest ";
                    if (Category < (NumIntervals - 1))
                        start = "Up to " + Limits[Category] + " : ";
                    nextline += start + NumClustersinCategory[Category].ToString() + "(" + NumPeaksinCategory[Category].ToString() + ") Width-x " 
                        + Width_xReal[Category].ToString("E4") + " Width-y " + Width_yReal[Category].ToString("E4") + " * ";
                }
                nextline += "\n";
                int ActualMaxCount = MaxCount + 1;
                nextline += TrimHistograms(HistCounts, ref ActualMaxCount, 0, 1);
                for (int histloop = 0; histloop < ActualMaxCount; histloop++)
                    nextline += HistCounts[histloop] + ", ";
                DAVectorUtility.SALSAPrint(0, nextline);
            }

        }   // End Statistics

        public static string TrimHistograms(int[] CountList, ref int NumCounts, double start, double interval)
        {   // Remove zeros at end of Histograms

            string summary = "";
            int LocalSize = NumCounts;
            for (int CountBackwards = (LocalSize - 1); CountBackwards >= 0; CountBackwards--)
            {
                if (CountList[CountBackwards] > 0)
                    break;
                --NumCounts;
            }
            if(NumCounts == 0)
            {
                summary = "Empty Histogram ";
                NumCounts = 1;
                return summary;
            }
            summary = "Start " + start.ToString("F4") + " Interval " + interval.ToString("F4") + " End " + (start + (NumCounts - 1) * interval).ToString("F4") + " * ";
            return summary;

        }   // End TrimHistograms(int[] CountList, ref int NumCounts, double start, double interval)

        public static string TrimHistograms(int[] CountList, ref int NumCounts, int start, int interval)
        {   // Remove zeros at end of Histograms

            string summary = "";
            int LocalSize = NumCounts;
            for (int CountBackwards = (LocalSize - 1); CountBackwards >= 0; CountBackwards--)
            {
                if (CountList[CountBackwards] > 0)
                    break;
                --NumCounts;
            }
            if (NumCounts == 0)
            {
                summary = "Empty Histogram ";
                NumCounts = 1;
                return summary;
            }
            summary = "Start " + start.ToString() + " Interval " + interval.ToString() + " End " + (start + (NumCounts - 1) * interval).ToString() + " * ";
            return summary;

        }   // End TrimHistograms(int[] CountList, ref int NumCounts, int start, int interval)

        public void SetWidths(int ClusterNumber, out double Width_x, out double Width_y)
        {   // Calculate Contribution to the x and y widths of this cluster. This is NOT divided by Occupation Count

            Width_x = 0.0;
            Width_y = 0.0;
            double[] Center = new double[2];
            double[] Sigma = new double[2];
            int PointsinCluster = 0;

            for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
            {
                int ClusterforPoint = this.PointstoClusterIDs[GlobalPointIndex];
                if (ClusterforPoint != ClusterNumber)
                    continue;
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    Center[VectorIndex] += GoldenExamination.PeakPosition[GlobalPointIndex][VectorIndex];
                ++PointsinCluster;
            }
            if (PointsinCluster == 0)
                return;

            Center[0] = Center[0] / PointsinCluster;
            Center[1] = Center[1] / PointsinCluster;
            Sigma[0] = Program.SigmaVectorParameters_i_[0] * Center[0];
            Sigma[1] = Program.SigmaVectorParameters_i_[1];

            for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
            {
                int ClusterforPoint = this.PointstoClusterIDs[GlobalPointIndex];
                if (ClusterforPoint != ClusterNumber)
                    continue;
                double tmp = (GoldenExamination.PeakPosition[GlobalPointIndex][0] - Center[0]) / Sigma[0];
                Width_x += tmp * tmp;
                tmp = (GoldenExamination.PeakPosition[GlobalPointIndex][1] - Center[1]) / Sigma[1];
                Width_y += tmp * tmp;
            }

        }   // End SetWidths

        public void Difference(ArbitraryClustering BaseClusters)
        {
            int Sizelimit = 0;
            for (int skipcount = 0; skipcount < 6; skipcount++)
            {
                int NumberTargetPeaksinMajorClusters = 0;
                int NumberBasePeaks = 0;
                int NumberBaseClustersUsed = 0;
                if (skipcount == 1)
                    Sizelimit = 2;
                if (skipcount == 2)
                    Sizelimit = 5;
                if (skipcount == 3)
                    Sizelimit = 10;
                if (skipcount == 4)
                    Sizelimit = 20;
                if (skipcount == 5)
                    Sizelimit = 50;

                // UsedTargetClusters = 0 This cluster has no overlap with Golden;  = 1 has overlap but not largest cluster; = 2 matched
                int[] UsedTargetClusters = new int[this.MaxIndependentClusterIndex];
                for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
                    UsedTargetClusters[targetindex] = 0;

                //  Loop over Clusters in Base Cluster set
                for (int BaseClusterLoop = 0; BaseClusterLoop < BaseClusters.MaxIndependentClusterIndex; BaseClusterLoop++)
                {
                    if (BaseClusters.ClusterCountsbyIndex[BaseClusterLoop] <= Sizelimit)
                        continue;
                    ++NumberBaseClustersUsed;

                    //  UsedTargetClusters1 is number of golden peaks of this cluster in given other cluster 
                    int RawClusterLabel = BaseClusters.ClusterIndextoID[BaseClusterLoop];
                    int[] UsedTargetClusters1 = new int[this.MaxIndependentClusterIndex];
                    for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
                        UsedTargetClusters1[targetindex] = 0;
                    NumberBasePeaks += BaseClusters.ClusterCountsbyIndex[BaseClusterLoop];

                    for (int GlobalPointIndex = 0; GlobalPointIndex < BaseClusters.NumberofPoints; GlobalPointIndex++)
                    {
                        if (BaseClusters.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
                            continue;
                        int RawtargetCluster = this.PointstoClusterIDs[GlobalPointIndex];
                        ++UsedTargetClusters1[this.ClusterIDtoIndex[RawtargetCluster]];
                    }

                    int MaxMatch = 0;
                    int TargetClusterMatch = -1;
                    for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
                    {
                        if (UsedTargetClusters1[targetindex] == 0)
                            continue;
                        if (MaxMatch < UsedTargetClusters1[targetindex])
                        {
                            MaxMatch = UsedTargetClusters1[targetindex];
                            TargetClusterMatch = targetindex;
                        }
                        UsedTargetClusters[targetindex] = Math.Max(1, UsedTargetClusters[targetindex]);
                    }
                    NumberTargetPeaksinMajorClusters += MaxMatch;
                    UsedTargetClusters[TargetClusterMatch] = 2;
                }
                int NumberNotOne = 0;
                int NumberOne = 0;
                int NumberMajor = 0;
                for (int targetindex = 0; targetindex < this.MaxIndependentClusterIndex; targetindex++)
                {
                    if (UsedTargetClusters[targetindex] == 0)
                        continue;
                    if (UsedTargetClusters[targetindex] == 2)
                        ++NumberMajor;
                    if (UsedTargetClusters[targetindex] == 1)
                    {
                        if (this.ClusterCountsbyIndex[targetindex] == 1)
                            ++NumberOne;
                        else
                            ++NumberNotOne;
                    }
                }

                DAVectorUtility.SALSAPrint(0, "Base= " + BaseClusters.ClusterString + " Compared to " + this.ClusterString + " Count Limit " + Sizelimit.ToString() + " Base Peaks " + NumberBasePeaks + " In Major Cluster " + NumberTargetPeaksinMajorClusters.ToString() +
                    " Base Clusters " + NumberBaseClustersUsed.ToString() + " " + " Mapped to Major Clusters " + NumberMajor.ToString() +
                    " # Stray >1 " + NumberNotOne.ToString() + " # Stray Clusters 1 member " + NumberOne.ToString());
            }

        }   // end Difference

        public int[][] GeneralClusterOneDHistogram; // [x,y][Bin] Histogram of General Clusters for 1D Projections
        public int[] GeneralClusterTwoDHistogram; // [Bin] Histogram of General Clusters for radial distance

        public void HistogramPeaks(int ClusterCountCut)
        {
            GeneralClusterOneDHistogram = new int[Program.ParameterVectorDimension][];
            GeneralClusterTwoDHistogram = new int[GoldenExamination.TwoDHistogramSize];
            for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                GeneralClusterOneDHistogram[VectorIndex] = new int[GoldenExamination.OneDHistogramSize];

            int NumberHistogrammed = 0;
            for (int ThisClusterLoop = 0; ThisClusterLoop < this.MaxIndependentClusterIndex; ThisClusterLoop++)
            {   // Loop over all clusters skipping small ones
                if (this.ClusterCountsbyIndex[ThisClusterLoop] <= ClusterCountCut)
                    continue;

                int RawClusterLabel = this.ClusterIndextoID[ThisClusterLoop];
                int TotalThisPeaks = this.ClusterCountsbyIndex[ThisClusterLoop];
                ++NumberHistogrammed;

                //  Extract Points for this cluster
                double[][] ThisPeakPositions = new double[this.ClusterCountsbyIndex[ThisClusterLoop]][];
                for (int LocalPointIndex = 0; LocalPointIndex < this.ClusterCountsbyIndex[ThisClusterLoop]; LocalPointIndex++)
                    ThisPeakPositions[LocalPointIndex] = new double[Program.ParameterVectorDimension];
                int localcount = 0;
                for (int GlobalPointIndex = 0; GlobalPointIndex < this.NumberofPoints; GlobalPointIndex++)
                {
                    if (this.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
                        continue;
                    ThisPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
                    ThisPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
                    ++localcount;
                }

                //  Set Sigmas
                double[] Sigma = new double[Program.ParameterVectorDimension];
                double[] Center = new double[Program.ParameterVectorDimension];
                this.SetCenter_Sigma(ThisPeakPositions, TotalThisPeaks, ref Center, ref Sigma);

                GoldenExamination.SetHistograms(ThisPeakPositions, TotalThisPeaks, Center, ref Sigma,
                    GeneralClusterOneDHistogram, GeneralClusterTwoDHistogram);

            }   // End loop over clusters

            DAVectorUtility.SALSAPrint(0, "\n" + this.ClusterString + " Point-Center Histograms for Occupation Count Cuts: 1D Hist Size " + GoldenExamination.OneDHistogramSize.ToString() + " 1D Hist Interval " + GoldenExamination.OneDHistogramInterval.ToString("F3")
                + " Max Val " + GoldenExamination.OneDHistogramMaxVal.ToString("F3") + " 2D Hist Size " + GoldenExamination.TwoDHistogramSize.ToString()
                + " 2D Hist Interval " + GoldenExamination.TwoDHistogramInterval.ToString("F3") + " Max Val " + GoldenExamination.TwoDHistogramMaxVal.ToString("F3"));
            string tempstring = this.ClusterString + " " + NumberHistogrammed.ToString() + " Clusters with occupation count cut Greater Than " + ClusterCountCut.ToString() + "\n";
            int ActualHistSize = GoldenExamination.TwoDHistogramSize;
            tempstring += " 2D " + this.ClusterString + " Histogram "+ TrimHistograms(GeneralClusterTwoDHistogram, ref ActualHistSize, 0.0, GoldenExamination.TwoDHistogramInterval); 
            for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                tempstring += GeneralClusterTwoDHistogram[HistBin].ToString() + ", ";
            DAVectorUtility.SALSAPrint(0, tempstring);
            for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
            {
                tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + this.ClusterString + " " + " Histogram ";
                ActualHistSize = GoldenExamination.OneDHistogramSize;
                tempstring += TrimHistograms(GeneralClusterOneDHistogram[VectorIndex], ref ActualHistSize, 0.0, GoldenExamination.OneDHistogramInterval); 
                for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                    tempstring += GeneralClusterOneDHistogram[VectorIndex][HistBin].ToString() + ", ";
                DAVectorUtility.SALSAPrint(0, tempstring);
            }

        }   // End HistogramPeaks(int ClusterCountCut)

        public void SetCenter_Sigma(double[][] PointPositions, int NumberofPoints, ref double[] Center, ref double[] Sigma)
        {
            for (int ClusterLoop = 0; ClusterLoop < NumberofPoints; ClusterLoop++)
            {
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    Center[VectorIndex] += PointPositions[ClusterLoop][VectorIndex];
            }
            Center[0] = Center[0] / NumberofPoints;
            Center[1] = Center[1] / NumberofPoints;
            Sigma[0] = Program.SigmaVectorParameters_i_[0] * Center[0];
            Sigma[1] = Program.SigmaVectorParameters_i_[1];

        }   // End SetCenterSigma

    }   // End ArbitraryClustering

    public class GoldenExamination
    {
        public static int[] GoldenID;   // Integer ID of Golden records; -1 if not a Golden Point
        public static string[] GoldenLabel;   // String ID of Golden records; "" if not a Golden Point
        public static double[][] PeakPosition;  // [GlobalIndex][x,y] Peak Positions read in
        public static int NumberClusteringMethods;    // Number of Clustering Methods
        public static int OneDHistogramSize;   // Size of 1D Point-Cut Center Histograms including overflow
        public static int TwoDHistogramSize;   // Size of 2D Point-Cut Center Histograms including overflow
        public static double OneDHistogramInterval;   // Interval of 1D Point-Cut Center Histograms 
        public static double TwoDHistogramInterval;   // Interval of 2D Point-Cut Center Histograms
        public static double OneDHistogramMaxVal; // Upper Limit of 1D Point-Cut Center Histograms
        public static double TwoDHistogramMaxVal; // Upper Limit of 1D Point-Cut Center Histograms
        public static string[] CoordLabels = { "m/z", "RT " };

        public static int OneD_Center_HistogramSize;   // Size of 1D Method Center - Golden Center Histograms including overflow
        public static int TwoD_Center_HistogramSize;   // Size of 2D Method Center - Golden Center Histograms including overflow
        public static double OneD_Center_HistogramInterval;   // Interval of 1D Method Center - Golden Center Histograms 
        public static double TwoD_Center_HistogramInterval;   // Interval of 2D Method Center - Golden Center Histograms
        public static double OneD_Center_HistogramMaxVal; // Upper Limit of 1D Method Center - Golden Center Histograms
        public static double TwoD_Center_HistogramMaxVal; // Upper Limit of 1D Method Center - Golden Center Histograms

        public double CutValue;   // Cut on this
        public int MinimumforAveraging;  // Minimum count for an average to be performed
        public int CutMinimumforAveraging;  // Minimum count for an average to be performed
        public int HandleGoldenCluster = 0;  // How to set target Center =0 Average of all peaks
        public int NumberGoldenAccumulations = 0;
        public int NumberGoldenAccumulationsCenterTotalCalc = 0;
        public int NumberGoldenAccumulationsCenterCutCalc = 0;
        public int[] NumberOtherAccumulationsCenterTotalCalc;
        public int[] NumberOtherAccumulationsCenterCutCalc;
        public int[][][] Other_CenterOneDHistogram; // [Method][x,y][Bin] Histogram of Centers for 1D Projections
        public int[][] Other_CenterTwoDHistogram; // [Method][Bin] Histogram of Centers for radial distance

        public GoldenClusterRecord Accumulation;

        public GoldenExamination()
        {
            NumberGoldenAccumulations = 0;
            NumberGoldenAccumulationsCenterTotalCalc = 0;
            NumberGoldenAccumulationsCenterCutCalc = 0;
            NumberOtherAccumulationsCenterTotalCalc = new int[NumberClusteringMethods];
            NumberOtherAccumulationsCenterCutCalc = new int[NumberClusteringMethods];
            for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
            {
                NumberOtherAccumulationsCenterTotalCalc[ClusteringMethod] = 0;
                NumberOtherAccumulationsCenterCutCalc[ClusteringMethod] = 0;
            }
            Accumulation = new GoldenClusterRecord();
            Accumulation.Label = "Accumulation";
            Accumulation.ID = -1;
            Accumulation.TotalGoldenPeaks = 0;
            Accumulation.TotalGoldenPeaksinCut = 0;
            Accumulation.TotalGoldenPeaksoutofCut = 0;

            Other_CenterOneDHistogram = new int[NumberClusteringMethods][][];
            Other_CenterTwoDHistogram = new int[NumberClusteringMethods][];

            for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
            {
                Other_CenterOneDHistogram[ClusteringMethod] = new int[Program.ParameterVectorDimension][];
                Other_CenterTwoDHistogram[ClusteringMethod] = new int[TwoD_Center_HistogramSize];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    Other_CenterOneDHistogram[ClusteringMethod][VectorIndex] = new int[OneD_Center_HistogramSize];
            }
        }

        public class GoldenClusterRecord
        {
            public string Label;
            public int ID;
            public int TotalGoldenPeaks;   // Number of Golden peaks input for this cluster
            public int TotalGoldenPeaksinCut;  // Number of Golden Peaks in Cut
            public int TotalGoldenPeaksoutofCut;   // Number of Golden Peaks out of cut
            public int[] Other_TotalGoldenPeaks;   // [Method] Number of Golden peaks for each Method for this cluster
            public int[] Other_TotalGoldenPeaksinCut;  // [Method] Number of Golden Peaks for each Method in Cut
            public int[] Other_TotalGoldenPeaksoutofCut;   // [Method] Number of Golden Peaks for each Method out of cut but in other cluster
            public int[] Other_TotalGoldenPeaksoutofCluster;   // [Method] Number of Golden Peaks for each Method out of cluster
            public int[] Other_TotalNonGoldenPeaksinCluster;   // [Method] Number of Non Golden Peaks for each Method in cluster
            public int[] Other_TotalNonGoldenPeaksinCut;   // [Method] Number of Non Golden Peaks for each Method in cut
            public double[] GoldenClusterSigma;    // [x,y] Standard Deviation -- NOT squared
            public double[] GoldenClusterCenterTotal;  // [x,y] Golden Center averaging all peaks
            public double[] GoldenClusterCenterCut;  // [x,y] Golden Center averaging peaks in Cut
            public double[][] Other_ClusterCenterDifferenceTotal;  // [Method][x,y] Center for Method in Golden Peaks averaging all peaks -- normalized ABSSOLUTE difference
            public double[][] Other_ClusterCenterDifferenceCut;  // [Method][x,y] Center Method in Golden Peaks averaging peaks in Cut -- normalized ABSOLUTE difference
            public int[][] GoldenClusterOneDHistogram; // [x,y][Bin] Histogram of Golden Clusters for 1D Projections
            public int[] GoldenClusterTwoDHistogram; // [Bin] Histogram of Golden Clusters for radial distance
            public int[][][] Other_ClusterOneDHistogram; // [Method][x,y][Bin] Histogram of Golden Clusters for 1D Projections
            public int[][] Other_ClusterTwoDHistogram; // [Method][Bin] Histogram of Golden Clusters for radial distance


            public GoldenClusterRecord()
            {

                Other_TotalGoldenPeaks = new int[NumberClusteringMethods];
                Other_TotalGoldenPeaksinCut = new int[NumberClusteringMethods];
                Other_TotalGoldenPeaksoutofCut = new int[NumberClusteringMethods];
                Other_TotalGoldenPeaksoutofCluster = new int[NumberClusteringMethods];
                Other_TotalNonGoldenPeaksinCluster = new int[NumberClusteringMethods];
                Other_TotalNonGoldenPeaksinCut = new int[NumberClusteringMethods];

                // Next 3 not used in an accumulation
                GoldenClusterSigma = new double[Program.ParameterVectorDimension];
                GoldenClusterCenterTotal = new double[Program.ParameterVectorDimension];
                GoldenClusterCenterCut = new double[Program.ParameterVectorDimension];

                Other_ClusterCenterDifferenceTotal = new double[NumberClusteringMethods][];
                Other_ClusterCenterDifferenceCut = new double[NumberClusteringMethods][];

                GoldenClusterOneDHistogram = new int[Program.ParameterVectorDimension][];
                GoldenClusterTwoDHistogram = new int[TwoDHistogramSize];
                Other_ClusterOneDHistogram = new int[NumberClusteringMethods][][];
                Other_ClusterTwoDHistogram = new int[NumberClusteringMethods][];

                for (int ClusteringMethod = 0; ClusteringMethod < NumberClusteringMethods; ClusteringMethod++)
                {
                    Other_ClusterCenterDifferenceTotal[ClusteringMethod] = new double[Program.ParameterVectorDimension];
                    Other_ClusterCenterDifferenceCut[ClusteringMethod] = new double[Program.ParameterVectorDimension];
                    Other_ClusterOneDHistogram[ClusteringMethod] = new int[Program.ParameterVectorDimension][];
                    Other_ClusterTwoDHistogram[ClusteringMethod] = new int[TwoDHistogramSize];
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                        Other_ClusterOneDHistogram[ClusteringMethod][VectorIndex] = new int[OneDHistogramSize];
                }
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    GoldenClusterOneDHistogram[VectorIndex] = new int[OneDHistogramSize];

            }   // End Initializer for GoldenClusterRecord

        }   // End GoldenClusterRecord

        // Look at Individual Golden Clusters and compares with 3 Clustering Methods
        //  Outputs results for each cluster and accumulations (averages)
        public void GoldenComparisons(ArbitraryClustering GoldenBase, ArbitraryClustering[] Methods)
        {
            for (int BaseClusterLoop = 0; BaseClusterLoop < GoldenBase.MaxIndependentClusterIndex; BaseClusterLoop++)
            {   // Loop over Golden Clusters skipping small ones
                if (GoldenBase.ClusterCountsbyIndex[BaseClusterLoop] <= this.MinimumforAveraging)
                    continue;

                GoldenClusterRecord GoldenRecordforthisCluster = new GoldenClusterRecord();
                int RawClusterLabel = GoldenBase.ClusterIndextoID[BaseClusterLoop];
                GoldenRecordforthisCluster.TotalGoldenPeaks = GoldenBase.ClusterCountsbyIndex[BaseClusterLoop];
                ++NumberGoldenAccumulations;

                //  Extract Golden Points for this cluster
                int[] ListofGoldenPoints = new int[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]];
                int[] AlwaysGolden = new int[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]];
                double[][] GoldenPeakPositions = new double[GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]][];
                for (int LocalPointIndex = 0; LocalPointIndex < GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]; LocalPointIndex++)
                    GoldenPeakPositions[LocalPointIndex] = new double[Program.ParameterVectorDimension];
                int localcount = 0;
                for (int GlobalPointIndex = 0; GlobalPointIndex < GoldenBase.NumberofPoints; GlobalPointIndex++)
                {
                    if (GoldenBase.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
                        continue;
                    AlwaysGolden[localcount] = 1;
                    ListofGoldenPoints[localcount] = GlobalPointIndex;
                    GoldenPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
                    GoldenPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
                    GoldenRecordforthisCluster.ID = GoldenExamination.GoldenID[GlobalPointIndex];   // Same for all entries
                    GoldenRecordforthisCluster.Label = GoldenExamination.GoldenLabel[GlobalPointIndex];
                    ++localcount;
                }

                //  Set Sigmas
                double[] GoldenSigma = new double[Program.ParameterVectorDimension];
                this.SetSigma(GoldenPeakPositions, GoldenRecordforthisCluster.TotalGoldenPeaks, ref GoldenSigma);

                //  Set Cut for Golden Peaks
                int dum1;
                int dum2;
                this.SetCutStatus(GoldenPeakPositions, AlwaysGolden, GoldenRecordforthisCluster.TotalGoldenPeaks, ref GoldenSigma,
                    out GoldenRecordforthisCluster.TotalGoldenPeaksinCut,
                    out GoldenRecordforthisCluster.TotalGoldenPeaksoutofCut, out dum1, out dum2,
                    ref GoldenRecordforthisCluster.GoldenClusterCenterTotal, ref GoldenRecordforthisCluster.GoldenClusterCenterCut);

                //  Basic Accumulations
                ++NumberGoldenAccumulationsCenterTotalCalc;
                this.Accumulation.TotalGoldenPeaks += GoldenRecordforthisCluster.TotalGoldenPeaks;
                this.Accumulation.TotalGoldenPeaksinCut += GoldenRecordforthisCluster.TotalGoldenPeaksinCut;
                this.Accumulation.TotalGoldenPeaksoutofCut += GoldenRecordforthisCluster.TotalGoldenPeaksoutofCut;

                if (GoldenRecordforthisCluster.TotalGoldenPeaksinCut <= this.CutMinimumforAveraging)
                    continue;

                SetHistograms(GoldenPeakPositions, GoldenRecordforthisCluster.TotalGoldenPeaks, GoldenRecordforthisCluster.GoldenClusterCenterCut, ref GoldenSigma,
                    GoldenRecordforthisCluster.GoldenClusterOneDHistogram, GoldenRecordforthisCluster.GoldenClusterTwoDHistogram);

                //  Golden Histogram Accumulations
                ++NumberGoldenAccumulationsCenterCutCalc;
                for (int Histbin = 0; Histbin < TwoDHistogramSize; Histbin++)
                    this.Accumulation.GoldenClusterTwoDHistogram[Histbin] += GoldenRecordforthisCluster.GoldenClusterTwoDHistogram[Histbin];
                for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                {
                    for (int Histbin = 0; Histbin < OneDHistogramSize; Histbin++)
                        this.Accumulation.GoldenClusterOneDHistogram[VectorIndex][Histbin] += GoldenRecordforthisCluster.GoldenClusterOneDHistogram[VectorIndex][Histbin];
                }

                //  Analyze Other Clustering Methods
                for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
                {
                    ArbitraryClustering ClusterMethod = Methods[MethodLoop];
                    int[] UsedTargetClusters1 = new int[ClusterMethod.MaxIndependentClusterIndex];
                    for (int targetindex = 0; targetindex < ClusterMethod.MaxIndependentClusterIndex; targetindex++)
                        UsedTargetClusters1[targetindex] = 0;
                    for (int LocalPointIndex = 0; LocalPointIndex < GoldenBase.ClusterCountsbyIndex[BaseClusterLoop]; LocalPointIndex++)
                    {
                        int GlobalPointIndex = ListofGoldenPoints[LocalPointIndex];
                        int RawtargetCluster = ClusterMethod.PointstoClusterIDs[GlobalPointIndex];
                        ++UsedTargetClusters1[ClusterMethod.ClusterIDtoIndex[RawtargetCluster]];
                    }

                    //  MaxMatch is number of Matching peaks
                    //  TargetClusterMatch is index of Method Cluster matching this Golden Cluster
                    int MaxMatch = 0;
                    int TargetClusterMatch = -1;
                    for (int targetindex = 0; targetindex < ClusterMethod.MaxIndependentClusterIndex; targetindex++)
                    {
                        if (UsedTargetClusters1[targetindex] == 0)
                            continue;
                        if (MaxMatch < UsedTargetClusters1[targetindex])
                        {
                            MaxMatch = UsedTargetClusters1[targetindex];
                            TargetClusterMatch = targetindex;
                        }
                    }
                    if (MaxMatch <= this.CutMinimumforAveraging)
                        continue;

                    // Set Points of Clustering Method in this Cluster
                    int OtherRawClusterLabel = ClusterMethod.ClusterIndextoID[TargetClusterMatch];
                    int num = ClusterMethod.ClusterCountsbyIndex[TargetClusterMatch];
                    int[] ListofOtherPoints = new int[num];
                    int[] GoldenorNot = new int[num];
                    double[][] OtherPeakPositions = new double[num][];
                    for (int LocalPointIndex = 0; LocalPointIndex < num; LocalPointIndex++)
                        OtherPeakPositions[LocalPointIndex] = new double[Program.ParameterVectorDimension];
                    localcount = 0;
                    int WrongGoldenCluster = 0;
                    for (int GlobalPointIndex = 0; GlobalPointIndex < GoldenBase.NumberofPoints; GlobalPointIndex++)
                    {
                        if (ClusterMethod.PointstoClusterIDs[GlobalPointIndex] != OtherRawClusterLabel)
                            continue;
                        GoldenorNot[localcount] = 1;
                        if (GoldenBase.PointstoClusterIDs[GlobalPointIndex] != RawClusterLabel)
                        {
                            ++WrongGoldenCluster;
                            GoldenorNot[localcount] = 0;
                        }
                        ListofOtherPoints[localcount] = GlobalPointIndex;
                        OtherPeakPositions[localcount][0] = GoldenExamination.PeakPosition[GlobalPointIndex][0];
                        OtherPeakPositions[localcount][1] = GoldenExamination.PeakPosition[GlobalPointIndex][1];
                        ++localcount;
                    }
                    GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCluster[MethodLoop] = WrongGoldenCluster;
                    GoldenRecordforthisCluster.Other_TotalGoldenPeaks[MethodLoop] = MaxMatch;
                    GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCluster[MethodLoop] = GoldenRecordforthisCluster.TotalGoldenPeaks - MaxMatch;

                    //  Set Cut for Other Method
                    int Totin = 0;
                    int Totout = 0;
                    double[] FullCenterfromClustering = new double[Program.ParameterVectorDimension];
                    double[] CutCenterfromClustering = new double[Program.ParameterVectorDimension];
                    this.SetCutStatus(OtherPeakPositions, GoldenorNot, num, ref GoldenSigma, out Totin, out Totout,
                        out GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop], out GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCut[MethodLoop],
                        ref FullCenterfromClustering, ref CutCenterfromClustering);
                    GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCut[MethodLoop] = Totin - GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop];

                    //  Accumulate
                    ++NumberOtherAccumulationsCenterTotalCalc[MethodLoop];

                    this.Accumulation.Other_TotalGoldenPeaks[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaks[MethodLoop];
                    this.Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksinCut[MethodLoop];
                    this.Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCut[MethodLoop];
                    this.Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop] += GoldenRecordforthisCluster.Other_TotalGoldenPeaksoutofCluster[MethodLoop];
                    this.Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop] += GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCluster[MethodLoop];
                    this.Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop] += GoldenRecordforthisCluster.Other_TotalNonGoldenPeaksinCut[MethodLoop];

                    if (Totin <= this.CutMinimumforAveraging)
                        continue;

                    double distance = 0.0;
                    int idist;
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        double tmp = (CutCenterfromClustering[VectorIndex] - GoldenRecordforthisCluster.GoldenClusterCenterCut[VectorIndex]) / GoldenSigma[VectorIndex];
                        GoldenRecordforthisCluster.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] = Math.Abs(tmp);
                        distance += tmp * tmp;
                        idist = (int)Math.Floor(Math.Abs(tmp) / GoldenExamination.OneD_Center_HistogramInterval);
                        idist = Math.Min(idist, GoldenExamination.OneD_Center_HistogramSize - 1);
                        ++this.Other_CenterOneDHistogram[MethodLoop][VectorIndex][idist];

                        double tmp1 = (FullCenterfromClustering[VectorIndex] - GoldenRecordforthisCluster.GoldenClusterCenterTotal[VectorIndex]) / GoldenSigma[VectorIndex];
                        GoldenRecordforthisCluster.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] = Math.Abs(tmp1);
                    }
                    idist = (int)Math.Floor(Math.Abs(distance) / GoldenExamination.TwoD_Center_HistogramInterval);
                    idist = Math.Min(idist, GoldenExamination.TwoD_Center_HistogramSize - 1);
                    ++this.Other_CenterTwoDHistogram[MethodLoop][idist];

                    SetHistograms(OtherPeakPositions, num, CutCenterfromClustering, ref GoldenSigma,
                        GoldenRecordforthisCluster.Other_ClusterOneDHistogram[MethodLoop], GoldenRecordforthisCluster.Other_ClusterTwoDHistogram[MethodLoop]);

                    //  Final accumulations
                    ++NumberOtherAccumulationsCenterCutCalc[MethodLoop];
                    for (int VectorIndex = 0; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                    {
                        this.Accumulation.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] += GoldenRecordforthisCluster.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex];
                        this.Accumulation.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] += GoldenRecordforthisCluster.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex];
                        for (int Histbin = 0; Histbin < OneDHistogramSize; Histbin++)
                            this.Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][Histbin] += GoldenRecordforthisCluster.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][Histbin];
                    }
                    for (int Histbin = 0; Histbin < TwoDHistogramSize; Histbin++)
                        this.Accumulation.Other_ClusterTwoDHistogram[MethodLoop][Histbin] += GoldenRecordforthisCluster.Other_ClusterTwoDHistogram[MethodLoop][Histbin];


                }   // End loop over MethodLoop

            }   // End loop over Golden Clusters

        }   // End GoldenComparisons

        public void SetSigma(double[][] PointPositions, int NumberofPoints, ref double[] Sigma)
        {
            double[] center = new double[2];
            for (int ClusterLoop = 0; ClusterLoop < NumberofPoints; ClusterLoop++)
            {
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    center[VectorIndex] += PointPositions[ClusterLoop][VectorIndex];
            }
            center[0] = center[0] / NumberofPoints;
            center[1] = center[1] / NumberofPoints;
            Sigma[0] = Program.SigmaVectorParameters_i_[0] * center[0];
            Sigma[1] = Program.SigmaVectorParameters_i_[1];

        }   // End SetSigma

        //  Analyze a cluster
        //  Goldenstatus =1 if point Golden
        public void SetCutStatus(double[][] PointPositions, int[] GoldenStatus, int NumberofPoints, ref double[] Sigma, out int NumberinCut, out int NumberoutofCut,
            out int NumberGoldeninCut, out int NumberGoldenoutofCut, ref double[] FullCenter, ref double[] CutCenter)
        {
            NumberoutofCut = 0;
            NumberinCut = 0;
            NumberGoldeninCut = 0;
            NumberGoldenoutofCut = 0;

            for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                FullCenter[VectorIndex] = 0.0;
            int TotalGoldeninCluster = 0;
            for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
            {
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    FullCenter[VectorIndex] += PointPositions[ClusterPointLoop][VectorIndex];
                if (GoldenStatus[ClusterPointLoop] == 1)
                    ++TotalGoldeninCluster;
            }
            for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
            {
                FullCenter[VectorIndex] = FullCenter[VectorIndex] / NumberofPoints;
                CutCenter[VectorIndex] = FullCenter[VectorIndex];
            }

            // Iterate points in cut
            double[] NewCenter = new double[Program.ParameterVectorDimension];

            for (int IterationLoop = 0; IterationLoop < 10; IterationLoop++)
            {
                double distance;
                NumberinCut = 0;
                NumberGoldeninCut = 0;
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    NewCenter[VectorIndex] = 0.0;
                for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
                {
                    double tmp;
                    distance = 0.0;
                    int iout = 0;
                    for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    {
                        tmp = (PointPositions[ClusterPointLoop][VectorIndex] - CutCenter[VectorIndex]) / Sigma[VectorIndex];
                        if (Math.Abs(tmp) > this.CutValue)
                            iout = 1;
                        distance += tmp * tmp;
                    }
                    if ((this.HandleGoldenCluster == 1) && (iout == 1))
                        continue;
                    if ((this.HandleGoldenCluster == 0) && (distance > this.CutValue))
                        continue;
                    ++NumberinCut;
                    if (GoldenStatus[ClusterPointLoop] == 1)
                        ++NumberGoldeninCut;
                    for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                        NewCenter[VectorIndex] += PointPositions[ClusterPointLoop][VectorIndex];
                }
                if (NumberinCut <= 0)
                    break;
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                    CutCenter[VectorIndex] = NewCenter[VectorIndex] / NumberinCut;
            }
            NumberoutofCut = NumberofPoints - NumberinCut;
            NumberGoldenoutofCut = TotalGoldeninCluster - NumberGoldeninCut;
            return;

        }   // End SetCutStatus

        //  Set Histograms -- each point is put in three histograms
        public static void SetHistograms(double[][] PointPositions, int NumberofPoints, double[] CutCenter, ref double[] Sigma, int[][] OneDHistogram, int[] TwoDHistogram)
        {
            for (int ClusterPointLoop = 0; ClusterPointLoop < NumberofPoints; ClusterPointLoop++)
            {
                double tmp;
                double distance = 0.0;
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                {
                    tmp = (PointPositions[ClusterPointLoop][VectorIndex] - CutCenter[VectorIndex]) / Sigma[VectorIndex];
                    int idist = (int)Math.Floor(Math.Abs(tmp) / GoldenExamination.OneDHistogramInterval);
                    idist = Math.Min(idist, GoldenExamination.OneDHistogramSize - 1);
                    ++OneDHistogram[VectorIndex][idist];
                    distance += tmp * tmp;
                }
                int idist2D = (int)Math.Floor(distance / GoldenExamination.TwoDHistogramInterval);
                idist2D = Math.Min(idist2D, GoldenExamination.TwoDHistogramSize - 1);
                ++TwoDHistogram[idist2D];
            }

        }   // End SetHistograms

        public void PrintAccumulation()
        {
            string[] ClusteringLabels = { "DAVS  ", "Medea ", "Mclust", "Golden" };

            DAVectorUtility.SALSAPrint(0, "\n*********************** Golden Cluster Analysis Cut Value "
                + CutValue.ToString("F3") + " Minimum needed to Average " + MinimumforAveraging.ToString() + " After Cut " + CutMinimumforAveraging.ToString() + " Handling Option " + HandleGoldenCluster.ToString());

            DAVectorUtility.SALSAPrint(0, "1D Hist Size " + OneDHistogramSize.ToString() + " 1D Hist Interval " + OneDHistogramInterval.ToString("F3") + " Max Val " + OneDHistogramMaxVal.ToString("F3")
                + " 2D Hist Size " + TwoDHistogramSize.ToString() + " 2D Hist Interval " + TwoDHistogramInterval.ToString("F3") + " Max Val " + TwoDHistogramMaxVal.ToString("F3"));

            DAVectorUtility.SALSAPrint(0, " Records Accumulated " + NumberGoldenAccumulations.ToString());
            string tempstring = "";
            for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
                tempstring += " " + ClusteringLabels[MethodLoop] + " " + NumberOtherAccumulationsCenterTotalCalc[MethodLoop].ToString();
            DAVectorUtility.SALSAPrint(0, " Records with Good Centers and Accumulated " + NumberGoldenAccumulationsCenterTotalCalc.ToString() + tempstring);
            tempstring = "";
            for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
                tempstring += " " + ClusteringLabels[MethodLoop] + " " + NumberOtherAccumulationsCenterCutCalc[MethodLoop].ToString();
            DAVectorUtility.SALSAPrint(0, " Records with Good Cut Centers and Histograms " + NumberGoldenAccumulationsCenterCutCalc.ToString() + tempstring);

            double divisor1 = 1.0 / ((double)NumberGoldenAccumulations);
            double tmp1 = Accumulation.TotalGoldenPeaksinCut * divisor1;
            double tmp2 = Accumulation.TotalGoldenPeaksoutofCut * divisor1;
            DAVectorUtility.SALSAPrint(0, " Golden Cut Analysis Cut Value " + CutValue.ToString("F3") + " Minimum needed to Average " + MinimumforAveraging.ToString()
                + " After Cut " + CutMinimumforAveraging.ToString() + " Handling Option " + HandleGoldenCluster.ToString()
                + " Full Set of Golden Cluster Peaks " + Accumulation.TotalGoldenPeaksinCut.ToString()
                + " In Cut " + Accumulation.TotalGoldenPeaksinCut.ToString() + " Per Cluster " + tmp1.ToString("F4")
                + " Out of Cut " + Accumulation.TotalGoldenPeaksoutofCut.ToString() + " Per Cluster " + tmp2.ToString("F4"));
            tempstring = "\n2D Pure Golden Histogram ";
            int ActualHistSize = TwoDHistogramSize;
            tempstring += ArbitraryClustering.TrimHistograms(Accumulation.GoldenClusterTwoDHistogram, ref ActualHistSize, 0.0, TwoDHistogramInterval);
            for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                tempstring += Accumulation.GoldenClusterTwoDHistogram[HistBin].ToString() + ", ";
            DAVectorUtility.SALSAPrint(0, tempstring);
            for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
            {
                tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " Golden Histogram ";
                ActualHistSize = OneDHistogramSize;
                tempstring += ArbitraryClustering.TrimHistograms(Accumulation.GoldenClusterOneDHistogram[VectorIndex], ref ActualHistSize, 0.0, OneDHistogramInterval);
                for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                    tempstring += Accumulation.GoldenClusterOneDHistogram[VectorIndex][HistBin].ToString() + ", ";
                DAVectorUtility.SALSAPrint(0, tempstring);
            }

            //  Loop over Methods
            for (int MethodLoop = 0; MethodLoop < NumberClusteringMethods; MethodLoop++)
            {
                double divisor2 = 1.0 / ((double)NumberOtherAccumulationsCenterTotalCalc[MethodLoop]);
                tmp1 = NumberOtherAccumulationsCenterTotalCalc[MethodLoop] * divisor1;
                tmp2 = NumberOtherAccumulationsCenterCutCalc[MethodLoop] * divisor1;
                DAVectorUtility.SALSAPrint(0, " Golden Cut Analysis Cut Value " + CutValue.ToString("F3") + " Minimum needed to Average " + MinimumforAveraging.ToString()
                    + " After Cut " + CutMinimumforAveraging.ToString() + " Handling Option " + HandleGoldenCluster.ToString()
                    + "\nMethod " + ClusteringLabels[MethodLoop] + " Accumulations " + NumberOtherAccumulationsCenterTotalCalc[MethodLoop].ToString() + " "
                    + tmp1.ToString("F4") + " with Cut and Histogram " + NumberOtherAccumulationsCenterCutCalc[MethodLoop].ToString() + " " + tmp2.ToString("F4"));
                int totalGPeaksinthislot = Accumulation.Other_TotalGoldenPeaks[MethodLoop] + Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop];
                double divisor3 = 1.0 / ((double)totalGPeaksinthislot);
                tmp1 = Accumulation.Other_TotalGoldenPeaks[MethodLoop] * divisor3;
                tmp2 = Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop] * divisor3;
                double tmp3 = Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop] * divisor3;
                DAVectorUtility.SALSAPrint(0, "Method " + ClusteringLabels[MethodLoop] + " Golden Peaks in Clusters "
                    + Accumulation.Other_TotalGoldenPeaks[MethodLoop].ToString() + " " + tmp1.ToString("F4")
                    + " Golden Peaks NOT in Method Cluster at all " + Accumulation.Other_TotalGoldenPeaksoutofCluster[MethodLoop].ToString() + " " + tmp2.ToString("F4") + " False Peaks in Cluster "
                    + Accumulation.Other_TotalNonGoldenPeaksinCluster[MethodLoop].ToString() + " " + tmp3.ToString("F4") + " NOT in divisor");
                tmp1 = Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop] * divisor3;
                tmp2 = Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop] * divisor3;
                tmp3 = Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop] * divisor3;
                DAVectorUtility.SALSAPrint(0, "Method " + ClusteringLabels[MethodLoop] + " Golden Peaks in Cut "
                    + Accumulation.Other_TotalGoldenPeaksinCut[MethodLoop].ToString() + " " + tmp1.ToString("F4")
                    + " Golden Peaks out of cut " + Accumulation.Other_TotalGoldenPeaksoutofCut[MethodLoop].ToString() + " " + tmp2.ToString("F4") + " False Peaks in Cut "
                    + Accumulation.Other_TotalNonGoldenPeaksinCut[MethodLoop].ToString() + " " + tmp3.ToString("F4") + " NOT in divisor");

                double divisor4 = 1.0 / ((double)NumberOtherAccumulationsCenterTotalCalc[MethodLoop]);
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                {
                    tmp1 = Accumulation.Other_ClusterCenterDifferenceTotal[MethodLoop][VectorIndex] * divisor4;
                    tmp2 = Accumulation.Other_ClusterCenterDifferenceCut[MethodLoop][VectorIndex] * divisor4;
                    DAVectorUtility.SALSAPrint(0, GoldenExamination.CoordLabels[VectorIndex] + " Method " + ClusteringLabels[MethodLoop] + " Average discrepancy for Full Cluster " + tmp1.ToString("F4")
                        + " Average discrepancy for Cut Cluster " + tmp2.ToString("F4"));
                }

                tempstring = "\n2D " + ClusteringLabels[MethodLoop] + " Abs(Point - Center) Histogram ";
                ActualHistSize = TwoDHistogramSize;
                tempstring += ArbitraryClustering.TrimHistograms(Accumulation.Other_ClusterTwoDHistogram[MethodLoop], ref ActualHistSize, 0.0, TwoDHistogramInterval);
                for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                    tempstring += Accumulation.Other_ClusterTwoDHistogram[MethodLoop][HistBin].ToString() + ", ";
                DAVectorUtility.SALSAPrint(0, tempstring);
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                {
                    tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + ClusteringLabels[MethodLoop] + " " + " Abs(Point - Center) Histogram ";
                    ActualHistSize = OneDHistogramSize;
                    tempstring += ArbitraryClustering.TrimHistograms(Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex], ref ActualHistSize, 0.0, OneDHistogramInterval);
                    for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                        tempstring += Accumulation.Other_ClusterOneDHistogram[MethodLoop][VectorIndex][HistBin].ToString() + ", ";
                    DAVectorUtility.SALSAPrint(0, tempstring);
                }

                tempstring = "\n2D " + ClusteringLabels[MethodLoop] + " Method-Golden Center Histogram ";
                ActualHistSize = TwoD_Center_HistogramSize;
                tempstring += ArbitraryClustering.TrimHistograms(this.Other_CenterTwoDHistogram[MethodLoop], ref ActualHistSize, 0.0, TwoD_Center_HistogramInterval);
                for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                    tempstring += this.Other_CenterTwoDHistogram[MethodLoop][HistBin].ToString() + ", ";
                DAVectorUtility.SALSAPrint(0, tempstring);
                for (int VectorIndex = 0; VectorIndex < 2; VectorIndex++)
                {
                    tempstring = "1D " + GoldenExamination.CoordLabels[VectorIndex] + " " + ClusteringLabels[MethodLoop] + " " + " Method-Golden Center Histogram ";
                    ActualHistSize = TwoD_Center_HistogramSize;
                    tempstring += ArbitraryClustering.TrimHistograms(this.Other_CenterOneDHistogram[MethodLoop][VectorIndex], ref ActualHistSize, 0.0, OneD_Center_HistogramInterval);
                    for (int HistBin = 0; HistBin < ActualHistSize; HistBin++)
                        tempstring += this.Other_CenterOneDHistogram[MethodLoop][VectorIndex][HistBin].ToString() + ", ";
                    DAVectorUtility.SALSAPrint(0, tempstring);
                }

            }   // End output Loop over Methods

            DAVectorUtility.SALSAPrint(0, "\n***************** End Average Status of Golden Clusters Cut Value " + CutValue.ToString("F3")
                + " Minimum needed to Average " + MinimumforAveraging.ToString() + " After Cut " + CutMinimumforAveraging.ToString() + " Handling Option " + HandleGoldenCluster.ToString() + "\n");

        }   // End PrintAccumulation()

    }   // End GoldenExamination


}   // End Namespace Salsa.DAVectorSponge