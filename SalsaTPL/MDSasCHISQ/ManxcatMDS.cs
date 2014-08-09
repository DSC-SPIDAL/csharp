using System;
using System.Collections;
using System.IO;
using System.Threading.Tasks;
using Manxcat;
using SALSALibrary;

namespace MDS
{
    public class ManxcatMDS
    {
        //  Input Variables
        public static int Chisqnorm = 0; // Option for Chisq Norm
        public static int DistanceFormula = 1; // =2 Square of Euclidean Distance, = 1 EuclideanDistance
        public static int WeightingOption = 0; // Weight Option
        public static int NumLinksHistogramBins = 100; // Bins in Links Histogram
        public static int LinkCutforCenter = 20; // Ensure Center has at least this many links

        //  Calculated Variables
        public static int[] PointStatus;
                            // Set status of points = -1 Deleted and Fixed = 0 normal, = 1 Center, =2 x axis, =3 in x,y plane NOT DECOMPOSED

        public static double[] PointWeights; // Intrinsic weight of each point
        public static double SystemRadius = 0.0; //  Radius of System
        public static double SystemMax = 0.0; // Maximumum value of distance
        public static double SystemAverage = 0.0; //  Average Distance in System
        public static double SystemSigma = 0.0; // Standard Deviation in System
        public static double SQRTSystemAverage = 0.0; //  SQRT Average Distance in System
        public static double MinimumDistance; // Minimum Distance for weight calculation
        public static double SQRTMinimumDistance; // Minimum Distance for weight calculation
        public static double SQUAREMinimumDistance; // Minimum Distance for weight calculation
        public static bool FindPointstoFix = true; // If true find key points to fix

        public static void SetupHotsunforMDS()
        {
            int NumberofPoints = SALSAUtility.PointCount_Global;
            Hotsun.ParameterVectorDimension = ManxcatCentral.Configuration.LocalVectorDimension;
            Hotsun.ndata = NumberofPoints*(NumberofPoints - 1L)/2L;
            Hotsun.Number_VectorParameters = NumberofPoints;
            Hotsun.npar = Hotsun.Number_VectorParameters*Hotsun.ParameterVectorDimension;
            SALSAUtility.sequentialBLAS = false;
            Hotsun.DecomposeParameters = true;
            Hotsun.fullmatrixset = false;

            return;
        }

        // End SetupHotsunforMDS

        //  Set Chisq version of MDS
        public static void SetupMDSasChisq()
        {
            int NumberofPoints = SALSAUtility.PointCount_Global;
            Hotsun.ndata = NumberofPoints*(NumberofPoints - 1L)/2L;

            // Initialize Linear Algebra
            MDSLinearAlgebra.Initialize();

            // Read in distance data
            SALSAUtility.DistanceProcessingOption = ManxcatCentral.Configuration.DistanceProcessingOption;
            SALSAParallelism.ReadDataFromFile(ManxcatCentral.ActualDataFileName);

            // If FindPointstoFix true, all fixed parameters are set to 0 UNLESS initialized differently
            FindPointstoFix = true;

            if (SALSAUtility.NumberFixedPoints > 0)
                FindPointstoFix = false;

            if (ManxcatCentral.Configuration.InitializationOption > 0)
                ManxcatCentral.Configuration.InitializationLoops = 1; // Can't loop if fixed initial positions


            // set up Fixed and Deleted Parameters
            PointStatus = new int[SALSAUtility.PointCount_Global];
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                PointStatus[GlobalPointIndex] = 0;
            }

            SALSAUtility.SetupDistanceWeights();

            if (SALSAUtility.TransformMethod == 10 || (SALSAUtility.NumberDistanceWeightCuts > 0))
            {
                for (int cutloop = 0; cutloop < SALSAUtility.NumberDistanceWeightCuts; cutloop++)
                {
                    SALSAUtility.DistanceWeightCuts[cutloop] = Math.Pow(SALSAUtility.DistanceWeightCuts[cutloop],
                                                                        SALSAUtility.TransformParameter);
                }
            }
            else if (SALSAUtility.TransformMethod == 8 || SALSAUtility.TransformMethod == 9)
            {
                SALSAUtility.NumberDistanceWeightCuts = 0;
            }


            // Loop over initial analysis until no more deleted points
            SALSAUtility.SALSAPrint(1,
                                    "\nInitial Processing Parameters\nDistance Cut " +
                                    SALSAUtility.DistanceCut.ToString("F3") + " Link Cut " +
                                    SALSAUtility.LinkCut.ToString()
                                    + " Allowed Deleted Fraction " + SALSAUtility.AllowedDeletedFraction.ToString("F3") +
                                    " Undefined Distance Value " + SALSAUtility.UndefinedDistanceValue.ToString("F3")
                                    + "\nDistance Transformation Method " + SALSAUtility.TransformMethod.ToString() +
                                    " with Parameter " + SALSAUtility.TransformParameter.ToString("F4"));
            double MissingDistances = 0.0;
            int DisconnectedPoints = 0;
            int DisconnectedLoopCount = 1;
            while (true)
            {
                int oldDisconnectedPoints = DisconnectedPoints;
                SALSAUtility.SALSAPrint(1,
                                        "\n ******* Loop over Identification of Disconnected Points " +
                                        DisconnectedLoopCount.ToString());
                ManxcatMDSBasicDataProcessing.IntialDistanceAnalysis(ref SystemAverage, ref SystemMax, ref SystemSigma,
                                                                     out DisconnectedPoints, out MissingDistances);
                SQRTSystemAverage = Math.Sqrt(SystemAverage);
                if (oldDisconnectedPoints >= DisconnectedPoints)
                    break;
                int AllowedDisconnectedNumber = (int) SALSAUtility.AllowedDeletedFraction +
                                                SALSAUtility.PointCount_Global;
                if (DisconnectedPoints > AllowedDisconnectedNumber)
                {
                    Exception e =
                        SALSAUtility.SALSAError("Must stop as Number of Disconnected Points " +
                                                DisconnectedPoints.ToString() + " Exceeds Cut " +
                                                AllowedDisconnectedNumber.ToString());
                    throw (e);
                }
                ++DisconnectedLoopCount;
            }
            int deletedpoints = 0;
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                if (PointStatus[GlobalPointIndex] == -1)
                {
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        Hotsun.FixedParameter[GlobalPointIndex][LocalVectorIndex] = true;
                    }
                    ++deletedpoints;
                }
            }
            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
                SALSAUtility.SALSAPrint(1, "Deleted Points due to bad link count " + deletedpoints.ToString());

            // If necessary Clean up and Transform Distances overriding previous averages
            ManxcatMDSBasicDataProcessing.CleanandTransformDistances(ref SystemAverage, ref SystemMax, ref SystemSigma,
                                                                     out DisconnectedPoints, out MissingDistances);

            // If ManxcatCentral.MetadataforRun.MinimumDistance  is positive, it is absolute Minimum Distance
            // If ManxcatCentral.MetadataforRun.MinimumDistance  is negative, it is - multiplier of System Average
            if (ManxcatCentral.Configuration.MinimumDistance < 0)
                MinimumDistance = -ManxcatCentral.Configuration.MinimumDistance*SystemAverage;
            else
                MinimumDistance = ManxcatCentral.Configuration.MinimumDistance;
            SQRTMinimumDistance = Math.Sqrt(MinimumDistance);
            SQUAREMinimumDistance = MinimumDistance*MinimumDistance;

            //  Now find centers and histogram distances
            int PointsinDistanceHistogram = ManxcatCentral.Configuration.HistogramBinCount;
            double Histmin = 0.0;
            double Histmax = SystemMax;
            ManxcatMDSBasicDataProcessing.SetUpHistogramRange(PointsinDistanceHistogram, ref Histmin, ref Histmax);
            var Bincounts = new double[2 + PointsinDistanceHistogram];
            int Center = 0;
            ManxcatMDSBasicDataProcessing.FindCenter(ref Center, ref SystemRadius, Histmin, Histmax,
                                                     PointsinDistanceHistogram, ref Bincounts);

            PointStatus[Center] = 1;
            if (FindPointstoFix)
            {
                for (int LocalVectorIndex = 0; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
                {
                    Hotsun.FixedParameter[Center][LocalVectorIndex] = true;
                }
            }

            //  Set up Weights
            PointWeights = new double[SALSAUtility.PointCount_Global];
            WeightingOption = ManxcatCentral.Configuration.WeightingOption;
            SetupWeightings(PointWeights);

            double DistancesNearEachOther = 0.0;
            int NotLonelyPoints = 0;
            int xAxis = 0;
            double xAxisExtent = 0.0;
            ManxcatMDSBasicDataProcessing.FindxAxis(Center, ref xAxis, ref xAxisExtent, ref DistancesNearEachOther,
                                                    ref NotLonelyPoints);
            PointStatus[xAxis] = 2;


            if (FindPointstoFix)
            {
                for (int LocalVectorIndex = 1; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
                {
                    Hotsun.FixedParameter[xAxis][LocalVectorIndex] = true;
                }
            }

            int xyPlane = 0;
            double xyPlaneExtent = 0.0;
            ManxcatMDSBasicDataProcessing.FindxyPlane(Center, xAxis, ref xyPlane, ref xyPlaneExtent);
            PointStatus[xyPlane] = 3;

            if (FindPointstoFix)
            {
                for (int LocalVectorIndex = 2; LocalVectorIndex < Hotsun.ParameterVectorDimension; LocalVectorIndex++)
                {
                    Hotsun.FixedParameter[xyPlane][LocalVectorIndex] = true;
                }
            }

            if (!FindPointstoFix)
            {
                if (SALSAUtility.NumberFixedPoints < (Hotsun.ParameterVectorDimension - 1))
                {
                    Exception e =
                        SALSAUtility.SALSAError("Not enough fixed points " + SALSAUtility.NumberFixedPoints.ToString());
                    throw (e);
                }
            }

            if (SALSAUtility.NumberFixedPoints > 0)
            {
                for (int FixedIndex = 0; FixedIndex < SALSAUtility.NumberFixedPoints; FixedIndex++)
                {
                    int ParameterIndex = SALSAUtility.FixedPointOriginal[FixedIndex];
                    ParameterIndex = SALSAUtility.OriginalPointtoUsedPointMap[ParameterIndex];
                    // This is Used point number
                    if (ParameterIndex < 0)
                    {
                        Exception e = SALSAUtility.SALSAError("Illegal Fixed Point " + FixedIndex.ToString());
                        throw (e);
                    }
                    for (int LocalVectorIndex = 0;
                         LocalVectorIndex < Hotsun.ParameterVectorDimension;
                         LocalVectorIndex++)
                    {
                        Hotsun.FixedParameter[ParameterIndex][LocalVectorIndex] = true;
                    }
                }
            }

            // Normalize Chisq for MDS
            ManxcatCentral.ChisqPrintConstant = 0.5*ManxcatCentral.ChisqPrintConstant/
                                                (Hotsun.ndata - 0.5*MissingDistances - Hotsun.npar);
                // Print Chisq per point (factor of .5 as each point doubled)
            Chisqnorm = ManxcatCentral.Configuration.Chisqnorm;
            DistanceFormula = ManxcatCentral.Configuration.DistanceFormula;

            if (DistanceFormula == 0)
                DistanceFormula = 2;

            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
            {
                double tmp = SALSAUtility.NumberVariedPoints;
                double fraction = MissingDistances/(tmp*(tmp - 1.0));
                double EstimatedDimension = 2.0*SystemAverage*SystemAverage/(SystemSigma*SystemSigma);
                SALSAUtility.SALSAPrint(1,
                                        "\nAFTER CLEAN Disconnected Points " + DisconnectedPoints.ToString() +
                                        " Missing Distances " + MissingDistances.ToString("F0") + " Fraction " +
                                        fraction.ToString("F4"));
                SALSAUtility.SALSAPrint(1,
                                        "AFTER TRANSFORM Max " + SystemMax.ToString("E4") + " Average " +
                                        SystemAverage.ToString("E4") + " Sigma " + SystemSigma.ToString("E4") +
                                        " Estimated Dimension " + EstimatedDimension.ToString("F2") +
                                        "\n Center " + Center.ToString() + " With Link Cut " +
                                        LinkCutforCenter.ToString() + " Radius " + SystemRadius.ToString("E4") +
                                        " xAxis " + xAxis.ToString() + " " + xAxisExtent.ToString("E4") + " xyPlane " +
                                        xyPlane.ToString()
                                        + " " + xyPlaneExtent.ToString("E4") + "\n Minimum Distance " +
                                        MinimumDistance.ToString("E4") + " Distances Less than this " +
                                        DistancesNearEachOther.ToString("F0") + " Points Affected " +
                                        NotLonelyPoints.ToString());
                string histogramcounts = "";
                for (int binloop = 0; binloop < (2 + PointsinDistanceHistogram); binloop++)
                {
                    histogramcounts += Bincounts[binloop].ToString("F0") + ", ";
                }
                double BinSize = (Histmax - Histmin)/PointsinDistanceHistogram;
                SALSAUtility.SALSAPrint(1,
                                        "\nDistance Histogram Min " + Histmin.ToString("E4") + " Max " +
                                        Histmax.ToString("E4") + " Binsize " + BinSize.ToString("E4")
                                        + " #Counts with under/overflow\n" + histogramcounts);
            }
            return;
        }

        //  End SetupChisq

        public static void Sequel()
        {
            /* Sequel of ManxcatMDS work */

            SALSAUtility.SALSAPrint(1, "\nStarting density and histogram generation");

            var orignalPnumToCnumTable = new Hashtable();
            if (SALSAUtility.IsClustersSelected)
            {
                ReadClusterFile(orignalPnumToCnumTable);
            }

            double xminWhole = double.MaxValue;
            double xmaxWhole = double.MinValue;
            double yminWhole = double.MaxValue;
            double ymaxWhole = double.MinValue;

            double xminSelected = double.MaxValue;
            double xmaxSelected = double.MinValue;
            double yminSelected = double.MaxValue;
            double ymaxSelected = double.MinValue;

            double xminSelectedInter = double.MaxValue;
            double xmaxSelectedInter = double.MinValue;
            double yminSelectedInter = double.MaxValue;
            double ymaxSelectedInter = double.MinValue;

            SALSAUtility.SALSAPrint(1, "\n\tFinding min/max values");

            ManxcatMDSBasicDataProcessing.FindMinMaxForDensity(ref xminWhole, ref xmaxWhole, ref yminWhole,
                                                               ref ymaxWhole, ref xminSelected, ref xmaxSelected,
                                                               ref yminSelected, ref ymaxSelected, ref xminSelectedInter,
                                                               ref xmaxSelectedInter, ref yminSelectedInter,
                                                               ref ymaxSelectedInter,
                                                               orignalPnumToCnumTable);
            SALSAUtility.SALSAPrint(1, string.Format("\n\t\txminWhole: {0}\txmaxWhole: {1}", xminWhole, xmaxWhole));
            SALSAUtility.SALSAPrint(1, string.Format("\t\tyminWhole: {0}\tymaxWhole: {1}", yminWhole, ymaxWhole));
            SALSAUtility.SALSAPrint(1,
                                    string.Format("\n\t\txminSelected: {0}\txmaxSelected: {1}", xminSelected,
                                                  xmaxSelected));
            SALSAUtility.SALSAPrint(1,
                                    string.Format("\t\tyminSelected: {0}\tymaxSelected: {1}", yminSelected, ymaxSelected));
            SALSAUtility.SALSAPrint(1,
                                    string.Format("\n\t\txminSelectedInter: {0}\txmaxSelectedInter: {1}",
                                                  xminSelectedInter, xmaxSelectedInter));
            SALSAUtility.SALSAPrint(1,
                                    string.Format("\t\tyminSelectedInter: {0}\tymaxSelectedInter: {1}",
                                                  yminSelectedInter, ymaxSelectedInter));

            SALSAUtility.SALSAPrint(1, "\n\tDone");


            double deltaxWhole = (xmaxWhole - xminWhole)/SALSAUtility.Xres;
            double deltayWhole = (ymaxWhole - yminWhole)/SALSAUtility.Yres;
            double deltasWhole = deltaxWhole*deltayWhole;

            double deltaxSelected = (xmaxSelected - xminSelected)/SALSAUtility.Xres;
            double deltaySelected = (ymaxSelected - yminSelected)/SALSAUtility.Yres;
            double deltasSelected = deltaxSelected*deltaySelected;

            double deltaxSelectedInter = (xmaxSelectedInter - xminSelectedInter)/SALSAUtility.Xres;
            double deltaySelectedInter = (ymaxSelectedInter - yminSelectedInter)/SALSAUtility.Yres;
            double deltasSelectedInter = deltaxSelectedInter*deltaySelectedInter;

            double[][] densityMatrixWhole, densityMatrixSelected, densityMatrixSelectedInter;
            double[] xHistogramWhole, yHistogramWhole;
            double[] xHistogramSelected, yHistogramSelected;
            double[] xHistogramSelectedInter, yHistogramSelectedInter;

            double countWhole;
            double countSelected;
            double countSelectedInter;

            SALSAUtility.SALSAPrint(1, "\n\tGenerating density matrices and histograms");
            ManxcatMDSBasicDataProcessing.GenerateDensityMatrix(out densityMatrixWhole, out xHistogramWhole,
                                                                out yHistogramWhole,
                                                                out densityMatrixSelected, out xHistogramSelected,
                                                                out yHistogramSelected, out densityMatrixSelectedInter,
                                                                out xHistogramSelectedInter, out yHistogramSelectedInter,
                                                                xminWhole, xmaxWhole, yminWhole,
                                                                ymaxWhole, xminSelected, xmaxSelected, yminSelected,
                                                                ymaxSelected,
                                                                xminSelectedInter, xmaxSelectedInter, yminSelectedInter,
                                                                ymaxSelectedInter,
                                                                out countWhole, out countSelected,
                                                                out countSelectedInter,
                                                                orignalPnumToCnumTable);

            SALSAUtility.SALSAPrint(1, "\tDone");
            if (SALSAUtility.MPI_Rank == 0)
            {
                /* Density matrices and histograms should be good here */
                string outDir = Path.GetDirectoryName(ManxcatCentral.Configuration.ReducedVectorOutputFileName);
                SALSAUtility.SALSAPrint(1, "\n\tGenerating files for whole sample");
                GenerateDensityDataFile(densityMatrixWhole, xHistogramWhole, yHistogramWhole, xminWhole, xmaxWhole,
                                        yminWhole, ymaxWhole, deltaxWhole, deltayWhole, deltasWhole, countWhole,
                                        outDir, "whole");
                SALSAUtility.SALSAPrint(1, "\tDone");
                if (SALSAUtility.IsClustersSelected)
                {
                    SALSAUtility.SALSAPrint(1, "\n\tGenerating files for selected cluster (intra)");
                    GenerateDensityDataFile(densityMatrixSelected, xHistogramSelected, yHistogramSelected,
                                            xminSelected, xmaxSelected, yminSelected, ymaxSelected, deltaxSelected,
                                            deltaySelected, deltasSelected, countSelected, outDir, "selected");
                    SALSAUtility.SALSAPrint(1, "\tDone");

                    SALSAUtility.SALSAPrint(1, "\n\tGenerating files for selected cluster (inter)");
                    GenerateDensityDataFile(densityMatrixSelectedInter, xHistogramSelectedInter, yHistogramSelectedInter,
                                            xminSelectedInter, xmaxSelectedInter, yminSelectedInter, ymaxSelectedInter,
                                            deltaxSelectedInter, deltaySelectedInter, deltasSelectedInter,
                                            countSelectedInter, outDir, "selected-inter");
                    SALSAUtility.SALSAPrint(1, "\tDone");
                }
            }
        }

        private static void GenerateDensityDataFile(double[][] densityMatrx, double[] xHist, double[] yHist,
                                                    double xmin, double xmax, double ymin, double ymax,
                                                    double deltax, double deltay, double deltas, double count,
                                                    string outDir, string prefix)
        {
            double cellmax = 0.0;
            for (int i = 0; i < SALSAUtility.Yres; i++)
            {
                for (int j = 0; j < SALSAUtility.Xres; j++)
                {
                    if (densityMatrx[i][j] > cellmax)
                    {
                        cellmax = densityMatrx[i][j];
                    }
                }
            }
            double cellmean = count/(SALSAUtility.Xres*SALSAUtility.Yres);
            double power = cellmax < (SALSAUtility.Alpha*cellmean)
                               ? 1.0
                               : (Math.Log(SALSAUtility.Alpha)/Math.Log(cellmax/cellmean));
            // Constant value by which the number of points in a 2D square is multiplied.
            // The resulting value is independent of the total number of points as well as 
            // the x,y resolution. The mult value is a factor changing the z value scale.
            double c = 1.0/cellmax;

            // Output density values
            Console.WriteLine(new string('*', 40));
            Console.WriteLine("DataSet\t" + prefix);
            Console.WriteLine("Count\t" + count);
            Console.WriteLine("Deltas\t" + deltas);
            Console.WriteLine("CellMean\t" + cellmean);
            Console.WriteLine("CellMax\t" + cellmax);
            Console.WriteLine("Power\t" + power);
            Console.WriteLine("Const\t" + c);
            for (int i = 0; i < 10; i++)
            {
                double density = i/10.0;
                double densityToCount = Math.Pow(density, (1/power))/c;
                Console.WriteLine(density + "\t" + densityToCount);
            }
            Console.WriteLine(new string('*', 40));


            int xpointcount = 2*SALSAUtility.Xres;
            int ypointcount = 2*SALSAUtility.Yres;

            string densityFile = Path.Combine(outDir, prefix + "-density.txt");

            string xHistFile = Path.Combine(outDir, prefix + "-xHist.txt");

            string yHistFile = Path.Combine(outDir, prefix + "-yHist.txt");

            string scriptFile = Path.Combine(outDir, prefix + "-plot.txt");

            using (StreamWriter writer = new StreamWriter(densityFile),
                                xHistWriter = new StreamWriter(xHistFile),
                                yHistWriter = new StreamWriter(yHistFile),
                                scriptWriter = new StreamWriter(scriptFile))
            {
                writer.WriteLine("#xcoord\tycoord\thistogramValue");
                xHistWriter.WriteLine("#xval\thistogramvalue");
                yHistWriter.WriteLine("#yval\thistogramvalue");

                // Generating x histogram
                double xoffset = xmin + 0.5*deltax;
                for (int i = 0; i < SALSAUtility.Xres; ++i)
                {
                    double xcoord = xoffset + i*deltax;
                    if (SALSAUtility.Normalize) xcoord /= xmax;
                    xHistWriter.WriteLine(xcoord + "\t" + xHist[i]);
                }

                // Generating y histogram
                double yoffset = ymin + 0.5*deltay;
                for (int i = 0; i < SALSAUtility.Yres; ++i)
                {
                    double ycoord = yoffset + i*deltay;
                    if (SALSAUtility.Normalize) ycoord /= ymax;
                    yHistWriter.WriteLine(ycoord + "\t" + yHist[i]);
                }

                for (int i = 0; i < xpointcount; i++)
                {
                    double x = xmin + ((IsOdd(i) ? (i + 1)/2 : i/2)*deltax);
                    int cellx = IsOdd(i) ? (i - 1)/2 : i/2;

                    for (int j = 0; j < ypointcount; j++)
                    {
                        double y = ymin + ((IsOdd(j) ? (j + 1)/2 : j/2)*deltay);
                        int celly = IsOdd(j) ? (j - 1)/2 : j/2;

                        double cellvalue = Math.Pow((densityMatrx[celly][cellx]*c), power);

                        cellvalue = cellvalue > SALSAUtility.Pcutf ? SALSAUtility.Pcutf : cellvalue;

                        if (SALSAUtility.Normalize)
                        {
                            writer.WriteLine(x/xmax + "\t" + y/ymax + "\t" + cellvalue);
                        }
                        else
                        {
                            writer.WriteLine(x + "\t" + y + "\t" + cellvalue);
                        }
                    }
                    writer.WriteLine();
                }

                // Fill up the remaining region from beyond x=xmax and y=ymax as zero 
                writer.WriteLine();
                if (SALSAUtility.Normalize)
                {
                    writer.WriteLine(xmin/xmax + "\t" + ymax/ymax + "\t" + 0.0);
                    writer.WriteLine(xmin/xmax + "\t" + SALSAUtility.Xmaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(xmax/xmax + "\t" + ymax/ymax + "\t" + 0.0);
                    writer.WriteLine(xmax/xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(xmax/xmax + "\t" + ymin/ymax + "\t" + 0.0);
                    writer.WriteLine(xmax/xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(SALSAUtility.Xmaxbound + "\t" + ymin/ymax + "\t" + 0.0);
                    writer.WriteLine(SALSAUtility.Xmaxbound + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                }
                else
                {
                    writer.WriteLine(xmin + "\t" + ymax + "\t" + 0.0);
                    writer.WriteLine(xmin + "\t" + SALSAUtility.Xmaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(xmax + "\t" + ymax + "\t" + 0.0);
                    writer.WriteLine(xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(xmax + "\t" + ymin + "\t" + 0.0);
                    writer.WriteLine(xmax + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                    writer.WriteLine();
                    writer.WriteLine(SALSAUtility.Xmaxbound + "\t" + ymin + "\t" + 0.0);
                    writer.WriteLine(SALSAUtility.Xmaxbound + "\t" + SALSAUtility.Ymaxbound + "\t" + 0.0);
                }


                WriteGnuplotScript(prefix, densityFile, xHistFile, yHistFile, scriptWriter);
            }
        }

        private static void WriteGnuplotScript(string prefix, string densityFile, string xHistFile, string yHistFile,
                                               StreamWriter scriptWriter)
        {
            scriptWriter.WriteLine("set terminal png truecolor nocrop font arial 14 size 1200,1200 xffffff");
            scriptWriter.WriteLine();

            string pngfile = prefix + (SALSAUtility.Normalize ? "-normalized" : string.Empty) + "-plot.png";
            scriptWriter.WriteLine("set output '" + pngfile + "'");

            scriptWriter.WriteLine("set size 1.0, 1.0");
            scriptWriter.WriteLine("set multiplot");

            scriptWriter.WriteLine();

            // Title box
            scriptWriter.WriteLine("set origin 0.0, 0.85");
            scriptWriter.WriteLine("set size 0.95, 0.1");
            scriptWriter.WriteLine("set border linecolor rgbcolor \"white\"");
            scriptWriter.WriteLine("unset key");
            string title = SALSAUtility.ManxcatRunName + " (" + prefix +
                           (SALSAUtility.Normalize ? "-normalized" : string.Empty) + ")";
            scriptWriter.WriteLine("set title \"" + title + "\" textcolor rgbcolor \"black\"");
            scriptWriter.WriteLine("plot [0:1] [0:1] 0.0 lt rgb \"white\"");

            scriptWriter.WriteLine("set border linecolor rgbcolor \"black\"");
            scriptWriter.WriteLine("set dummy u,v");
            scriptWriter.WriteLine("unset key");
            scriptWriter.WriteLine("set size ratio 1.0");
            scriptWriter.WriteLine("set style fill  solid 0.85 noborder");
            scriptWriter.WriteLine("set style line 1 lt 1 lw 4");
            scriptWriter.WriteLine("set pm3d map");
            scriptWriter.WriteLine("set palette rgbformulae 30,31,32 model RGB negative");

            scriptWriter.WriteLine();

            // Y histogram (rotated)
            scriptWriter.WriteLine("set origin 0.0, 0.45");
            scriptWriter.WriteLine("set size 0.45, 0.45");
            scriptWriter.WriteLine("set xtics rotate by -90");
            string xlabel = "Count";
            scriptWriter.WriteLine("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
            string ylabel = "Euclidean Distance";
            scriptWriter.WriteLine("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
            title = "Histogram (rotated) of " + ylabel;
            scriptWriter.WriteLine("set title \"" + title + "\" textcolor rgbcolor \"black\"");
            scriptWriter.WriteLine("plot [][:" + SALSAUtility.Ymaxbound + "] '" + Path.GetFileName(yHistFile) +
                                   "' using 2:1 with filledcurves y1 lt rgb \"black\"");

            scriptWriter.WriteLine("set xtics rotate by 0");
            scriptWriter.WriteLine();

            // Density plot
            scriptWriter.WriteLine("set origin 0.45, 0.45");
            scriptWriter.WriteLine("set size 0.5, 0.5");
            xlabel = "Original Distance";
            scriptWriter.WriteLine("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
            ylabel = "Euclidena Distance";
            scriptWriter.WriteLine("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
            title = "Heat Map of " + ylabel + " vs " + xlabel;
            scriptWriter.WriteLine("set title \"" + title + "\" textcolor rgbcolor \"black\"");
            scriptWriter.WriteLine("splot [:" + SALSAUtility.Xmaxbound + "] [:" + SALSAUtility.Ymaxbound + "] '" +
                                   Path.GetFileName(densityFile) + "'");


            scriptWriter.WriteLine();

            // Y histogram (unrotated)
            scriptWriter.WriteLine("set origin 0.0, 0.0");
            scriptWriter.WriteLine("set size 0.45, 0.45");
            xlabel = "Euclidean Distance";
            scriptWriter.WriteLine("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
            ylabel = "Count";
            scriptWriter.WriteLine("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
            title = "Histogram of " + xlabel;
            scriptWriter.WriteLine("set title \"" + title + "\" textcolor rgbcolor \"black\"");
            scriptWriter.WriteLine("plot [:" + SALSAUtility.Ymaxbound + "] []'" + Path.GetFileName(yHistFile) +
                                   "' with filledcurves x1 lt rgb \"black\"");


            scriptWriter.WriteLine();

            // X histogram
            scriptWriter.WriteLine("set origin 0.45, 0.0");
            scriptWriter.WriteLine("set size 0.45, 0.45");
            xlabel = "Original Distance";
            scriptWriter.WriteLine("set xlabel \"" + xlabel + "\" textcolor rgbcolor \"black\"");
            ylabel = "Count";
            scriptWriter.WriteLine("set ylabel \"" + ylabel + "\" textcolor rgbcolor \"black\"");
            title = "Histogram of " + xlabel;
            scriptWriter.WriteLine("set title \"" + title + "\" textcolor rgbcolor \"black\"");
            scriptWriter.WriteLine("plot [:" + SALSAUtility.Xmaxbound + "] []'" + Path.GetFileName(xHistFile) +
                                   "' with filledcurves x1 lt rgb \"black\"");

            scriptWriter.WriteLine();

            scriptWriter.WriteLine("unset multiplot");
        }

        private static bool IsOdd(int value)
        {
            return (value & 1) == 1;
        }

        private static void ReadClusterFile(Hashtable orignalPnumToCnumTable)
        {
            using (var reader = new StreamReader(SALSAUtility.ClusterFile))
            {
                var sep = new[] {'\t', ' '};
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    if (!string.IsNullOrEmpty(line))
                    {
                        string[] splits = line.Trim().Split(sep);
                        int pnum = int.Parse(splits[0]);
                        int cnum = int.Parse(splits[1]);
                        if (orignalPnumToCnumTable.Contains(pnum))
                        {
                            throw new Exception("Point numbers in the cluster file should be unique");
                        }
                        orignalPnumToCnumTable.Add(pnum, cnum);
                    }
                }
            }
        }

        public static bool Calcfg(Desertwind Solution)
        {
            // Assume zerocr and first/DiagonalofMatrix are zeroed before call
            // Here we only calculate diagonal elements of Chisq matrix

            bool violat = false;
            Hotsun.succ = true;
            Hotsun.zerocr = 0.0;
            Hotsun.idata = 0; // Not used in MDS

            var Findzerocr = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 double WeightFunction1, WeightFunction2;
                                 double localzerocr = 0.0;
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LocalToProcessIndex1 = beginpoint;
                                      LocalToProcessIndex1 < indexlen + beginpoint;
                                      LocalToProcessIndex1++)
                                 {
                                     int GlobalIndex1 = LocalToProcessIndex1 + SALSAUtility.PointStart_Process;
                                     int OriginalPointIndex1 = SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex1];
                                     bool vary1 = true;
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         vary1 = false;

                                     // Skip Fixed rows if StoredDistanceOption =3
                                     if ((SALSAUtility.StoredDistanceOption == 3) && !vary1)
                                         continue;

                                     for (int GlobalIndex2 = 0;
                                          GlobalIndex2 < SALSAUtility.PointCount_Global;
                                          GlobalIndex2++)
                                     {
                                         if (GlobalIndex1 == GlobalIndex2)
                                             continue;

                                         int OriginalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex2];
                                         bool vary2 = true;
                                         if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             vary2 = false;

                                         // Varied (row) - Fixed (column) Contribution to Chisq Doubled if StoredDistanceOption 3 and so no row for fixed parameters
                                         // vary1 forced to be true in this scenario
                                         bool symmetrize = false;
                                         if ((SALSAUtility.StoredDistanceOption == 3) && !vary2)
                                             symmetrize = true;

                                         // If requested skip chisq for fixed cross fixed points
                                         if (!SALSAUtility.CalcFixedCrossFixed)
                                         {
                                             if (!vary1 && !vary2)
                                                 continue;
                                         }

                                         if (
                                             !SetChisqWeights(GlobalIndex1, GlobalIndex2, out WeightFunction1,
                                                              out WeightFunction2))
                                             continue;

                                         // Calculate contribution to Chisq
                                         double SquaredDistance = 0.0;

                                         for (int LocalVectorIndex1 = 0;
                                              LocalVectorIndex1 < Hotsun.ParameterVectorDimension;
                                              LocalVectorIndex1++)
                                         {
                                             double tmp = Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex1] -
                                                          Hotsun.GlobalParameter[GlobalIndex2][LocalVectorIndex1];
                                             SquaredDistance += tmp*tmp;
                                         }
                                         double DistanceFudge1 = 1.0;
                                         double DistanceFudge2 = 1.0;
                                         double ActualDistance = SquaredDistance;

                                         if (DistanceFormula == 1)
                                         {
                                             SquaredDistance += SQUAREMinimumDistance;
                                             ActualDistance = Math.Sqrt(SquaredDistance);
                                             DistanceFudge1 = 0.5/ActualDistance;
                                             DistanceFudge2 = DistanceFudge1/SquaredDistance;
                                         }

                                         double funcl = WeightFunction1 - ActualDistance*WeightFunction2;
                                         double symmetryweight = 1.0;
                                         if (symmetrize)
                                             symmetryweight = 2.0;
                                         double increment = funcl*funcl*symmetryweight;
                                         if (SALSAUtility.chisqcomponent > 0)
                                         {
                                             if ((SALSAUtility.chisqcomponent == 1) && (!vary1 || !vary2))
                                                 // Varied Varied
                                                 continue;
                                             if ((SALSAUtility.chisqcomponent == 2) && (!vary1 || vary2))
                                                 // Varied Fixed
                                                 continue;
                                             if ((SALSAUtility.chisqcomponent == 3) && (vary1 || !vary2))
                                                 // Fixed Varied
                                                 continue;
                                             if ((SALSAUtility.chisqcomponent == 4) && (vary1 || vary2)) // Fixed Fixed
                                                 continue;
                                         }
                                         localzerocr += increment;

                                         if (SALSAUtility.chisqcomponent > 0)
                                             continue;
                                         if (!vary1)
                                             continue;
                                         for (int LocalVectorIndex1 = 0;
                                              LocalVectorIndex1 < Hotsun.ParameterVectorDimension;
                                              LocalVectorIndex1++)
                                         {
                                             double tmp1 = (Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex1] -
                                                            Hotsun.GlobalParameter[GlobalIndex2][LocalVectorIndex1]);

                                             //  Calculate First Derivative
                                             if (!Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1])
                                             {
                                                 Solution.first[LocalToProcessIndex1][LocalVectorIndex1] += -4.0*tmp1*
                                                                                                            funcl*
                                                                                                            DistanceFudge1*
                                                                                                            WeightFunction2;
                                             }

                                             // Calculate Diagonal Matrix Contributions to Second Derivative
                                             for (int LocalVectorIndex2 = 0;
                                                  LocalVectorIndex2 < Hotsun.ParameterVectorDimension;
                                                  LocalVectorIndex2++)
                                             {
                                                 double tmp2 =
                                                     (Hotsun.GlobalParameter[GlobalIndex1][LocalVectorIndex2] -
                                                      Hotsun.GlobalParameter[GlobalIndex2][LocalVectorIndex2]);
                                                 double ApproximateCalc = 8.0*tmp1*tmp2*DistanceFudge1*DistanceFudge1*
                                                                          WeightFunction2*WeightFunction2;
                                                 Solution.DiagonalofMatrix[LocalToProcessIndex1][
                                                     LocalVectorIndex1, LocalVectorIndex2] += ApproximateCalc;

                                                 double correction = 0.0;

                                                 if ((DistanceFormula == 2) && (LocalVectorIndex1 == LocalVectorIndex2))
                                                     correction = -4.0*funcl*WeightFunction2;

                                                 if ((DistanceFormula == 1) && (LocalVectorIndex1 == LocalVectorIndex2))
                                                     correction = -4.0*funcl*WeightFunction2*DistanceFudge1;

                                                 if (DistanceFormula == 1)
                                                     correction += 4.0*funcl*tmp1*tmp2*DistanceFudge2*WeightFunction2;
                                                 Solution.ExactDiagonalofMatrix[LocalToProcessIndex1][
                                                     LocalVectorIndex1, LocalVectorIndex2] += ApproximateCalc +
                                                                                              correction;
                                             } // // End LocalVectorIndex2
                                         } // End LocalVectorIndex1
                                     } // End GlobalPointIndex2
                                 } // End LocaltoProcessIndex1

                                 Findzerocr.addapoint(ThreadNo, localzerocr);
                             }); // End loop over Point dependent quantities

            Findzerocr.sumoverthreadsandmpi();
            Hotsun.zerocr = Findzerocr.Total;
            Solution.Chisquared = Hotsun.zerocr;
            return violat;
        }

        // End Calcfg

        //  DistanceProcessingOption = 1 and DistanceFormula = 1
        //  Chisqnorm = 0 ** (Euclidean Distance - Observed Distance) / Observed Distance
        //  Chisqnorm = 1 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.75
        //  Chisqnorm = 2 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.5     Sammon
        //  Chisqnorm = 3 ** (Euclidean Distance - Observed Distance)                               SMACOF
        //  Chisqnorm = 4 ** (Euclidean Distance - Observed Distance) * (Observed Distance)^0.5
        //
        //  DistanceProcessingOption = 2 and DistanceFormula = 2
        //  Chisqnorm = 0 ** (Euclidean Distance - Observed Distance) / Observed Distance
        //  Chisqnorm = 1 ** (Euclidean Distance - Observed Distance) / (Observed Distance)^0.5     Sammon
        //  Chisqnorm = 2 ** (Euclidean Distance - Observed Distance)                               SMACOF
        //  Chisqnorm = 3 ** (Euclidean Distance - Observed Distance) * Observed Distance
        //  Chisqnorm = 4 ** (Euclidean Distance - Observed Distance) * (Observed Distance)^2
        //
        //  DistanceProcessingOption = 1 and DistanceFormula = 2
        //  Chisqnorm = 0 ** (Euclidean Distance - SQRT(Observed Distance)) /  Observed Distance^0.5
        //  Chisqnorm = 1 ** (Euclidean Distance - SQRT(Observed Distance)) / (Observed Distance)^0.25
        //  Chisqnorm = 2 ** (Euclidean Distance - SQRT(Observed Distance))
        //  Chisqnorm = 3 ** (Euclidean Distance - SQRT(Observed Distance)) * (Observed Distance)^0.5
        //  Chisqnorm = 4 ** (Euclidean Distance - SQRT(Observed Distance)) *  Observed Distance

        //  Return false if interpoint distance not set
        public static bool SetChisqWeights(int GlobalPointIndex1, int GlobalPointIndex2, out double WeightFunction1,
                                           out double WeightFunction2)
        {
            double InterPointDistance = SALSAParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
            if (InterPointDistance < -0.5)
            {
                WeightFunction1 = 0.0;
                WeightFunction2 = 0.0;
                return false;
            }

            double Weighting = 1.0;
            if (SALSAUtility.NumberDistanceWeightCuts > 0)
            {
                int distpos = SALSAUtility.NumberDistanceWeightCuts;
                for (int weightbin = 0; weightbin < SALSAUtility.NumberDistanceWeightCuts; weightbin++)
                {
                    if (SALSAUtility.DistanceWeightCuts[weightbin] > InterPointDistance)
                    {
                        distpos = weightbin;
                        break;
                    }
                }
                Weighting = SALSAUtility.ActualWeightCuts[distpos];
            }

            if (WeightingOption > 0)
                Weighting = Math.Sqrt(PointWeights[GlobalPointIndex1]*PointWeights[GlobalPointIndex2]);

            if (Chisqnorm == 0)
            {
                if (InterPointDistance > MinimumDistance)
                {
                    WeightFunction1 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier;
                    WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/InterPointDistance;
                    return true;
                }
                WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/MinimumDistance;
                WeightFunction1 = WeightFunction2*InterPointDistance;
                return true;
            }

            if (Chisqnorm == 1)
            {
                // Sammon if use Distance Squared with DistanceFormula == 2
                if (InterPointDistance > MinimumDistance)
                {
                    double SQRTDistance = Math.Sqrt(InterPointDistance);
                    WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/
                                      (SQRTDistance*Math.Sqrt(SQRTSystemAverage*SQRTDistance));
                    WeightFunction1 = WeightFunction2*InterPointDistance;
                    return true;
                }
                WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/
                                  (SQRTMinimumDistance*Math.Sqrt(SQRTSystemAverage*SQRTMinimumDistance));
                WeightFunction1 = WeightFunction2*InterPointDistance;
                return true;
            }

            if (Chisqnorm == 2)
            {
                // SMACOF if use Distance Squared  with DistanceFormula == 2
                // Sammon if use NonSquared Distance with DistanceFormula == 1

                if (InterPointDistance > MinimumDistance)
                {
                    double SQRTDistance = Math.Sqrt(InterPointDistance);
                    WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/
                                      (SQRTDistance*SQRTSystemAverage);
                    WeightFunction1 = WeightFunction2*InterPointDistance;
                    return true;
                }
                WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/
                                  (SQRTMinimumDistance*SQRTSystemAverage);
                WeightFunction1 = WeightFunction2*InterPointDistance;
                return true;
            }

            if (Chisqnorm == 3)
            {
                // SMACOF if use Nonsquared Distance   with DistanceFormula == 1
                WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier/SystemAverage;
                WeightFunction1 = WeightFunction2*InterPointDistance;
                return true;
            }

            if (Chisqnorm == 4)
            {
                // Even more large distance weight than SMACOF if DistanceFormula == 1
                WeightFunction2 = Weighting*ManxcatCentral.ChisqFunctionCalcMultiplier*Math.Sqrt(InterPointDistance)/
                                  (SQRTSystemAverage*SystemAverage);
                WeightFunction1 = WeightFunction2*InterPointDistance;
                return true;
            }

            //  Default
            WeightFunction1 = Weighting;
            WeightFunction2 = Weighting;
            return true;
        }

        // End SetChisqWeights

        //  Set initial value of param in Solution
        //  InitializationOption =0 None (except fixed points) =1 Old File Type Model =2 New File Type
        //  InitializationOption 1 must have a file with length identical to number of used points
        //  For fixed points, any values read in with InitializationOption 1 or 2 overwrite those read in SALSAProcessVariedandFixed
        //  Such overwriting must happen for InitializationOption 1 but only happens for InitializationOption 2 if fixed points read
        public static void InitializeParameters(Desertwind Solution, int CountStartingPoints)
        {
            double RadiusUsed = Math.Sqrt(SystemRadius);
            var Randobject = new Random();
            var InitialMDSString = new string[SALSAUtility.PointCount_Global];
            string OriginalMDSFileName = ManxcatCentral.Configuration.InitializationFileName;

            if (!OriginalMDSFileName.Contains(":") && !OriginalMDSFileName.Contains("$"))
                OriginalMDSFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + OriginalMDSFileName;
            int InitializationNumberofPoints = 0;

            if (ManxcatCentral.Configuration.InitializationOption == 1)
            {
                // This file must be same length as number of used points

                var ColorValue = new int[SALSAUtility.PointCount_Global];
                if (ManxcatMDSDataProcessing.ReadMDSCluster_File(OriginalMDSFileName, InitialMDSString, ColorValue,
                                                                 ref InitializationNumberofPoints))
                    SALSAUtility.SALSAPrint(0,
                                            OriginalMDSFileName + " Read Successfully Points " +
                                            InitializationNumberofPoints.ToString());
                if (SALSAUtility.PointCount_Global != InitializationNumberofPoints)
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Inconsistent Initialization File Point Counts " +
                                                InitializationNumberofPoints.ToString() +
                                                " Expected is " + SALSAUtility.PointCount_Global.ToString());
                    throw (e);
                }
            }

            if (ManxcatCentral.Configuration.InitializationOption == 2)
            {
                int InitializationFileType = -1;
                var InitializationFileProperties = new SALSAFileProperties();
                var InitializationPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];
                SALSA_Properties.ReadDataPointFile(OriginalMDSFileName, ref InitializationFileType,
                                                   InitializationFileProperties, ref InitializationPointProperties,
                                                   ref InitializationNumberofPoints);

                if ((SALSAUtility.NumberOriginalPoints < InitializationNumberofPoints) ||
                    (InitializationFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints))
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Inconsistent Initialization File Point Counts " +
                                                InitializationNumberofPoints.ToString() + " or "
                                                + InitializationFileProperties.NumberOriginalPoints.ToString() +
                                                " Expected is " + SALSAUtility.NumberOriginalPoints.ToString());
                    throw (e);
                }

                for (int InitialIndex = 0; InitialIndex < InitializationNumberofPoints; InitialIndex++)
                {
                    int OriginalIndex = InitializationPointProperties[InitialIndex].OriginalPointNumber;
                    if (!InitializationPointProperties[InitialIndex].valuesset)
                        continue;

                    // Consider all used points
                    int UsedIndex = SALSAUtility.OriginalPointtoUsedPointMap[OriginalIndex];
                    if (UsedIndex < 0) continue;

                    SALSAUtility.GlobalPointProperties[UsedIndex].x = InitializationPointProperties[InitialIndex].x;
                    SALSAUtility.GlobalPointProperties[UsedIndex].y = InitializationPointProperties[InitialIndex].y;
                    SALSAUtility.GlobalPointProperties[UsedIndex].z = InitializationPointProperties[InitialIndex].z;
                    SALSAUtility.GlobalPointProperties[UsedIndex].valuesset = true;
                    SALSAUtility.GlobalPointProperties[UsedIndex].source =
                        InitializationPointProperties[InitialIndex].source;

                    if (InitializationPointProperties[InitialIndex].errorsset)
                    {
                        SALSAUtility.GlobalPointProperties[UsedIndex].xerr =
                            InitializationPointProperties[InitialIndex].xerr;
                        SALSAUtility.GlobalPointProperties[UsedIndex].yerr =
                            InitializationPointProperties[InitialIndex].yerr;
                        SALSAUtility.GlobalPointProperties[UsedIndex].zerr =
                            InitializationPointProperties[InitialIndex].zerr;
                        SALSAUtility.GlobalPointProperties[UsedIndex].errorsset = true;
                    }
                    SALSAUtility.GlobalPointProperties[UsedIndex].family1 =
                        InitializationPointProperties[InitialIndex].family1;
                    SALSAUtility.GlobalPointProperties[UsedIndex].family2 =
                        InitializationPointProperties[InitialIndex].family2;
                    SALSAUtility.GlobalPointProperties[UsedIndex].cluster =
                        InitializationPointProperties[InitialIndex].cluster;
                    SALSAUtility.GlobalPointProperties[UsedIndex].familylabel1 =
                        InitializationPointProperties[InitialIndex].familylabel1;
                    SALSAUtility.GlobalPointProperties[UsedIndex].familylabel2 =
                        InitializationPointProperties[InitialIndex].familylabel2;
                    SALSAUtility.GlobalPointProperties[UsedIndex].clusterlabel =
                        InitializationPointProperties[InitialIndex].clusterlabel;
                    SALSAUtility.GlobalPointProperties[UsedIndex].pointlabel =
                        InitializationPointProperties[InitialIndex].pointlabel;
                    SALSAUtility.GlobalPointProperties[UsedIndex].PointType =
                        InitializationPointProperties[InitialIndex].PointType;
                }
            }

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LongIndex = beginpoint; LongIndex < indexlen + beginpoint; LongIndex++)
                                 {
                                     int GlobalIndex = LongIndex + SALSAUtility.PointStart_Process;

                                     if (ManxcatCentral.Configuration.InitializationOption == 1)
                                     {
                                         for (int LocalVectorIndex = 0;
                                              LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                              LocalVectorIndex++)
                                         {
                                             string[] split =
                                                 InitialMDSString[SALSAUtility.ActualtoNaiveUsedOrder[GlobalIndex]].
                                                     Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                                             Solution.param[LongIndex][LocalVectorIndex] =
                                                 Convert.ToDouble(split[LocalVectorIndex + 1]);
                                         }
                                     }

                                     if (ManxcatCentral.Configuration.InitializationOption == 2)
                                     {
                                         if (!SALSAUtility.GlobalPointProperties[GlobalIndex].valuesset)
                                         {
                                             Exception e =
                                                 SALSAUtility.SALSAError("Error in Initialized Points -- Unset point " +
                                                                         GlobalIndex.ToString());
                                         }
                                         Solution.param[LongIndex][0] =
                                             SALSAUtility.GlobalPointProperties[GlobalIndex].x;
                                         Solution.param[LongIndex][1] =
                                             SALSAUtility.GlobalPointProperties[GlobalIndex].y;

                                         if (Hotsun.ParameterVectorDimension > 2)
                                             Solution.param[LongIndex][2] =
                                                 SALSAUtility.GlobalPointProperties[GlobalIndex].z;
                                     }

                                     if (ManxcatCentral.Configuration.InitializationOption == 0)
                                     {
                                         double bigorsmall = 2.0*Randobject.NextDouble();

                                         //                    if (bigorsmall < 0.5)
                                         //                        bigorsmall = 0.25;
                                         //                    else
                                         //                        bigorsmall = 4.0;

                                         //                     bigorsmall = 1.0; // NEW

                                         for (int LocalVectorIndex = 0;
                                              LocalVectorIndex < Hotsun.ParameterVectorDimension;
                                              LocalVectorIndex++)
                                         {
                                             if (Hotsun.FixedParameter[GlobalIndex][LocalVectorIndex])
                                             {
                                                 if (FindPointstoFix)
                                                 {
                                                     // Built in fixing of Points to remove rotation ambiguity
                                                     Solution.param[LongIndex][LocalVectorIndex] = 0.0;
                                                 }
                                                 else
                                                 {
                                                     // Explicit fixing of points in fit. These were set in GlobalPointProperties in SALSAProcessVariedandFixed
                                                     if (!SALSAUtility.GlobalPointProperties[GlobalIndex].valuesset)
                                                     {
                                                         Exception e =
                                                             SALSAUtility.SALSAError(
                                                                 "Error in Fixed Points -- Unset point " +
                                                                 GlobalIndex.ToString());
                                                     }

                                                     if (LocalVectorIndex == 0)
                                                         Solution.param[LongIndex][LocalVectorIndex] =
                                                             SALSAUtility.GlobalPointProperties[GlobalIndex].x;

                                                     if (LocalVectorIndex == 1)
                                                         Solution.param[LongIndex][LocalVectorIndex] =
                                                             SALSAUtility.GlobalPointProperties[GlobalIndex].y;

                                                     if (LocalVectorIndex == 2)
                                                         Solution.param[LongIndex][LocalVectorIndex] =
                                                             SALSAUtility.GlobalPointProperties[GlobalIndex].z;
                                                 }
                                             }
                                             else
                                             {
                                                 // Point is NOT fixed
                                                 Solution.param[LongIndex][LocalVectorIndex] = bigorsmall*
                                                                                               Randobject.NextDouble()*
                                                                                               RadiusUsed;
                                             }
                                         }
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // End InitializeParameters(Desertwind Solution)

        public static void SetupWeightings(double[] WeightsasRead)
        {
            if (WeightingOption == 0)
            {
                for (int Globalindex = 0; Globalindex < SALSAUtility.PointCount_Global; Globalindex++)
                {
                    WeightsasRead[Globalindex] = 1.0;
                }
                return;
            }

            string WeightFileName = ManxcatCentral.Configuration.WeightingFileName;

            if (!WeightFileName.Contains(":") && !WeightFileName.Contains("$"))
                WeightFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + WeightFileName;

            double sumofweights = 0.0;
            int NumberofLines = 0;
            int NumberUsedPoints = 0;

            try
            {
                // Check if file exists
                if (!File.Exists(WeightFileName))
                {
                    Exception e = SALSAUtility.SALSAError("File " + WeightFileName + " does not exists.");

                    throw (e);
                }

                // Create a new stream to read from a file
                using (StreamReader sr = File.OpenText(WeightFileName))
                {
                    // Read contents of a file, line by line, into a string
                    String inputLineStr;

                    while ((inputLineStr = sr.ReadLine()) != null)
                    {
                        inputLineStr = inputLineStr.Trim();

                        if (inputLineStr.Length < 1)
                            continue; //replace empty line

                        int fred = SALSAUtility.OriginalPointtoUsedPointMap[NumberofLines];
                        if (fred >= 0)
                        {
                            WeightsasRead[fred] = Convert.ToDouble(inputLineStr);
                            sumofweights += WeightsasRead[NumberofLines];
                            ++NumberUsedPoints;
                        }
                        ++NumberofLines;
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed to read data from " + WeightFileName + " " + e);

                throw (e);
            }

            if (NumberUsedPoints != SALSAUtility.PointCount_Global)
            {
                SALSAUtility.SALSAError("Incorrect Weight count read " + NumberUsedPoints + " Expected " +
                                        SALSAUtility.PointCount_Global);
            }
            double AveragetoOne = NumberUsedPoints/sumofweights;
            double minweight = sumofweights;
            double maxweight = 0.0;

            for (int Globalindex = 0; Globalindex < SALSAUtility.PointCount_Global; Globalindex++)
            {
                WeightsasRead[Globalindex] *= AveragetoOne;
                minweight = Math.Min(minweight, WeightsasRead[Globalindex]);
                maxweight = Math.Max(maxweight, WeightsasRead[Globalindex]);
            }
            SALSAUtility.SALSAPrint(1,
                                    "File " + WeightFileName + " Non trivial Point Weights Maximum " +
                                    maxweight.ToString("F3") + " Minimum " + minweight.ToString("E4"));
        }

        // End SetupWeightings()

        //  Calculate Matrix Global Vector product storing as a distributed vector
        public static void GlobalMatrixVectorProduct(double[][] DistributedVector, Desertwind Solution, bool useexact,
                                                     double[][] GlobalxVector, double[][] GlobalVectoronRight)
        {
            double[][,] MatrixDiagonals;

            if (useexact)
                MatrixDiagonals = Solution.ExactDiagonalofMatrix;
            else
                MatrixDiagonals = Solution.DiagonalofMatrix;

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 double WeightFunction1, WeightFunction2;
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int LocalProcessIndex1 = beginpoint;
                                      LocalProcessIndex1 < indexlen + beginpoint;
                                      LocalProcessIndex1++)
                                 {
                                     int GlobalIndex1 = LocalProcessIndex1 + SALSAUtility.PointStart_Process;

                                     for (int LocalVectorIndex1 = 0;
                                          LocalVectorIndex1 < Hotsun.ParameterVectorDimension;
                                          LocalVectorIndex1++)
                                     {
                                         if (Hotsun.FixedParameter[GlobalIndex1][LocalVectorIndex1])
                                             DistributedVector[LocalProcessIndex1][LocalVectorIndex1] = 0.0;
                                         else
                                         {
                                             double tmp = 0.0;

                                             for (int GlobalIndex2 = 0;
                                                  GlobalIndex2 < SALSAUtility.PointCount_Global;
                                                  GlobalIndex2++)
                                             {
                                                 int OriginalPointIndex2 =
                                                     SALSAUtility.UsedPointtoOriginalPointMap[GlobalIndex2];
                                                 if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                                     SALSAUtility.SALSASHIFT)
                                                     continue;

                                                 if (
                                                     !SetChisqWeights(GlobalIndex1, GlobalIndex2, out WeightFunction1,
                                                                      out WeightFunction2))
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

                                                     if (DistanceFormula == 1)
                                                     {
                                                         SquaredDistance += SQUAREMinimumDistance;
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
                                                             MatrixDiagonals[LocalProcessIndex1][
                                                                 LocalVectorIndex1, LocalVectorIndex2];

                                                         if (Hotsun.UseDiagonalScaling)
                                                             MatrixElement *=
                                                                 Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1]*
                                                                 Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex2];

                                                         if (Hotsun.AddMarquardtQDynamically &&
                                                             (LocalVectorIndex1 == LocalVectorIndex2))
                                                             MatrixElement += Hotsun.Q;
                                                     }
                                                     else
                                                     {
                                                         // Off Diagonal Term
                                                         double correction = 0.0;
                                                         double VectorCrossProduct =
                                                             (GlobalxVector[GlobalIndex1][LocalVectorIndex1] -
                                                              GlobalxVector[GlobalIndex2][LocalVectorIndex1])*
                                                             (GlobalxVector[GlobalIndex1][LocalVectorIndex2] -
                                                              GlobalxVector[GlobalIndex2][LocalVectorIndex2]);
                                                         MatrixElement = -8.0*VectorCrossProduct*DistanceFudge1*
                                                                         DistanceFudge1*WeightFunction2*WeightFunction2;

                                                         if (Hotsun.FullSecondDerivative)
                                                         {
                                                             if ((DistanceFormula == 2) &&
                                                                 (LocalVectorIndex1 == LocalVectorIndex2))
                                                                 correction = 4.0*funcl*WeightFunction2;

                                                             if ((DistanceFormula == 1) &&
                                                                 (LocalVectorIndex1 == LocalVectorIndex2))
                                                                 correction = 4.0*funcl*WeightFunction2*DistanceFudge1;

                                                             if (DistanceFormula == 1)
                                                                 correction += -4.0*funcl*VectorCrossProduct*
                                                                               WeightFunction2*DistanceFudge2;
                                                             MatrixElement += correction;
                                                         }

                                                         if (Hotsun.UseDiagonalScaling)
                                                             MatrixElement *=
                                                                 Hotsun.sqdginv[GlobalIndex1][LocalVectorIndex1]*
                                                                 Hotsun.sqdginv[GlobalIndex2][LocalVectorIndex2];
                                                     }
                                                     tmp += MatrixElement*
                                                            GlobalVectoronRight[GlobalIndex2][LocalVectorIndex2];
                                                 }
                                             }
                                             DistributedVector[LocalProcessIndex1][LocalVectorIndex1] = tmp;
                                         } // End of Varied parameter
                                     }
                                 }
                             }); // End loop over Point dependent quantities
        }

        // End MatrixGlobalVector(double[][] DistributedVector, double[][,] MatrixDiagonals,double[][] GlobalxVector, double[][] GlobalVectoronRight)


        // Solve Matrix Equations
        //  The solution is rescaled so correct for native matrix
        //  However Diagonally scaled matrix used in solver as long as Hotsun.UseDiagonalScalinginSolvers = true
        public static bool SolveMatrix(double[][] Answer, Desertwind Solution)
        {
            // Scale RHS Vector
            MDSLinearAlgebra.DiagScaleVector(Hotsun.UtilityLocalVector1, Solution.first, Hotsun.sqdginv);

            //  Solve Scaled Matrix Equations
            int RealMatrixSize = Hotsun.npar - ((Hotsun.ParameterVectorDimension + 1)*Hotsun.ParameterVectorDimension)/2;
            RealMatrixSize = Math.Min(RealMatrixSize, Hotsun.CGIterationLimit);
            bool matrixsuccess;
            matrixsuccess = MDSLinearAlgebra.ConjugateGradientSolver(Answer, Solution, Hotsun.FullSecondDerivative,
                                                                     Hotsun.GlobalParameter, Hotsun.UtilityLocalVector1,
                                                                     ref Hotsun.NumberofCGIterations, RealMatrixSize,
                                                                     Hotsun.CGResidualLimit);

            if (!matrixsuccess)
                ++Hotsun.TotalCGFailures;
            Hotsun.TotalCGIterations += Hotsun.NumberofCGIterations;


            //  Correct answer for diagonal scaling
            MDSLinearAlgebra.DiagScaleVector(Answer, Answer, Hotsun.sqdginv);
            return matrixsuccess;
        }

        // End SolveMatrix(double[][] Answer, Desertwind Solution)

        public static void FindQlimits(Desertwind Solution, ref double Qhigh, ref double Qlow, ref int ReasontoStop1,
                                       ref int ReasontoStop2)
        {
            if (Hotsun.FullSecondDerivative)
                MDSLinearAlgebra.FindTraceandNorm(Solution.ExactDiagonalofMatrix, Hotsun.GlobalParameter,
                                                  ref Hotsun.ChisqMatrixTrace, ref Hotsun.ChisqMatrixNorm);
            else
                MDSLinearAlgebra.FindTraceandNorm(Solution.DiagonalofMatrix, Hotsun.GlobalParameter,
                                                  ref Hotsun.ChisqMatrixTrace, ref Hotsun.ChisqMatrixNorm);

            Qhigh = 0.0;
            Qlow = 0.0;
            double PowerEigenvalue1 = 0.0;
            ReasontoStop2 = -3;

            ReasontoStop1 = MDSLinearAlgebra.PowerIterate(Solution, 0, 0.0, out PowerEigenvalue1);

            if (ReasontoStop1 > 0)
                Hotsun.TotalPowerIterations += ReasontoStop1;
            else
                ++Hotsun.TotalPowerFailures;

            if (ReasontoStop1 > 0)
            {
                Qhigh = PowerEigenvalue1 - Hotsun.addonforQcomputation;
            }
            else
            {
                Qhigh = Math.Min(Hotsun.ChisqMatrixNorm, Hotsun.ChisqMatrixTrace);
                Qlow = Hotsun.Qrange*Qhigh;
                return;
            }

            double PowerEigenvalue2 = 0.0;
            ReasontoStop2 = MDSLinearAlgebra.PowerIterate(Solution, 1, Qhigh, out PowerEigenvalue2);

            if (ReasontoStop2 > 0)
                Hotsun.TotalPowerIterations += ReasontoStop2;
            else
                Hotsun.TotalPowerFailures += 1;

            if (ReasontoStop2 > 0)
            {
                Qlow = Qhigh - PowerEigenvalue2;
            }
            else
            {
                Qlow = Hotsun.Qrange*Qhigh;
            }
            bool Qlowtoosmall = Qlow < (Hotsun.Qrange*Qhigh);
            SALSAUtility.SynchronizeMPIvariable(ref Qlowtoosmall);

            if (Qlowtoosmall)
                Qlow = Hotsun.Qrange*Qhigh;
            return;
        }

        // End FindQlimits(Desertwind Solution, ref double Qhigh, ref double Qlow, ref int ReasontoStop1, ref int ReasontoStop2)

        public static void FillupHotsun()
        {
            if (ManxcatCentral.Configuration.InitializationOption != 1 ||
                !File.Exists(ManxcatCentral.Configuration.InitializationFileName))
            {
                throw new Exception("Simple initialization file necessary with the specified processing option");
            }

            using (
                var reader =
                    new StreamReader(File.Open(ManxcatCentral.Configuration.InitializationFileName, FileMode.Open,
                                               FileAccess.Read, FileShare.Read)))
            {
                var sep = new[] {' ', '\t'};
                string[] splits;
                Hotsun.GlobalParameter = new double[SALSAUtility.PointCount_Global][];
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    if (!string.IsNullOrEmpty(line))
                    {
                        splits = line.Trim().Split(sep);
                        if (splits.Length == 5)
                        {
                            int originalIndex = int.Parse(splits[0]);
                            var vector = new[]
                                             {
                                                 double.Parse(splits[1]), double.Parse(splits[2]),
                                                 double.Parse(splits[3])
                                             };
                            if (originalIndex >= 0 && originalIndex < SALSAUtility.OriginalPointtoUsedPointMap.Length)
                            {
                                int globalIndex = SALSAUtility.OriginalPointtoUsedPointMap[originalIndex];
                                int usedIndex = SALSAUtility.NaivetoActualUsedOrder[globalIndex];
                                Hotsun.GlobalParameter[usedIndex] = vector;
                            }
                        }
                    }
                }
            }
        }
    }

    // End class ManxcatMDS
}

// End Namespace MDS