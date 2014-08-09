using System;
using System.Collections;
using System.Threading.Tasks;
using MPI;
using Manxcat;
using SALSALibrary;
using Salsa.Core;

namespace MDS
{
    public class ManxcatMDSBasicDataProcessing
    {
        // Initial distance analysis looped over till no more deletions
        public static void IntialDistanceAnalysis(ref double AverageDistance, ref double MaxDistance,
                                                  ref double SigmaDistance, out int TotalDisconnectedPoints,
                                                  out double TotalMissingDistances)
        {
            var FindMissingDistances = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
            var FindDisconnectedPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
            var FindSystemMax = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);
            var FindSystemMeanSigma = new GlobalReductions.FindMeanSigma(SALSAUtility.ThreadCount);
            var FindLinkswithoutCut = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount); // Links without Cut
            var FindLinkswithCut = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount); // Links with cut

            int HistogramSize = ManxcatMDS.NumLinksHistogramBins;
            var FindLinkHistogramBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                     HistogramSize);

            var DeletedPoint = new int[SALSAUtility.PointCount_Global]; // Set to zero if Point used = 1 if deleted
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                DeletedPoint[GlobalPointIndex] = 0;
            }

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 FindLinkHistogramBinCounts.startthread(ThreadNo);

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
                                     int OriginalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;
                                     int Countlinks = 0;
                                     for (int GlobalPointIndex2 = 0;
                                          GlobalPointIndex2 < SALSAUtility.PointCount_Global;
                                          GlobalPointIndex2++)
                                     {
                                         int OriginalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;
                                         if (GlobalPointIndex1 == GlobalPointIndex2)
                                             continue;
                                         double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1,
                                                                                         GlobalPointIndex2);
                                         if (tmp2 < -0.5)
                                         {
                                             FindMissingDistances.addapoint(ThreadNo, 1.0);
                                             continue;
                                         }
                                         // Unset Distances if either point on deleted list
                                         if ((ManxcatMDS.PointStatus[GlobalPointIndex1] == -1) ||
                                             (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1))
                                         {
                                             FindMissingDistances.addapoint(ThreadNo, 1.0);
                                             SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                               -1.0);
                                             continue;
                                         }
                                         FindSystemMeanSigma.addapoint(ThreadNo, tmp2);
                                         FindSystemMax.addapoint(ThreadNo, tmp2);
                                         ++Countlinks;
                                     }

                                     if (Countlinks <= SALSAUtility.LinkCut)
                                     {
                                         FindDisconnectedPoints.addapoint(ThreadNo, 1);
                                         DeletedPoint[GlobalPointIndex1] = 1;
                                     }
                                     if (Countlinks < HistogramSize)
                                     {
                                         // Histogram low Count Links
                                         FindLinkHistogramBinCounts.addapoint(ThreadNo, Countlinks);
                                     }
                                     FindLinkswithoutCut.addapoint(ThreadNo, Countlinks);
                                     if (Countlinks > SALSAUtility.LinkCut)
                                     {
                                         FindLinkswithCut.addapoint(ThreadNo, Countlinks);
                                     }
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindMissingDistances.sumoverthreadsandmpi();
            FindDisconnectedPoints.sumoverthreadsandmpi();
            FindSystemMax.sumoverthreadsandmpi();
            FindSystemMeanSigma.sumoverthreadsandmpi();
            FindLinkswithCut.sumoverthreadsandmpi();
            FindLinkswithoutCut.sumoverthreadsandmpi();
            FindLinkHistogramBinCounts.sumoverthreadsandmpi();

            // Reconcile Deleted points across all processes
            if (SALSAUtility.MPI_Size > 1)
            {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIREDUCETiming);
                DeletedPoint = SALSAUtility.MPI_communicator.Allreduce(DeletedPoint, Operation<int>.Add);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIREDUCETiming);
            }
            int countdeleted = 0;
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                if (DeletedPoint[GlobalPointIndex] == 0)
                    continue;
                if (DeletedPoint[GlobalPointIndex] != 1)
                {
                    Exception e =
                        SALSAUtility.SALSAError("System Error Inconsistent Deleted point " + GlobalPointIndex.ToString() +
                                                " Status " + DeletedPoint[GlobalPointIndex].ToString());
                    throw (e);
                }
                ManxcatMDS.PointStatus[GlobalPointIndex] = -1;
                ++countdeleted;
            }

            TotalDisconnectedPoints = FindDisconnectedPoints.TotalInt;
            TotalMissingDistances = FindMissingDistances.Total;

            double tmp = SALSAUtility.NumberVariedPoints;
            double pointsused = tmp*(tmp - 1) - TotalMissingDistances;
            double localtotalpoints = FindSystemMeanSigma.TotalNumberofPoints;
            if (Math.Abs(localtotalpoints - pointsused) > 2.0)
            {
                Exception e =
                    SALSAUtility.SALSAError("System Error and Must stop as Illegal Point Counts " +
                                            pointsused.ToString("F0") + " " + localtotalpoints.ToString("F0"));
                throw (e);
            }

            AverageDistance = FindSystemMeanSigma.Totalmean;
            SigmaDistance = FindSystemMeanSigma.Totalsigma;
            MaxDistance = FindSystemMax.TotalMax;


            // Output Link Statistics
            double TotalLinkswithoutCut = FindLinkswithoutCut.Total;
            double TotalLinkswithCut = FindLinkswithCut.Total;
            double TotalPointsbeforeLinkCut = FindLinkswithoutCut.TotalNumberofPoints;
            double TotalPointsafterLinkCut = FindLinkswithCut.TotalNumberofPoints;
            double AverageLinksperpoint = TotalLinkswithoutCut/TotalPointsbeforeLinkCut;
            double AverageLinksperpointaftercut = TotalLinkswithCut/TotalPointsafterLinkCut;

            if ((SALSAUtility.DebugPrintOption > 0) && (SALSAUtility.MPI_Rank == 0))
            {
                double fraction = TotalMissingDistances/(tmp*(tmp - 1.0));
                double EstimatedDimension = 2.0*AverageDistance*AverageDistance/(SigmaDistance*SigmaDistance);
                SALSAUtility.SALSAPrint(1,
                                        "\nAs INPUT Disconnected Points " + TotalDisconnectedPoints.ToString() +
                                        " Missing Distances " + TotalMissingDistances.ToString("F0")
                                        + " Fraction " + fraction.ToString("F4"));
                SALSAUtility.SALSAPrint(1,
                                        "As INPUT Max " + MaxDistance.ToString("E4") + " Average " +
                                        AverageDistance.ToString("E4") + " Sigma "
                                        + SigmaDistance.ToString("E4") + " Estimated Dimension " +
                                        EstimatedDimension.ToString("F2"));
                SALSAUtility.SALSAPrint(1,
                                        " Overall Average Links " + AverageLinksperpoint.ToString("F2") + " from " +
                                        TotalPointsbeforeLinkCut.ToString("F0")
                                        + " Points and After cut " + AverageLinksperpointaftercut.ToString("F2") +
                                        " From Points " + TotalPointsafterLinkCut.ToString("F0"));

                //  Output Link Histogram
                string histogramcounts = "";
                double total = 0.0;
                for (int binloop = 0; binloop < ManxcatMDS.NumLinksHistogramBins; binloop++)
                {
                    histogramcounts += FindLinkHistogramBinCounts.TotalSum[binloop].ToString("F0") + ", ";
                    total += FindLinkHistogramBinCounts.TotalSum[binloop];
                }
                SALSAUtility.SALSAPrint(1,
                                        "\nLink Histogram Total " + total.ToString("F0") + " upto links " +
                                        ManxcatMDS.NumLinksHistogramBins.ToString() + " #Counts starting at zero\n" +
                                        histogramcounts);

                string deletedpoints = "\nNumber Deleted Points " + countdeleted.ToString() + " Deleted Points ";
                for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
                {
                    if (ManxcatMDS.PointStatus[GlobalPointIndex] != -1)
                        continue;
                    deletedpoints += GlobalPointIndex.ToString() + ", ";
                }
                SALSAUtility.SALSAPrint(1, deletedpoints);
            }
            return;
        }

        //End IntialDistanceAnalysis

        // Second distance analysis designed to have same interface as IntialDistanceAnalysis
        public static void CleanandTransformDistances(ref double AverageDistance, ref double MaxDistce,
                                                      ref double SigmaDistance, out int TotalDisconnectedPoints,
                                                      out double TotalMissingDistances)
        {
            double initialAverageDistance = AverageDistance;
            double initialMaxDistance = MaxDistce;
            double initialSigmaDistance = SigmaDistance;

            double initialEstimatedDimension = 0.0;
            double initialScaleFactor = 0.0,
                   initialTransformedMaxDistance = 0.0;

            if (SALSAUtility.TransformMethod == 8 || SALSAUtility.TransformMethod == 9)
            {
                // 4D or SQRT(4D) transformation specific statistics
                initialEstimatedDimension = 2.0*initialAverageDistance*initialAverageDistance/
                                            (initialSigmaDistance*initialSigmaDistance);
                double initialIndividualSigma = Math.Sqrt(initialAverageDistance/initialEstimatedDimension);
                initialScaleFactor = 2.0*initialIndividualSigma*initialIndividualSigma;

                initialTransformedMaxDistance = initialScaleFactor*
                                                Transform4D(SpecialFunction.igamc(initialEstimatedDimension*0.5,
                                                                                  initialMaxDistance/initialScaleFactor));
                if (SALSAUtility.TransformMethod == 9)
                {
                    initialTransformedMaxDistance = Math.Sqrt(initialTransformedMaxDistance);
                }
            }

            if (SALSAUtility.NumberDistanceWeightCuts > 0)
            {
                for (int weightbin = 0; weightbin <= SALSAUtility.NumberDistanceWeightCuts; weightbin++)
                {
                    SALSAUtility.ActualWeightCuts[weightbin] = 0.0;
                }
            }
            var FindMissingDistances = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
            var FindDisconnectedPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
            var FindSystemMax = new GlobalReductions.FindDoubleMax(SALSAUtility.ThreadCount);
            var FindSystemMeanSigma = new GlobalReductions.FindMeanSigma(SALSAUtility.ThreadCount);
            var FindNewoldDistancestatistics = new GlobalReductions.FindCorrelation(SALSAUtility.ThreadCount);
            var FindDistanceWeightBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                      1 +
                                                                                      SALSAUtility.
                                                                                          NumberDistanceWeightCuts);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;
                                 FindDistanceWeightBinCounts.startthread(ThreadNo);

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
                                     int OriginalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;
                                     int Countlinks = 0;
                                     for (int GlobalPointIndex2 = 0;
                                          GlobalPointIndex2 < SALSAUtility.PointCount_Global;
                                          GlobalPointIndex2++)
                                     {
                                         int OriginalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;
                                         if (GlobalPointIndex1 == GlobalPointIndex2)
                                             continue;
                                         double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1,
                                                                                         GlobalPointIndex2);
                                         if (tmp2 < -0.5)
                                         {
                                             FindMissingDistances.addapoint(ThreadNo, 1.0);
                                             continue;
                                         }
                                         // Unset Distances if either point on deleted list
                                         if ((ManxcatMDS.PointStatus[GlobalPointIndex1] == -1) ||
                                             (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1))
                                         {
                                             FindMissingDistances.addapoint(ThreadNo, 1.0);
                                             SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                               -1.0);
                                         }
                                         else
                                         {
                                             double newvalue = tmp2;

                                             // Transform Distances if Requested
                                             if (SALSAUtility.TransformMethod == 12)
                                             {
                                                 // log(1-d) transformation
                                                 tmp2 = Math.Min(0.975, tmp2);
                                                 newvalue =
                                                     Math.Log(1 - Math.Pow(tmp2, SALSAUtility.TransformParameter))/
                                                     Math.Log(1 - Math.Pow(0.975, SALSAUtility.TransformParameter));
                                                 SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                                   newvalue);
                                             }
                                             else if (SALSAUtility.TransformMethod == 11)
                                             {
                                                 // (1-x) ** power with power given by TransformParameter
                                                 tmp2 = Math.Min(1.0, tmp2);
                                                 newvalue = 1.0 - Math.Pow(1.0 - tmp2, SALSAUtility.TransformParameter);
                                                 SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                                   newvalue);
                                             }
                                             else if (SALSAUtility.TransformMethod == 10)
                                             {
                                                 // (x ** power with power given by TransformParameter
                                                 tmp2 = Math.Min(1.0, tmp2);
                                                 newvalue = Math.Pow(tmp2, SALSAUtility.TransformParameter);
                                                 SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                                   newvalue);
                                             }
                                             else if (SALSAUtility.TransformMethod == 8)
                                             {
                                                 // 4D Transformation. 
                                                 // Note. Works only when no missing distances
                                                 tmp2 = initialScaleFactor*
                                                        Transform4D(SpecialFunction.igamc(
                                                            initialEstimatedDimension*0.5, tmp2/initialScaleFactor));
                                                 newvalue = tmp2/initialTransformedMaxDistance;
                                                 SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                                   newvalue);
                                             }
                                             else if (SALSAUtility.TransformMethod == 9)
                                             {
                                                 // SQRT(4D) Transformation. 
                                                 // Note. Works only when no missing distances
                                                 tmp2 = initialScaleFactor*
                                                        Transform4D(SpecialFunction.igamc(
                                                            initialEstimatedDimension*0.5, tmp2/initialScaleFactor));
                                                 newvalue = Math.Sqrt(tmp2)/initialTransformedMaxDistance;
                                                 SALSAParallelism.putDistanceValue(GlobalPointIndex1, GlobalPointIndex2,
                                                                                   newvalue);
                                             }

                                             FindNewoldDistancestatistics.addapoint(ThreadNo, tmp2, newvalue);
                                             FindSystemMax.addapoint(ThreadNo, newvalue);
                                             FindSystemMeanSigma.addapoint(ThreadNo, newvalue);
                                             ++Countlinks;

                                             if (SALSAUtility.NumberDistanceWeightCuts > 0)
                                             {
                                                 int distpos = SALSAUtility.NumberDistanceWeightCuts;
                                                 for (int weightbin = 0;
                                                      weightbin < SALSAUtility.NumberDistanceWeightCuts;
                                                      weightbin++)
                                                 {
                                                     if (SALSAUtility.DistanceWeightCuts[weightbin] > newvalue)
                                                     {
                                                         distpos = weightbin;
                                                         break;
                                                     }
                                                 }
                                                 FindDistanceWeightBinCounts.addapoint(ThreadNo, distpos);
                                             }
                                         }
                                     }

                                     if ((Countlinks <= 0) && (ManxcatMDS.PointStatus[GlobalPointIndex1] != -1))
                                     {
                                         FindDisconnectedPoints.addapoint(ThreadNo, 1);
                                     }
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            // Accumulate Threads and invoke MPI
            FindNewoldDistancestatistics.sumoverthreadsandmpi();
            FindSystemMeanSigma.sumoverthreadsandmpi();
            FindSystemMax.sumoverthreadsandmpi();
            FindDisconnectedPoints.sumoverthreadsandmpi();
            FindMissingDistances.sumoverthreadsandmpi();
            FindDistanceWeightBinCounts.sumoverthreadsandmpi();

            AverageDistance = FindSystemMeanSigma.Totalmean;
            SigmaDistance = FindSystemMeanSigma.Totalsigma;
            MaxDistce = FindSystemMax.TotalMax;
            TotalDisconnectedPoints = FindDisconnectedPoints.TotalInt;
            TotalMissingDistances = FindMissingDistances.Total;
            double tmp = SALSAUtility.NumberVariedPoints;
            double pointsused = tmp*(tmp - 1) - TotalMissingDistances;
            double localtotalpoints = FindSystemMeanSigma.TotalNumberofPoints;
            if (Math.Abs(localtotalpoints - pointsused) > 2.0)
            {
                Exception e =
                    SALSAUtility.SALSAError("System Error and Must stop as Illegal Point Counts " +
                                            pointsused.ToString("F0") + " " + localtotalpoints.ToString("F0"));
                throw (e);
            }

            // Print correlations
            string label = "";
            if (SALSAUtility.TransformMethod != 0)
            {
                label = "\nTransform Method " + SALSAUtility.TransformMethod.ToString() + " Parameter " +
                        SALSAUtility.TransformParameter.ToString("F4") + " ";
            }

            FindNewoldDistancestatistics.print(label + "\nOld and Transformed Distances\n", "F4");
            if (TotalDisconnectedPoints > 0)
            {
                Exception e =
                    SALSAUtility.SALSAError(
                        "System Error and Must stop as Disconnected Points with ZERO links after first round " +
                        TotalDisconnectedPoints.ToString());
                throw (e);
            }

            // Set Distance Weightings
            if (SALSAUtility.NumberDistanceWeightCuts <= 0)
                return;
            string distanceresults = "";
            double fudge = 1.0 + SALSAUtility.NumberDistanceWeightCuts;
            for (int weightbin = 0; weightbin <= SALSAUtility.NumberDistanceWeightCuts; weightbin++)
            {
                double weightforbin = 0.0;
                if (FindDistanceWeightBinCounts.TotalSum[weightbin] > 0.5)
                {
                    weightforbin = pointsused/(fudge*FindDistanceWeightBinCounts.TotalSum[weightbin]);
                }
                string labelreason = " Rest";
                if (weightbin < SALSAUtility.NumberDistanceWeightCuts)
                    labelreason = " Distce < " + SALSAUtility.DistanceWeightCuts[weightbin].ToString("F3");
                distanceresults += FindDistanceWeightBinCounts.TotalSum[weightbin].ToString("F0") + " Pts weight " +
                                   weightforbin.ToString("F3") + labelreason + " * ";
                SALSAUtility.ActualWeightCuts[weightbin] = Math.Sqrt(weightforbin);
            }
            SALSAUtility.SALSAPrint(1, "\nDistance Counts and Cuts " + distanceresults);
            return;
        }

        private static double Transform4D(double higherDimensionInput)
        {
            /* Saliya - copying code directly from Adam Hughes project. 
             * The comments except this one are from his code as well.*/

            double Einit = -Math.Log(higherDimensionInput);
            //         double Einit = HigherDimInput;
            double E = Einit;
            double Eold = E;
            // double diff;

            for (int recurse = 0; recurse < 50; recurse++)
            {
                E = Einit + Math.Log(1.0 + E);
                //      E = 1.0 + E - E / 2;
                /*      diff = E - Eold;
                      if (diff < 0)
                      {
                          diff = Eold - E;
                      }*/
                //                if (diff < 0.00001)
                if (Math.Abs(E - Eold) < 0.00001)
                    return E;
                Eold = E;
            }
            return E;
        }

//End CleanandTransformDistances

        // Find Center defined as Point with minimum mean distance to other points
        // Form Histogram of distance Bin counts
        public static void FindCenter(ref int Center, ref double Radius, double Histmin, double Histmax, int NumberBins,
                                      ref double[] Bincounts)
        {
            var FindCenterCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);

            int HistogramSize = 2 + NumberBins;
            double HistFudge = NumberBins/(Histmax - Histmin);
            var FindDistanceHistogramBinCounts = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                         HistogramSize);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;
                                 FindDistanceHistogramBinCounts.startthread(ThreadNo);

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     double distcemeanperpoint = 0.0;
                                     int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
                                     if (ManxcatMDS.PointStatus[GlobalPointIndex1] == -1)
                                         continue;
                                     int OriginalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;
                                     int Countlinks = 0;
                                     for (int GlobalPointIndex2 = 0;
                                          GlobalPointIndex2 < SALSAUtility.PointCount_Global;
                                          GlobalPointIndex2++)
                                     {
                                         if (ManxcatMDS.PointStatus[GlobalPointIndex2] == -1)
                                             continue;
                                         int OriginalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;
                                         if (GlobalPointIndex1 == GlobalPointIndex2)
                                             continue;
                                         double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1,
                                                                                         GlobalPointIndex2);
                                         if (tmp2 < -0.5)
                                             continue;
                                         distcemeanperpoint += tmp2;
                                         ++Countlinks;

                                         // Form histogram
                                         double Histvalue = (tmp2 - Histmin)*HistFudge;
                                         int HistPosition = Convert.ToInt32(Math.Floor(Histvalue));
                                         if (HistPosition == NumberBins)
                                         {
                                             HistPosition = Convert.ToInt32(Math.Floor(Histvalue - 0.00001));
                                         }
                                         if (HistPosition < 0)
                                             HistPosition = -1;
                                         if (HistPosition >= NumberBins)
                                             HistPosition = NumberBins;
                                         FindDistanceHistogramBinCounts.addapoint(ThreadNo, 1 + HistPosition);
                                     }
                                     if ((Countlinks >= ManxcatMDS.LinkCutforCenter) &&
                                         (ManxcatMDS.PointStatus[GlobalPointIndex1] == 0))
                                     {
                                         distcemeanperpoint = distcemeanperpoint/Countlinks;
                                         FindCenterCompute.addapoint(ThreadNo, GlobalPointIndex1, distcemeanperpoint);
                                     }
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindCenterCompute.sumoverthreadsandmpi();
            FindDistanceHistogramBinCounts.sumoverthreadsandmpi();
            Center = FindCenterCompute.TotalIndexValue;
            Radius = FindCenterCompute.TotalMaxOrMin;
            for (int binloop = 0; binloop < HistogramSize; binloop++)
                Bincounts[binloop] = FindDistanceHistogramBinCounts.TotalSum[binloop];
            return;
        }

        //End FindCenter

        /// <summary>
        /// Finds the minimum and maximum distances of original and Eculidean data
        /// for both whole and selected set of clusters.
        /// </summary>
        /// <param name="xminWhole"></param>
        /// <param name="xmaxWhole"></param>
        /// <param name="yminWhole"></param>
        /// <param name="ymaxWhole"></param>
        /// <param name="xminSelected"></param>
        /// <param name="xmaxSelected"></param>
        /// <param name="yminSelected"></param>
        /// <param name="ymaxSelected"></param>
        /// <param name="ymaxSelectedInter"></param>
        /// <param name="originalPnumToCnumTable"></param>
        /// <param name="xminSelectedInter"></param>
        /// <param name="xmaxSelectedInter"></param>
        /// <param name="yminSelectedInter"></param>
        public static void FindMinMaxForDensity(ref double xminWhole, ref double xmaxWhole, ref double yminWhole,
                                                ref double ymaxWhole,
                                                ref double xminSelected, ref double xmaxSelected,
                                                ref double yminSelected, ref double ymaxSelected,
                                                ref double xminSelectedInter, ref double xmaxSelectedInter,
                                                ref double yminSelectedInter, ref double ymaxSelectedInter,
                                                Hashtable originalPnumToCnumTable)
        {
            /* MinorMax objects with index though index is not used here */
            var FindXminWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
            var FindXminSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
            var FindXminSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
            var FindYminWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
            var FindYminSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);
            var FindYminSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 0);

            var FindXmaxWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
            var FindXmaxSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
            var FindXmaxSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
            var FindYmaxWhole = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
            var FindYmaxSelected = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);
            var FindYmaxSelectedInter = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int distributedPointIndex = beginpoint;
                                      distributedPointIndex < indexlen + beginpoint;
                                      distributedPointIndex++)
                                 {
                                     int globalPointIndex1 = distributedPointIndex + SALSAUtility.PointStart_Process;
                                     if (ManxcatMDS.PointStatus[globalPointIndex1] == -1)
                                         continue;

                                     int originalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex1];
                                     int cnum1 = SALSAUtility.IsClustersSelected
                                                     ? ((int) originalPnumToCnumTable[originalPointIndex1])
                                                     : -1;


                                     if (SALSAUtility.OriginalPointDisposition[originalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;

                                     int usedPointIndex1 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex1];

                                     for (int globalPointIndex2 = 0;
                                          globalPointIndex2 < SALSAUtility.PointCount_Global;
                                          globalPointIndex2++)
                                     {
                                         if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
                                             continue;

                                         int originalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[originalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;

                                         if (globalPointIndex1 == globalPointIndex2)
                                             continue;

                                         double xval = SALSAParallelism.getDistanceValue(globalPointIndex1,
                                                                                         globalPointIndex2);
                                         if (xval < -0.5)
                                             continue;

                                         int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
                                         double yval = GetEuclideanDistance(Hotsun.GlobalParameter[usedPointIndex1],
                                                                            Hotsun.GlobalParameter[usedPointIndex2]);

                                         /* At this point the xval should be the transformed distances (if specified using TransformMethod & TransformParameter) */

                                         // Todo (html+density) - see if any distance cut needs to be considered. Also check if any pair is to be removed if not exist in pnumToCnum table

                                         UpdateMinMax(xval, yval, FindXminWhole, FindXmaxWhole, FindYminWhole,
                                                      FindYmaxWhole, ThreadNo);

                                         int cnum2 = SALSAUtility.IsClustersSelected
                                                         ? ((int) originalPnumToCnumTable[originalPointIndex2])
                                                         : -1;
                                         if (cnum1 != -1 && cnum2 != -1 &&
                                             SALSAUtility.SelectedClusters.Contains(cnum1) &&
                                             SALSAUtility.SelectedClusters.Contains(cnum2))
                                         {
                                             if (cnum1 == cnum2)
                                             {
                                                 // Intra cluster pairs (p1,p2) where both p1,p2 belong to one cluster.
                                                 // So check if p1,p2 belong to the same cluster in our set of selected clusters
                                                 UpdateMinMax(xval, yval, FindXminSelected, FindXmaxSelected,
                                                              FindYminSelected,
                                                              FindYmaxSelected, ThreadNo);
                                             }
                                             else
                                             {
                                                 // Inter cluster pairs (p1,p2) where both p1,p2 does NOT belong to one cluster.
                                                 // So check if p1,p2 does NOT belong to the same cluster in our set of selected clusters
                                                 UpdateMinMax(xval, yval, FindXminSelectedInter, FindXmaxSelectedInter,
                                                              FindYminSelectedInter,
                                                              FindYmaxSelectedInter, ThreadNo);
                                             }
                                         }
                                     }
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindXminWhole.sumoverthreadsandmpi();
            FindXminSelected.sumoverthreadsandmpi();
            FindXminSelectedInter.sumoverthreadsandmpi();
            FindYminWhole.sumoverthreadsandmpi();
            FindYminSelected.sumoverthreadsandmpi();
            FindYminSelectedInter.sumoverthreadsandmpi();

            FindXmaxWhole.sumoverthreadsandmpi();
            FindXmaxSelected.sumoverthreadsandmpi();
            FindXmaxSelectedInter.sumoverthreadsandmpi();
            FindYmaxWhole.sumoverthreadsandmpi();
            FindYmaxSelected.sumoverthreadsandmpi();
            FindYmaxSelectedInter.sumoverthreadsandmpi();

            xminWhole = FindXminWhole.TotalMaxOrMin;
            xminSelected = FindXminSelected.TotalMaxOrMin;
            xminSelectedInter = FindXminSelectedInter.TotalMaxOrMin;
            yminWhole = FindYminWhole.TotalMaxOrMin;
            yminSelected = FindYminSelected.TotalMaxOrMin;
            yminSelectedInter = FindYminSelectedInter.TotalMaxOrMin;

            xmaxWhole = FindXmaxWhole.TotalMaxOrMin;
            xmaxSelected = FindXmaxSelected.TotalMaxOrMin;
            xmaxSelectedInter = FindXmaxSelectedInter.TotalMaxOrMin;
            ymaxWhole = FindYmaxWhole.TotalMaxOrMin;
            ymaxSelected = FindYmaxSelected.TotalMaxOrMin;
            ymaxSelectedInter = FindYmaxSelectedInter.TotalMaxOrMin;
        }

        /// <summary>
        /// Given x and y values along their respective current minimum and maximum values, this will update
        /// the minimum and maximum values.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xmin"></param>
        /// <param name="xmax"></param>
        /// <param name="ymin"></param>
        /// <param name="ymax"></param>
        /// <param name="threadNo"></param>
        private static void UpdateMinMax(double x, double y, GlobalReductions.FindMinorMaxValuewithIndex xmin,
                                         GlobalReductions.FindMinorMaxValuewithIndex xmax,
                                         GlobalReductions.FindMinorMaxValuewithIndex ymin,
                                         GlobalReductions.FindMinorMaxValuewithIndex ymax, int threadNo)
        {
            // Ignore the index value
            xmin.addapoint(threadNo, 0, x);
            xmax.addapoint(threadNo, 0, x);
            ymin.addapoint(threadNo, 0, y);
            ymax.addapoint(threadNo, 0, y);
        }

        /// <summary>
        /// Generates density matrices of Euclidean Vs original distances for the whole
        /// sample and the points in selected set of clusters. In selected clusters we
        /// consider only the intra-cluster (within the cluster) pairs of points for distances.
        /// </summary>
        /// <param name="densityMatrixWhole">2D double array of size SALSAUtlity.Xres x SALSAUtility.Yres</param>
        /// <param name="yHistogramSelected"></param>
        /// <param name="densityMatrixSelected">2D double array of size SALSAUtlity.Xres x SALSAUtility.Yres</param>
        /// <param name="yHistogramSelectedInter"></param>
        /// <param name="xminWhole">minimum original distance of a pair in the whole sample</param>
        /// <param name="xmaxWhole">maximum original distance of a pair in the whole sample</param>
        /// <param name="yminWhole">minimum Euclidean distance of a pair in the whole sample</param>
        /// <param name="ymaxWhole">maximum Euclidean distance of a pair in the whole sample</param>
        /// <param name="xminSelected">minimum original distance of a intra-cluster pair in the selected clusters</param>
        /// <param name="xmaxSelected">maximum original distance of a intra-cluster pair in the selected clusters</param>
        /// <param name="yminSelected">minimum Euclidean distance of a intra-cluster pair in the selected clusters</param>
        /// <param name="ymaxSelected">maximum Euclidean distance of a intra-cluster pair in the selected clusters</param>
        /// <param name="countWhole"></param>
        /// <param name="countSelected"></param>
        /// <param name="originalPnumToCnumTable">mapping from original point number to cluster number</param>
        /// <param name="xHistogramWhole"></param>
        /// <param name="yHistogramWhole"></param>
        /// <param name="xHistogramSelected"></param>
        /// <param name="densityMatrixSelectedInter"></param>
        /// <param name="xHistogramSelectedInter"></param>
        public static void GenerateDensityMatrix(out double[][] densityMatrixWhole, out double[] xHistogramWhole,
                                                 out double[] yHistogramWhole,
                                                 out double[][] densityMatrixSelected, out double[] xHistogramSelected,
                                                 out double[] yHistogramSelected,
                                                 out double[][] densityMatrixSelectedInter,
                                                 out double[] xHistogramSelectedInter,
                                                 out double[] yHistogramSelectedInter,
                                                 double xminWhole, double xmaxWhole, double yminWhole, double ymaxWhole,
                                                 double xminSelected, double xmaxSelected, double yminSelected,
                                                 double ymaxSelected,
                                                 double xminSelectedInter, double xmaxSelectedInter,
                                                 double yminSelectedInter, double ymaxSelectedInter,
                                                 out double countWhole, out double countSelected,
                                                 out double countSelectedInter, Hashtable originalPnumToCnumTable)
        {
            var FindDensityMatrixWhole = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                   SALSAUtility.Yres, SALSAUtility.Xres);
            var FindDensityMatrixSelected = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                      SALSAUtility.Yres,
                                                                                      SALSAUtility.Xres);
            var FindDensityMatrixSelectedInter = new GlobalReductions.Find2DDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                           SALSAUtility.Yres,
                                                                                           SALSAUtility.Xres);

            var FindXHistogramWhole = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                              SALSAUtility.Xres);
            var FindXHistogramSelected = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                 SALSAUtility.Xres);
            var FindXHistogramSelectedInter = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                      SALSAUtility.Xres);
            var FindYHistogramWhole = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                              SALSAUtility.Yres);
            var FindYHistogramSelected = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                 SALSAUtility.Yres);
            var FindYHistogramSelectedInter = new GlobalReductions.FindDoubleArraySum(SALSAUtility.ThreadCount,
                                                                                      SALSAUtility.Yres);

            double deltaxWhole = (xmaxWhole - xminWhole)/SALSAUtility.Xres;
            double deltayWhole = (ymaxWhole - yminWhole)/SALSAUtility.Yres;

            double deltaxSelected = (xmaxSelected - xminSelected)/SALSAUtility.Xres;
            double deltaySelected = (ymaxSelected - yminSelected)/SALSAUtility.Yres;

            double deltaxSelectedInter = (xmaxSelectedInter - xminSelectedInter)/SALSAUtility.Xres;
            double deltaySelectedInter = (ymaxSelectedInter - yminSelectedInter)/SALSAUtility.Yres;

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 FindDensityMatrixWhole.startthread(ThreadNo);
                                 FindDensityMatrixSelected.startthread(ThreadNo);
                                 FindDensityMatrixSelectedInter.startthread(ThreadNo);

                                 FindXHistogramWhole.startthread(ThreadNo);
                                 FindXHistogramSelected.startthread(ThreadNo);
                                 FindXHistogramSelectedInter.startthread(ThreadNo);
                                 FindYHistogramWhole.startthread(ThreadNo);
                                 FindYHistogramSelected.startthread(ThreadNo);
                                 FindYHistogramSelectedInter.startthread(ThreadNo);

                                 for (int distributedPointIndex = beginpoint;
                                      distributedPointIndex < indexlen + beginpoint;
                                      distributedPointIndex++)
                                 {
                                     int globalPointIndex1 = distributedPointIndex + SALSAUtility.PointStart_Process;
                                     if (ManxcatMDS.PointStatus[globalPointIndex1] == -1)
                                         continue;

                                     int originalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex1];
                                     int cnum1 = SALSAUtility.IsClustersSelected
                                                     ? ((int) originalPnumToCnumTable[originalPointIndex1])
                                                     : -1;


                                     if (SALSAUtility.OriginalPointDisposition[originalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;

                                     int usedPointIndex1 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex1];

                                     for (int globalPointIndex2 = 0;
                                          globalPointIndex2 < SALSAUtility.PointCount_Global;
                                          globalPointIndex2++)
                                     {
                                         if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
                                             continue;

                                         int originalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[originalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;

                                         if (globalPointIndex1 == globalPointIndex2)
                                             continue;

                                         double xval = SALSAParallelism.getDistanceValue(globalPointIndex1,
                                                                                         globalPointIndex2);
                                         if (xval < -0.5)
                                             continue;

                                         int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
                                         double yval = GetEuclideanDistance(Hotsun.GlobalParameter[usedPointIndex1],
                                                                            Hotsun.GlobalParameter[usedPointIndex2]);

                                         /* At this point the xval should be the transformed distances (if specified using TransformMethod & TransformParameter) */

                                         // Todo (html+density) - see if any distance cut needs to be considered. Also check if any pair is to be removed if not exist in pnumToCnum table

                                         UpdateCells(xval, yval, xmaxWhole, xminWhole, ymaxWhole, yminWhole, deltaxWhole,
                                                     deltayWhole,
                                                     FindDensityMatrixWhole, FindXHistogramWhole, FindYHistogramWhole,
                                                     ThreadNo);

                                         int cnum2 = SALSAUtility.IsClustersSelected
                                                         ? ((int) originalPnumToCnumTable[originalPointIndex2])
                                                         : -1;
                                         if (cnum1 != -1 && cnum2 != -1 && SALSAUtility.SelectedClusters.Contains(cnum1) &&
                                             SALSAUtility.SelectedClusters.Contains(cnum2))
                                         {
                                             if (cnum1 == cnum2)
                                             {
                                                 // Intra cluster pairs (p1,p2) where both p1,p2 belong to one cluster.
                                                 // So check if p1,p2 belong to the same cluster in our set of selected clusters
                                                 UpdateCells(xval, yval, xmaxSelected, xminSelected, ymaxSelected,
                                                             yminSelected,
                                                             deltaxSelected, deltaySelected, FindDensityMatrixSelected,
                                                             FindXHistogramSelected, FindYHistogramSelected, ThreadNo);
                                             }
                                             else
                                             {
                                                 // Inter cluster pairs (p1,p2) where both p1,p2 does NOT belong to one cluster.
                                                 // So check if p1,p2 does NOT belong to the same cluster in our set of selected clusters
                                                 UpdateCells(xval, yval, xmaxSelectedInter, xminSelectedInter,
                                                             ymaxSelectedInter, yminSelectedInter,
                                                             deltaxSelectedInter, deltaySelectedInter,
                                                             FindDensityMatrixSelectedInter,
                                                             FindXHistogramSelectedInter, FindYHistogramSelectedInter,
                                                             ThreadNo);
                                             }
                                         }
                                     }
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindDensityMatrixWhole.sumoverthreadsandmpi();
            FindDensityMatrixSelected.sumoverthreadsandmpi();
            FindDensityMatrixSelectedInter.sumoverthreadsandmpi();

            FindXHistogramWhole.sumoverthreadsandmpi();
            FindXHistogramSelected.sumoverthreadsandmpi();
            FindXHistogramSelectedInter.sumoverthreadsandmpi();
            FindYHistogramWhole.sumoverthreadsandmpi();
            FindYHistogramSelected.sumoverthreadsandmpi();
            FindYHistogramSelectedInter.sumoverthreadsandmpi();

            densityMatrixWhole = FindDensityMatrixWhole.TotalSum;
            densityMatrixSelected = FindDensityMatrixSelected.TotalSum;
            densityMatrixSelectedInter = FindDensityMatrixSelectedInter.TotalSum;

            xHistogramWhole = FindXHistogramWhole.TotalSum;
            xHistogramSelected = FindXHistogramSelected.TotalSum;
            xHistogramSelectedInter = FindXHistogramSelectedInter.TotalSum;
            yHistogramWhole = FindYHistogramWhole.TotalSum;
            yHistogramSelected = FindYHistogramSelected.TotalSum;
            yHistogramSelectedInter = FindYHistogramSelectedInter.TotalSum;

            countWhole = FindDensityMatrixWhole.TotalNumberofPoints;
            countSelected = FindDensityMatrixSelected.TotalNumberofPoints;
            countSelectedInter = FindDensityMatrixSelectedInter.TotalNumberofPoints;
        }

        private static void UpdateCells(double x, double y, double xmax, double xmin, double ymax, double ymin,
                                        double deltax, double deltay,
                                        GlobalReductions.Find2DDoubleArraySum FindDensityMatrix,
                                        GlobalReductions.FindDoubleArraySum FindXHistogram,
                                        GlobalReductions.FindDoubleArraySum FindYHistogram, int threadNo)
        {
            // cell number based on zero index from bottom left corner
            // if x is equal to xmax then it's placed in the last cell, which is xres-1 in zero based index
            // same is done for y when y == ymax
            int cellx = Math.Abs(x - xmax) < double.Epsilon
                            ? SALSAUtility.Xres - 1
                            : (int) Math.Floor((x - xmin)/deltax);
            int celly = Math.Abs(y - ymax) < double.Epsilon
                            ? SALSAUtility.Yres - 1
                            : (int) Math.Floor((y - ymin)/deltay);

            if (x > xmax || y > ymax || x < xmin || y < ymin)
            {
                // now this should never be reached
                throw new Exception("bad(1)-> x: " + x + " y: " + y + " xmax: " + xmax + " xmin: " + xmin + " ymax: " +
                                    ymax + " ymin: " + ymin);
            }

            if (cellx >= SALSAUtility.Xres || celly >= SALSAUtility.Yres)
            {
                // now this should never be reached
                throw new Exception("bad(2)-> x: " + x + " y:" + y + " xmax: " + xmax + " xmin: " + xmin + " ymax: " +
                                    ymax + " ymin: " + ymin + " cellx: " + cellx + " celly: " + celly);
            }

            FindDensityMatrix.addapoint(threadNo, celly, cellx);
            FindXHistogram.addapoint(threadNo, cellx);
            FindYHistogram.addapoint(threadNo, celly);
        }


        public static double GetEuclideanDistance(double[] vector1, double[] vector2)
        {
            if (vector1.Length != vector2.Length)
            {
                throw new Exception(ManxcatErrorMessages.UnequalVectorLengths);
            }

            int dimension = vector1.Length;

            double d = 0.0;
            for (int i = 0; i < dimension; ++i)
            {
                d += Math.Pow((vector2[i] - vector1[i]), 2);
            }
            return Math.Sqrt(d);
        }


        // Find x Axis position and number of nearby distances and number of points near another
        public static void FindxAxis(int Center, ref int xAxis, ref double MaxDistceGlobal,
                                     ref double DistancesNearEachOther, ref int NotLonelyPoints)
        {
            var FindNearbyPairs = new GlobalReductions.FindDoubleSum(SALSAUtility.ThreadCount);
            var FindCozyPoints = new GlobalReductions.FindIntSum(SALSAUtility.ThreadCount);
            var FindxAxisCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     int notlonely = 0;
                                     int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
                                     if (ManxcatMDS.PointStatus[GlobalPointIndex1] != 0)
                                         continue;
                                     int OriginalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;

                                     for (int GlobalPointIndex2 = 0;
                                          GlobalPointIndex2 < SALSAUtility.PointCount_Global;
                                          GlobalPointIndex2++)
                                     {
                                         if (GlobalPointIndex2 == GlobalPointIndex1)
                                             continue;
                                         int OriginalPointIndex2 =
                                             SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex2];
                                         if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex2] <
                                             SALSAUtility.SALSASHIFT)
                                             continue;
                                         double tmp = SALSAParallelism.getDistanceValue(GlobalPointIndex1,
                                                                                        GlobalPointIndex2);
                                         if (tmp < -0.5)
                                             continue;

                                         if (tmp < ManxcatMDS.MinimumDistance)
                                         {
                                             ++notlonely;
                                             FindNearbyPairs.addapoint(ThreadNo, 1.0);
                                         }
                                     }

                                     if (notlonely > 0)
                                         FindCozyPoints.addapoint(ThreadNo, 1);

                                     double MaxDistcePoint = SALSAParallelism.getDistanceValue(GlobalPointIndex1, Center);
                                     if (MaxDistcePoint < -0.5)
                                         continue;
                                     FindxAxisCompute.addapoint(ThreadNo, GlobalPointIndex1, MaxDistcePoint);
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindxAxisCompute.sumoverthreadsandmpi();
            FindCozyPoints.sumoverthreadsandmpi();
            FindNearbyPairs.sumoverthreadsandmpi();

            xAxis = FindxAxisCompute.TotalIndexValue;
            MaxDistceGlobal = FindxAxisCompute.TotalMaxOrMin;
            NotLonelyPoints = FindCozyPoints.TotalInt;
            DistancesNearEachOther = FindNearbyPairs.Total;
            return;
        }

        //End FindxAxis

        //  Find point defining xy Plane
        public static void FindxyPlane(int Center, int xAxis, ref int xyPlane, ref double MaxDistceGlobal)
        {
            var FindxyPlaneCompute = new GlobalReductions.FindMinorMaxValuewithIndex(SALSAUtility.ThreadCount, 1);

            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     int GlobalPointIndex1 = DistributedPointIndex + SALSAUtility.PointStart_Process;
                                     int OriginalPointIndex1 =
                                         SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex1];
                                     if (SALSAUtility.OriginalPointDisposition[OriginalPointIndex1] <
                                         SALSAUtility.SALSASHIFT)
                                         continue;

                                     if (ManxcatMDS.PointStatus[GlobalPointIndex1] != 0)
                                         continue;

                                     double tmp1 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, Center);
                                     if (tmp1 < -0.5)
                                         continue;
                                     double tmp2 = SALSAParallelism.getDistanceValue(GlobalPointIndex1, xAxis);
                                     if (tmp2 < -0.5)
                                         continue;
                                     double MaxDistcePoint = tmp1*tmp2;
                                     FindxyPlaneCompute.addapoint(ThreadNo, GlobalPointIndex1, MaxDistcePoint);
                                 } // End loop over points in this thread
                             }); // End loop over Point dependent quantities

            FindxyPlaneCompute.sumoverthreadsandmpi();
            xyPlane = FindxyPlaneCompute.TotalIndexValue;
            MaxDistceGlobal = FindxyPlaneCompute.TotalMaxOrMin;
            return;
        }

        //End FindxyPlane

        public static void SetUpHistogramRange(int PointsinHistogram, ref double Histmin, ref double Histmax)
        {
            // Choose good range to get rounded labels

            if (Histmax <= Histmin)
                return;
            if (Math.Abs(Histmin) < 0.000000001)
                Histmin = 0.0;
            else
            {
                double newminvalue = NewBinSize(Math.Abs(Histmin));
                if (Histmin < 0.0)
                    Histmin = -newminvalue;
                else
                    Histmin = newminvalue;
            }
            double binsize = (Histmax - Histmin)/PointsinHistogram;
            Histmax = Histmin + PointsinHistogram*NewBinSize(binsize);
            return;
        }

        // End SetUpHistogramRange

        public static double NewBinSize(double oldBinSize)
        {
            // Round Bin size up to a pretty value
            if (oldBinSize <= 0.000000001 || oldBinSize > 10000000000.0)
                return oldBinSize;

            double logvalue = Math.Log10(oldBinSize);
            int intlogvalue = Convert.ToInt32(Math.Floor(logvalue));
            double fudgepower = 1.0 - Convert.ToDouble(intlogvalue);
            double fudge = Math.Pow(10.0, fudgepower);
            double scaled = fudge*oldBinSize;
            scaled = Math.Min(scaled, 100.0);
            double Intversionofnewbinsize = Math.Ceiling(scaled)/fudge;
//            SALSAUtility.SALSAPrint(1, "Hist " + oldBinSize.ToString("F4") + " " + intlogvalue + " " + Intversionofnewbinsize.ToString("F4") + " scaled " + scaled.ToString("F2"));
            return Intversionofnewbinsize;
        }
    }
}