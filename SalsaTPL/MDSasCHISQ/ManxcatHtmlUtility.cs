using System;
using System.IO;
using System.Reflection;
using System.Text;
using MDS;
using SALSALibrary;
using Salsa.Core.Configuration.Sections;

namespace Manxcat
{
    public class ManxcatHtmlUtility
    {
        /// <summary>
        /// Generates a Web page linking to output files
        /// </summary>
        public static void WriteHTML()
        {
            ManxcatSection config = ManxcatCentral.Configuration;
            string directoryPath = Path.GetDirectoryName(config.ReducedVectorOutputFileName);
            directoryPath = string.IsNullOrEmpty(directoryPath) ? string.Empty : directoryPath;

            string stylesFile = Path.Combine(directoryPath, ManxcatConstants.StyleFileNameTag) + ".css";
            CopyStylesCSS(stylesFile);

            string htmlFile = Path.Combine(directoryPath, ManxcatConstants.HtmlFileNameTag) + ".html";
            GenerateHtmlContent(htmlFile, config);
        }

        private static void GenerateHtmlContent(string htmlFile, ManxcatSection config)
        {
            using (var writer = new StreamWriter(htmlFile))
            {
                using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.template")
                    )
                {
                    if (stream != null)
                    {
                        using (var reader = new StreamReader(stream))
                        {
                            while (!reader.EndOfStream)
                            {
                                string line = reader.ReadLine();
                                if (!string.IsNullOrEmpty(line))
                                {
                                    writer.WriteLine(!line.Contains(ManxcatConstants.Tag)
                                                         ? line
                                                         : ProcessTagLine(line, config));
                                }
                            }
                        }
                    }
                    else
                    {
                        throw new Exception(ManxcatErrorMessages.TemplateHtmlNotAvailable);
                    }
                }
            }
        }

        private static string ProcessTagLine(string line, ManxcatSection config)
        {
            int tagStartIdx = line.IndexOf(ManxcatConstants.Tag);
            int tagEndIdx = line.IndexOf(ManxcatConstants.Tag, tagStartIdx + ManxcatConstants.Tag.Length) +
                            (ManxcatConstants.Tag.Length - 1);
            string prefix = line.Substring(0, tagStartIdx);
            string suffix = line.Substring(tagEndIdx + 1);
            string meat = line.Substring(tagStartIdx + ManxcatConstants.Tag.Length,
                                         ((tagEndIdx + 1) - (tagStartIdx + 2*ManxcatConstants.Tag.Length)));
            switch (meat)
            {
                case ManxcatConstants.TagManxcatName:
                    return prefix + config.ManxcatRunName + suffix;
                case ManxcatConstants.TagManxcatDescription:
                    return prefix + config.ManxcatRunDescription + suffix;
                case ManxcatConstants.TagManxcatConfig:
                    return GenerateConfig(config);
                case ManxcatConstants.TagManxcatLinks:
                    return prefix + GenerateLinks(config) + suffix;
                default:
                    throw new Exception(ManxcatErrorMessages.UnidentifiedTagInHTMLGeneration);
            }
        }

        private static string GenerateConfig(ManxcatSection config)
        {
            var sb = new StringBuilder();
            sb.AppendLine("I/O");
            sb.AppendLine("\tCoordinateWriteFrequency:      " + config.CoordinateWriteFrequency);
            sb.AppendLine("\tDistanceMatrixFile:            " + config.DistanceMatrixFile);
            sb.AppendLine("\n\nManxcatCore");
            sb.AppendLine("\tAddonforQcomputation:          " + config.AddonforQcomputation);
            sb.AppendLine("\tCalcFixedCrossFixed:           " + config.CalcFixedCrossFixed);
            sb.AppendLine("\tCGResidualLimit:               " + config.CGResidualLimit);
            sb.AppendLine("\tChisqChangePerPoint:           " + config.ChisqChangePerPoint);
            sb.AppendLine("\tChisqnorm:                     " + config.Chisqnorm);
            sb.AppendLine("\tChisqPrintConstant:            " + config.ChisqPrintConstant);
            sb.AppendLine("\tConversionInformation:         " + config.ConversionInformation);
            sb.AppendLine("\tConversionOption:              " + config.ConversionOption);
            sb.AppendLine("\tDataPoints:                    " + config.DataPoints);
            sb.AppendLine("\tDerivtest:                     " + config.Derivtest);
            sb.AppendLine("\tDiskDistanceOption:            " + config.DiskDistanceOption);
            sb.AppendLine("\tDistanceCut:                   " + config.DistanceCut);
            sb.AppendLine("\tDistanceFormula:               " + config.DistanceFormula);
            sb.AppendLine("\tDistanceProcessingOption:      " + config.DistanceProcessingOption);
            sb.AppendLine("\tDistanceWeigthsCuts:           " + config.DistanceWeightsCuts);
            sb.AppendLine("\tEigenvaluechange:              " + config.Eigenvaluechange);
            sb.AppendLine("\tEigenvectorchange:             " + config.Eigenvectorchange);
            sb.AppendLine("\tExtraOption1:                  " + config.ExtraOption1);
            sb.AppendLine("\tExtraprecision:                " + config.Extraprecision);
            sb.AppendLine("\tFixedPointCriterion:           " + config.FixedPointCriterion);
            sb.AppendLine("\tFletcherRho:                   " + config.FletcherRho);
            sb.AppendLine("\tFletcherSigma:                 " + config.FletcherSigma);
            sb.AppendLine("\tFullSecondDerivativeOption:    " + config.FullSecondDerivativeOption);
            sb.AppendLine("\tFunctionErrorCalcMultiplier:   " + config.FunctionErrorCalcMultiplier);
            sb.AppendLine("\tHistogramBinCount              " + config.HistogramBinCount);
            sb.AppendLine("\tInitializationLoops:           " + config.InitializationLoops);
            sb.AppendLine("\tInitializationOption:          " + config.InitializationOption);
            sb.AppendLine("\tInitialSteepestDescents:       " + config.InitialSteepestDescents);
            sb.AppendLine("\tLinkCut:                       " + config.LinkCut);
            sb.AppendLine("\tLocalVectorDimension:          " + config.LocalVectorDimension);
            sb.AppendLine("\tMaxit:                         " + config.Maxit);
            sb.AppendLine("\tMinimumDistance:               " + config.MinimumDistance);
            sb.AppendLine("\tMPIIOStrategy:                 " + config.MPIIOStrategy);
            sb.AppendLine("\tNbadgo:                        " + config.Nbadgo);
            sb.AppendLine("\tOmega:                         " + config.Omega);
            sb.AppendLine("\tOmegaOption:                   " + config.OmegaOption);
            sb.AppendLine("\tPowerIterationLimit:           " + config.PowerIterationLimit);
            sb.AppendLine("\tProcessingOption:              " + config.ProcessingOption);
            sb.AppendLine("\tQgoodReductionFactor:          " + config.QgoodReductionFactor);
            sb.AppendLine("\tQHighInitialFactor:            " + config.QHighInitialFactor);
            sb.AppendLine("\tQLimiscalecalculationInterval: " + config.QLimitscalculationInterval);
            sb.AppendLine("\tRotationOption:                " + config.RotationOption);
            sb.AppendLine("\tSelectedfixedpoints:           " + config.Selectedfixedpoints);
            sb.AppendLine("\tSelectedvariedpoints:          " + config.Selectedvariedpoints);
            sb.AppendLine("\tStoredDistanceOption:          " + config.StoredDistanceOption);
            sb.AppendLine("\tTimeCutmillisec:               " + config.TimeCutmillisec);
            sb.AppendLine("\tTransformMethod:               " + config.TransformMethod);
            sb.AppendLine("\tTransformParameter:            " + config.TransformParameter);
            sb.AppendLine("\tUndefindDistanceValue:         " + config.UndefinedDistanceValue);
            sb.AppendLine("\tVariedPointCriterion:          " + config.VariedPointCriterion);
            sb.AppendLine("\tWeightingOption:               " + config.WeightingOption);
            sb.AppendLine("\tWrite2Das3D:                   " + config.Write2Das3D);
            sb.AppendLine("\n\nDensity");
            sb.AppendLine("\tAlpha:                         " + config.Alpha);
            sb.AppendLine("\tPcutf:                         " + config.Pcutf);
            sb.AppendLine("\tSelectedClusters:              " + config.SelectedClusters);
            sb.AppendLine("\tXmaxBound:                     " + config.XmaxBound);
            sb.AppendLine("\tXres:                          " + config.Xres);
            sb.AppendLine("\tYmaxBound:                     " + config.YmaxBound);
            sb.AppendLine("\tYres:                          " + config.Yres);
            return sb.ToString();
        }

        private static string GenerateLinks(ManxcatSection config)
        {
            const string liTemplate = "<li><a href=\"{0}\">{1}</a></li>";
            const string txtExt = ".txt";
            string dotslash = string.IsNullOrEmpty(config.ServerUrlPrefix)
                                  ? "./"
                                  : (config.ServerUrlPrefix + "/" + config.ManxcatRunName + "/");
            var sb = new StringBuilder("<ul>");

            string reducedVectorFileNameTag = Path.GetFileNameWithoutExtension(config.ReducedVectorOutputFileName);

            string link = dotslash + ManxcatConstants.SimplePointsPrefix + reducedVectorFileNameTag + txtExt;
            string name = "Simple Points";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + reducedVectorFileNameTag + ManxcatConstants.ColonPointsSuffix + txtExt;
            name = "Colon Points";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + reducedVectorFileNameTag + ManxcatConstants.GroupPointsSuffix + txtExt;
            name = "Group Points";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + Path.GetFileName(config.SummaryOutputFileName);
            name = "Summary File";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + Path.GetFileName(config.TimingOutputFileName);
            name = "Timing File";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + ManxcatConstants.OutFileNameTag + txtExt;
            name = "Manxcat Ouput";
            sb.AppendLine(string.Format(liTemplate, link, name));

            link = dotslash + "whole-plot.png";
            name = "Density Graph";
            sb.AppendLine(string.Format(liTemplate, link, name));

            if (SALSAUtility.IsClustersSelected)
            {
                link = dotslash + "selected-plot.png";
                name = "Density Graph for Selected Clusters (Intra Cluster Distances)";
                sb.AppendLine(string.Format(liTemplate, link, name));

                link = dotslash + "selected-inter-plot.png";
                name = "Density Graph for Selected Clusters (Inter Cluster Distances)";
                sb.AppendLine(string.Format(liTemplate, link, name));
            }

            return sb.ToString();
        }

        private static void CopyStylesCSS(string toPath)
        {
            using (var writer = new StreamWriter(toPath))
            {
                using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Manxcat.Web.style"))
                {
                    if (stream != null)
                    {
                        using (var reader = new StreamReader(stream))
                        {
                            while (!reader.EndOfStream)
                            {
                                writer.WriteLine(reader.ReadLine());
                            }
                        }
                    }
                    else
                    {
                        throw new Exception(ManxcatErrorMessages.TemplateStyleSheetNotAvailable);
                    }
                }
            }
        }

        public static void TestMethod()
        {
            string directoryName = Path.GetDirectoryName(ManxcatCentral.Configuration.ReducedVectorOutputFileName);
            string testfile = !string.IsNullOrEmpty(directoryName)
                                  ? Path.Combine(directoryName, "TestMethod_rank_" + SALSAUtility.MPI_Rank + ".txt")
                                  : "TestMethod_rank_" + SALSAUtility.MPI_Rank + ".txt";
            using (var writer = new StreamWriter(testfile))
            {
                for (int globalPointIndex = 0; globalPointIndex < SALSAUtility.PointCount_Global; globalPointIndex++)
                {
                    if (ManxcatMDS.PointStatus[globalPointIndex] == -1)
                        continue;
                    int originalPointIndex = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex];

                    string coordinates = "";
                    const int singleCluster = 1;
                    int usedPointIndex = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex];

                    /* for (int globalPointIndex2 = 0; globalPointIndex2 < SALSAUtility.PointCount_Global; globalPointIndex2++)
                    {
                        if (ManxcatMDS.PointStatus[globalPointIndex2] == -1)
                            continue;
                        int originalPointIndex2 = SALSAUtility.UsedPointtoOriginalPointMap[globalPointIndex2];
                        int usedPointIndex2 = SALSAUtility.NaivetoActualUsedOrder[globalPointIndex2];
                    }*/

                    for (int localVectorIndex = 0; localVectorIndex < 3; localVectorIndex++)
                    {
                        coordinates += Hotsun.GlobalParameter[usedPointIndex][localVectorIndex].ToString("E4") + "\t";
                    }

                    writer.WriteLine(String.Format(globalPointIndex + "\t" + coordinates + singleCluster));
                }
            }
        }
    }
}