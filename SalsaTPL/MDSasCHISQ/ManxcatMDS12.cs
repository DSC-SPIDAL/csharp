using System;
using System.IO;
using System.Text;
using Manxcat;
using SALSALibrary;
using LMatrix = MathNet.Numerics.LinearAlgebra.Matrix;
using LU = MathNet.Numerics.LinearAlgebra.LUDecomposition;

namespace MDS
{
    public class ManxcatMDSDataProcessing
    {
        public static void UpdateManxcatMDS_Option12()
        {
            if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("traditionaldirectory"))
            {
                UpdateManxcatMDS_Option12_TraditionalDirectory();
                return;
            }

            if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("cluster"))
            {
                UpdateManxcatMDS_Option12_Cluster();
                return;
            }

            if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("cutfile"))
            {
                UpdateManxcatMDS_Option12_FamilybyCuts();
                return;
            }

            if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("function"))
            {
                if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("statistics"))
                    UpdateManxcatMDS_Option12_Functions(true);
                else
                    UpdateManxcatMDS_Option12_Functions(false);
                return;
            }
        }

        // End  UpdateManxcatMDS_Option12()

        public static void UpdateManxcatMDS_Option12_Functions(bool statistics)
        {
            //  Setup weight: Only used in statistics
            var weights = new double[SALSAUtility.PointCount_Global];
            ManxcatMDS.WeightingOption = ManxcatCentral.Configuration.WeightingOption;
            ManxcatMDS.SetupWeightings(weights);
            string TypicalMDSFileName = "";

            //  Setup MDS: Only used in statistics
            var MDSPoints = new double[SALSAUtility.PointCount_Global,Hotsun.ParameterVectorDimension];

            if (statistics)
            {
                // Set up statistics by reading MDS file
                TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" +
                                     ManxcatCentral.Configuration.ReducedVectorOutputFileName;

                if (!File.Exists(TypicalMDSFileName))
                {
                    Exception e = SALSAUtility.SALSAError(" File " + TypicalMDSFileName + " Does Not Exist");

                    throw (e);
                }

                int NumberMDSPoints = 0;
                var InitialMDSString = new string[SALSAUtility.PointCount_Global];
                var ColorValue = new int[SALSAUtility.PointCount_Global];

                if (ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue, ref NumberMDSPoints))
                    SALSAUtility.SALSAPrint(0, "MDS File " + TypicalMDSFileName + " Read Successfully");

                for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
                {
                    // Extract MDS POints
                    string inputLineStr = InitialMDSString[PointIndex].Trim();
                    string[] split = inputLineStr.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

                    if (split.Length < Hotsun.ParameterVectorDimension)
                    {
                        Exception e =
                            SALSAUtility.SALSAError(" Point " + PointIndex.ToString() + " has wrong number of Entries "
                                                    + inputLineStr + " Entries" + split.Length);

                        throw (e);
                    }

                    for (int VectorIndex = 0; VectorIndex < Hotsun.ParameterVectorDimension; VectorIndex++)
                    {
                        MDSPoints[PointIndex, VectorIndex] = Convert.ToDouble(split[VectorIndex].Trim());
                    }
                }
            } // End set up of Statistics

            //  Set Mapping
            var ClusterLabel = new int[SALSAUtility.PointCount_Global];
            int NumDivisions = 10;
            int BasicSize = SALSAUtility.PointCount_Global/NumDivisions;
            int LeftOver = SALSAUtility.PointCount_Global - BasicSize*NumDivisions;
            int DivisionCount = 0;
            int Limit = 0;

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                if (PointIndex >= Limit)
                {
                    ++DivisionCount;
                    Limit += BasicSize;

                    if (LeftOver >= DivisionCount)
                        ++Limit;
                }
                ClusterLabel[PointIndex] = DivisionCount;
            }

            String FunctionFileName = ManxcatCentral.Configuration.ConversionInformation;

            if (!FunctionFileName.Contains(":"))
                FunctionFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + FunctionFileName;

            var functionkeys = new double[SALSAUtility.PointCount_Global];
            var functionlabels = new int[SALSAUtility.PointCount_Global];

            int pickout = -1;

            while (true)
            {
                ++pickout;
                bool endofdata = true;
                int NumberofLines = 0;
                double sumabs = 0.0;
                double totalweight = 0.0;
                double meangamma = 0.0;
                double vargamma = 0.0;
                var meanxyz = new double[Hotsun.ParameterVectorDimension];
                var varxyz = new double[Hotsun.ParameterVectorDimension,Hotsun.ParameterVectorDimension];
                var correlxyzgamma = new double[Hotsun.ParameterVectorDimension];

                if (statistics)
                {
                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        meanxyz[VectorIndex1] = 0.0;

                        for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
                        {
                            varxyz[VectorIndex1, VectorIndex2] = 0.0;
                        }
                        correlxyzgamma[VectorIndex1] = 0.0;
                    }
                }

                try
                {
                    // Check if file exists
                    if (!File.Exists(FunctionFileName))
                    {
                        Exception e = SALSAUtility.SALSAError("File " + FunctionFileName + " does not exists.");

                        throw (e);
                    }

                    // Create a new stream to read from a file
                    using (StreamReader sr = File.OpenText(FunctionFileName))
                    {
                        // Read contents of a file, line by line, into a string
                        String inputLineStr;

                        while ((inputLineStr = sr.ReadLine()) != null)
                        {
                            inputLineStr = inputLineStr.Trim();
                            string[] split = inputLineStr.Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);

                            if (split.Length <= pickout)
                                break;
                            endofdata = false;
                            string usethis = split[pickout].Trim();

                            if (usethis.Length < 1)
                                continue; //replace empty line
                            double gamma = Convert.ToDouble(usethis);
                            functionkeys[NumberofLines] = gamma;
                            functionlabels[NumberofLines] = NumberofLines;
                            sumabs += Math.Abs(gamma);

                            if (statistics)
                            {
                                double wgt = weights[NumberofLines];
                                meangamma += gamma*wgt;
                                vargamma += wgt*gamma*gamma;
                                totalweight += wgt;

                                for (int VectorIndex1 = 0;
                                     VectorIndex1 < Hotsun.ParameterVectorDimension;
                                     VectorIndex1++)
                                {
                                    meanxyz[VectorIndex1] += MDSPoints[NumberofLines, VectorIndex1]*wgt;
                                    correlxyzgamma[VectorIndex1] += MDSPoints[NumberofLines, VectorIndex1]*gamma*wgt;

                                    for (int VectorIndex2 = 0;
                                         VectorIndex2 < Hotsun.ParameterVectorDimension;
                                         VectorIndex2++)
                                    {
                                        varxyz[VectorIndex1, VectorIndex2] += MDSPoints[NumberofLines, VectorIndex1]*
                                                                              MDSPoints[NumberofLines, VectorIndex2]*wgt;
                                    }
                                }
                            }
                            ++NumberofLines;
                        }
                        sr.Close();
                    }
                }
                catch (Exception e)
                {
                    SALSAUtility.SALSAError("Failed to read data from " + FunctionFileName + " " + e);

                    throw (e);
                }

                if (endofdata)
                    break;

                if (NumberofLines != SALSAUtility.PointCount_Global)
                {
                    SALSAUtility.SALSAError("Incorrect Function count read "
                                            + NumberofLines + " Expected " + SALSAUtility.PointCount_Global);
                }

                // Set Statistics
                if (statistics)
                {
                    var alpha = new double[Hotsun.ParameterVectorDimension];
                    double alphanorm = 0.0;
                    meangamma = meangamma/totalweight;
                    vargamma = (vargamma/totalweight) - meangamma*meangamma;

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        meanxyz[VectorIndex1] = meanxyz[VectorIndex1]/totalweight;
                    }

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
                        {
                            varxyz[VectorIndex1, VectorIndex2] = (varxyz[VectorIndex1, VectorIndex2]/totalweight) -
                                                                 meanxyz[VectorIndex1]*meanxyz[VectorIndex2];
                        }
                        correlxyzgamma[VectorIndex1] = correlxyzgamma[VectorIndex1]/totalweight -
                                                       meangamma*meanxyz[VectorIndex1];
                        alpha[VectorIndex1] = correlxyzgamma[VectorIndex1]/varxyz[VectorIndex1, VectorIndex1];
                        alphanorm += alpha[VectorIndex1]*alpha[VectorIndex1];
                    }

                    // invert Matrix to find best alpha
                    var ConventionalFirst = new double[Hotsun.ParameterVectorDimension,1];
                    var ConventionalAnswer = new double[Hotsun.ParameterVectorDimension][];
                    var ConventionalMatrix = new double[Hotsun.ParameterVectorDimension,Hotsun.ParameterVectorDimension];

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        ConventionalAnswer[VectorIndex1] = new double[1];
                        ConventionalFirst[VectorIndex1, 0] = correlxyzgamma[VectorIndex1]/
                                                             Math.Sqrt(varxyz[VectorIndex1, VectorIndex1]);

                        for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
                        {
                            ConventionalMatrix[VectorIndex1, VectorIndex2] = varxyz[VectorIndex1, VectorIndex2]/
                                                                             Math.Sqrt(
                                                                                 varxyz[VectorIndex1, VectorIndex1]*
                                                                                 varxyz[VectorIndex2, VectorIndex2]);
                        }
                    }
                    LMatrix cMatrix = LMatrix.Create(ConventionalMatrix);
                    LMatrix RightHandSide = LMatrix.Create(ConventionalFirst);
                    LMatrix MatrixAnswer = cMatrix.Solve(RightHandSide);
                    ConventionalAnswer = MatrixAnswer;

                    alphanorm = 0.0;

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        alpha[VectorIndex1] = ConventionalAnswer[VectorIndex1][0]/
                                              Math.Sqrt(varxyz[VectorIndex1, VectorIndex1]);
                        alphanorm += alpha[VectorIndex1]*alpha[VectorIndex1];
                    }

                    double Fullcorrelation = 0.0;
                    double varalphaxyz = 0.0;

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        alpha[VectorIndex1] = alpha[VectorIndex1]/Math.Sqrt(alphanorm);
                        Fullcorrelation += alpha[VectorIndex1]*correlxyzgamma[VectorIndex1];
                    }

                    for (int VectorIndex1 = 0; VectorIndex1 < Hotsun.ParameterVectorDimension; VectorIndex1++)
                    {
                        for (int VectorIndex2 = 0; VectorIndex2 < Hotsun.ParameterVectorDimension; VectorIndex2++)
                        {
                            varalphaxyz += alpha[VectorIndex1]*alpha[VectorIndex2]*varxyz[VectorIndex1, VectorIndex2];
                        }
                    }
                    Console.WriteLine(varalphaxyz + " " + Fullcorrelation + " " + vargamma);
                    Fullcorrelation = Fullcorrelation/(Math.Sqrt(vargamma*varalphaxyz));
                    SALSAUtility.SALSAPrint(0, "Column "
                                               + pickout.ToString() + " File " + FunctionFileName + " MDS File " +
                                               TypicalMDSFileName +
                                               " Total Weight "
                                               + totalweight.ToString() + " Correlation " +
                                               Fullcorrelation.ToString("F4") + " Direction "
                                               + alpha[0].ToString("F4") + " " + alpha[1].ToString("F4") + " " +
                                               alpha[2].ToString("F4"));
                }

                // Set Divisions
                Array.Sort(functionkeys, functionlabels);
                var PointColors = new int[SALSAUtility.PointCount_Global];

                for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
                {
                    int UnsortedPosition = functionlabels[PointIndex];
                    PointColors[UnsortedPosition] = ClusterLabel[PointIndex];
                }
                string OutputFileExtension = ManxcatCentral.Configuration.ConversionInformation;
                OutputFileExtension = OutputFileExtension.Replace(".dat", "");
                OutputFileExtension = OutputFileExtension.Replace(".txt", "");
                OutputFileExtension = OutputFileExtension + "-" + pickout.ToString() + ".txt";
                string labelFileDirectory = ManxcatCentral.Configuration.ClusterDirectory;
                    // Directory for Colors attached to Points
                string OutputFileName = labelFileDirectory + "\\" + OutputFileExtension;
                WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);

                double test = 0.1*sumabs/NumberofLines;

                for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
                {
                    int UnsortedPosition = functionlabels[PointIndex];

                    if (Math.Abs(functionkeys[UnsortedPosition]) < test)
                        PointColors[UnsortedPosition] = 0;
                }

                OutputFileName = OutputFileName.Replace(".txt", "NZ.txt");
                WritePointCluster(OutputFileName, PointColors, SALSAUtility.PointCount_Global);
            } // End while over pickout (columns in data)
            return;
        }

        // End UpdateManxcatMDS_Option12_Functions

        public static void UpdateManxcatMDS_Option12_FamilybyCuts()
        {
            String CutFileName = ManxcatCentral.Configuration.ConversionInformation;

            if (!CutFileName.Contains(":"))
                CutFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + CutFileName;

            var inputdata = new string[SALSAUtility.PointCount_Global];
            int NumberofCuts = 0;

            string FamilyLabel = "Family";
            int FamilyNumber = 1;

            if (ManxcatCentral.Configuration.ConversionOption.ToLower().Contains("family2"))
                FamilyNumber = 2;

            ReadString_File(CutFileName, inputdata, ref NumberofCuts);

            if (NumberofCuts <= 0)
            {
                Exception e = SALSAUtility.SALSAError("Too few points "
                                                      + NumberofCuts.ToString() + " in Cut File " + CutFileName);

                throw (e);
            }

            var CutStart = new int[NumberofCuts];
            var CutLabel = new string[NumberofCuts];
            int CountFamilies = 0;
            int itmpold = -2;

            for (int CutPos = 0; CutPos < NumberofCuts; CutPos++)
            {
                string[] split = inputdata[CutPos].Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                int itmp = Convert.ToInt32(split[0]);

                if (itmp < 0)
                {
                    FamilyLabel = split[1];
                    continue;
                }

                if (itmp <= itmpold)
                    continue;
                CutStart[CountFamilies] = itmp;
                CutLabel[CountFamilies] = FamilyLabel + "-" + split[1];
                ++CountFamilies;
                itmpold = itmp;
            }

            if (FamilyNumber == 1)
                SALSAUtility.GlobalFileProperties.FamilyName1 = FamilyLabel + " from " + CutFileName;

            if (FamilyNumber == 2)
                SALSAUtility.GlobalFileProperties.FamilyName2 = FamilyLabel + " from " + CutFileName;

            int minlabel = 1;
            int CurrentFamily = 0;
            var CategoryIndex = new int[SALSAUtility.PointCount_Global];

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                if (CurrentFamily < CountFamilies - 1)
                {
                    if ((PointIndex + 1) >= CutStart[CurrentFamily + 1])
                    {
                        ++CurrentFamily;
                    }
                }
                CategoryIndex[PointIndex] = minlabel + CurrentFamily;
                continue;
            }

            var FamilyCounts = new int[CountFamilies];

            for (int Icount = 0; Icount < CountFamilies; Icount++)
            {
                FamilyCounts[Icount] = 0;
            }

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                ++FamilyCounts[CategoryIndex[PointIndex] - minlabel];
            }

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                int Icount = CategoryIndex[PointIndex];

                if (FamilyNumber == 1)
                {
                    SALSAUtility.GlobalPointProperties[PointIndex].family1 = Icount;
                    SALSAUtility.GlobalPointProperties[PointIndex].familylabel1 = CutLabel[Icount - minlabel] + " "
                                                                                  + Icount.ToString() + " out of "
                                                                                  + CountFamilies.ToString() + " Count " +
                                                                                  FamilyCounts[Icount - minlabel].
                                                                                      ToString();
                }

                if (FamilyNumber == 2)
                {
                    SALSAUtility.GlobalPointProperties[PointIndex].family2 = Icount;
                    SALSAUtility.GlobalPointProperties[PointIndex].familylabel2 = CutLabel[Icount - minlabel] + " "
                                                                                  + Icount.ToString() + " out of "
                                                                                  + CountFamilies.ToString() + " Count " +
                                                                                  FamilyCounts[Icount - minlabel].
                                                                                      ToString();
                }
            }

            string cleandate = SALSAUtility.startTime.ToLocalTime().ToString();
            cleandate = cleandate.Replace(":", ".");
            cleandate = cleandate.Replace(" ", "-");
            string OldComment = SALSAUtility.GlobalFileProperties.Comment;

            if (OldComment != "")
                OldComment += "\n";
            OldComment += "Family"
                          + FamilyNumber.ToString() + " " + FamilyLabel + " Information added from " + CutFileName +
                          " Time " + cleandate;
            SALSAUtility.GlobalFileProperties.Comment = OldComment;

            //  Write new label file
            string labelfilename = ManxcatCentral.Configuration.DataLabelsFileName;

            if (!labelfilename.Contains(":"))
                labelfilename = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + labelfilename;
            SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.Configuration.Write2Das3D,
                                                "colon,SameFileName ", SALSAUtility.GlobalFileProperties,
                                                SALSAUtility.GlobalPointProperties, SALSAUtility.NumberOriginalPoints);
            SALSAUtility.SALSAPrint(0, "Family"
                                       + FamilyNumber.ToString() + " " + FamilyLabel + " Info Added to " + labelfilename +
                                       " from " + CutFileName);
            ManxcatCentral.Configuration.Comment += "\nFamily"
                                                    + FamilyNumber.ToString() + " " + FamilyLabel + " Info Added to " +
                                                    labelfilename + " from " + CutFileName;
        }

        // End UpdateManxcatMDS_Option12_FamilybyCuts

        public static void UpdateManxcatMDS_Option12_Cluster()
        {
            //  Read Cluster Information with File specified in ConversionInformation
            string ClusterFileName = ManxcatCentral.Configuration.ConversionInformation;
            int NumberOfClusterLines = 0;
            var CategoryLabel = new string[SALSAUtility.PointCount_Global];
            ReadClusterLabel_File(ClusterFileName, CategoryLabel, ref NumberOfClusterLines);

            if (NumberOfClusterLines != SALSAUtility.PointCount_Global)
            {
                Exception e = SALSAUtility.SALSAError(" Illegal Count "
                                                      + NumberOfClusterLines.ToString() + " in File " + ClusterFileName);

                throw (e);
            }

            var CategoryIndex = new int[SALSAUtility.PointCount_Global];
            int minlabel = SALSAUtility.PointCount_Global;
            int maxlabel = 0;

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                CategoryIndex[PointIndex] = Convert.ToInt32(CategoryLabel[PointIndex]);

                if (minlabel > CategoryIndex[PointIndex])
                    minlabel = CategoryIndex[PointIndex];

                if (maxlabel < CategoryIndex[PointIndex])
                    maxlabel = CategoryIndex[PointIndex];
            }

            int TotalNumberClusters = maxlabel - minlabel + 1;
            var ClusterCounts = new int[TotalNumberClusters];

            for (int Icount = 0; Icount < TotalNumberClusters; Icount++)
            {
                ClusterCounts[Icount] = 0;
            }

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                ++ClusterCounts[CategoryIndex[PointIndex] - minlabel];
            }

            for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
            {
                int Icount = CategoryIndex[PointIndex];
                SALSAUtility.GlobalPointProperties[PointIndex].cluster = Icount;
                SALSAUtility.GlobalPointProperties[PointIndex].clusterlabel = "Cluster "
                                                                              + Icount.ToString() + " out of "
                                                                              + TotalNumberClusters.ToString() +
                                                                              " Count " +
                                                                              ClusterCounts[Icount - minlabel].ToString();
            }
            SALSAUtility.GlobalFileProperties.ClusterName = "Clusters from " + ClusterFileName;
            SALSAUtility.GlobalFileProperties.ClusterStartIndex = minlabel;

            string cleandate = SALSAUtility.startTime.ToLocalTime().ToString();
            cleandate = cleandate.Replace(":", ".");
            cleandate = cleandate.Replace(" ", "-");
            string OldComment = SALSAUtility.GlobalFileProperties.Comment;

            if (OldComment != "")
                OldComment += "\n";
            OldComment += "Cluster Information added from " + ClusterFileName + " Time " + cleandate;
            SALSAUtility.GlobalFileProperties.Comment = OldComment;

            //  Write new label file
            string labelfilename = ManxcatCentral.Configuration.DataLabelsFileName;

            if (!labelfilename.Contains(":") && !labelfilename.Contains("$"))
                labelfilename = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + labelfilename;
            SALSA_Properties.WriteDataPointFile(labelfilename, ManxcatCentral.Configuration.Write2Das3D,
                                                "colon,SameFileName ", SALSAUtility.GlobalFileProperties,
                                                SALSAUtility.GlobalPointProperties, SALSAUtility.NumberOriginalPoints);
            SALSAUtility.SALSAPrint(0, "Cluster Info Added to " + labelfilename + " from " + ClusterFileName);
            ManxcatCentral.Configuration.Comment += "\nCluster Info Added to "
                                                    + labelfilename + " from " + ClusterFileName;
        }

        // End UpdateManxcatMDS_Option12_Cluster()

        public static void UpdateManxcatMDS_Option12_TraditionalDirectory()
        {
            string labelFileDirectory = ManxcatCentral.Configuration.ClusterDirectory; // Input Directory
            string MDSandClusterDirectory = labelFileDirectory + "\\" + ManxcatCentral.Configuration.RunSetLabel +
                                            "-R" + ManxcatCentral.Configuration.RunNumber.ToString() + "-ManxcatMDS";
                // Output Directory

            if (Directory.Exists(MDSandClusterDirectory))
            {
                SALSAUtility.SALSAPrint(0, "The directory " + MDSandClusterDirectory + " exists");
            }
            else
            {
                DirectoryInfo di = Directory.CreateDirectory(MDSandClusterDirectory);
                SALSAUtility.SALSAPrint(0, "The directory " + MDSandClusterDirectory + " was created successfully at "
                                           + Directory.GetCreationTime(MDSandClusterDirectory));
            }

            string[] SMACOFfileEntries = Directory.GetFiles(MDSandClusterDirectory);
            string TypicalMDSFileName = ManxcatCentral.Configuration.ReducedVectorOutputFileName;

            if (!TypicalMDSFileName.Contains(":"))
                TypicalMDSFileName = ManxcatCentral.ResultDirectoryName + "\\SIMPLE" + TypicalMDSFileName;

            if (!File.Exists(TypicalMDSFileName))
            {
                Exception e = SALSAUtility.SALSAError(" File " + TypicalMDSFileName + " Does Not Exist");

                throw (e);
            }

            int NumberMDSPoints = 0;
            var InitialMDSString = new string[SALSAUtility.PointCount_Global];
            var ColorValue = new int[SALSAUtility.PointCount_Global];

            if (ReadMDSCluster_File(TypicalMDSFileName, InitialMDSString, ColorValue, ref NumberMDSPoints))
                Console.WriteLine(TypicalMDSFileName + " Read Successfully");

            string[] fileEntries = Directory.GetFiles(labelFileDirectory);

            foreach (string LabelFileName in fileEntries)
            {
                if (LabelFileName.Contains(".dot"))
                    continue;

                if (LabelFileName.Contains("Summary"))
                    continue;

                if (LabelFileName.Contains("Timing"))
                    continue;
                string LabelFileName1 = LabelFileName.Replace(labelFileDirectory + "\\", "");
                string coreFileName = MDSandClusterDirectory + "\\" + "MDSManxcat-" + LabelFileName1;

                if (File.Exists(coreFileName))
                    continue;
                int NumberOfLabels = 0;
                var CategoryLabel = new string[SALSAUtility.PointCount_Global];
                ReadClusterLabel_File(LabelFileName, CategoryLabel, ref NumberOfLabels);

                if (NumberOfLabels != SALSAUtility.PointCount_Global)
                {
                    Exception e = SALSAUtility.SALSAError(" Illegal Count "
                                                          + NumberOfLabels.ToString() + " in File " + LabelFileName);

                    throw (e);
                }

                var CategoryIndex = new int[SALSAUtility.PointCount_Global];

                for (int PointIndex = 0; PointIndex < SALSAUtility.PointCount_Global; PointIndex++)
                {
                    CategoryIndex[PointIndex] = Convert.ToInt32(CategoryLabel[PointIndex]);
                }
                WriteColor_Cluster(coreFileName, InitialMDSString, CategoryIndex, SALSAUtility.PointCount_Global, false);
                SALSAUtility.SALSAPrint(0, "Traditional Directory MDS Info Added to " + coreFileName);
            }
        }

        // End UpdateManxcatMDS_Option12_TraditionalDirectory

        public static bool ReadMDSCluster_File(string MDSClusterFileName, String[] InitialString, int[] ColorValue,
                                               ref int NumberofPoints)
        {
            NumberofPoints = 0;

            try
            {
                // Check if file exists
                if (!File.Exists(MDSClusterFileName))
                {
                    Exception e = SALSAUtility.SALSAError("File " + MDSClusterFileName + " does not exists");

                    throw (e);
                }

                // Create a new stream to read from a file
                using (StreamReader sr = File.OpenText(MDSClusterFileName))
                {
                    // Read contents of a file, line by line, into a string
                    String inputLineStr;
                    int newlabelnumber = -1;

                    while ((inputLineStr = sr.ReadLine()) != null)
                    {
                        if (inputLineStr.Length < 2)
                            continue; //replace empty line
                        inputLineStr = inputLineStr.Trim(new[] {' ', '\t'});

                        try
                        {
                            // Parse each record string
                            inputLineStr = inputLineStr.Replace("\t\t", "\t");
                            string[] split = inputLineStr.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                            newlabelnumber = Convert.ToInt32(split[0]);

                            if (split.Length != Hotsun.ParameterVectorDimension + 2)
                                Console.WriteLine(" Bad Line "
                                                  + split.Length + " " + NumberofPoints + " " + inputLineStr);

                            // Begin Changes saliya 3/22/11
                            // Note. This seems unnecessary. Let's agree on zero based indices
//                            if ((NumberofPoints + 1) != newlabelnumber)
//                            {
//                                Exception e = SALSAUtility.SALSAError("Unexpected Label Number "
//                                    + newlabelnumber + " Expected " + NumberofPoints + "+1");
//
//                                throw (e);
//                            }
                            // End Changes saliya 3/22/11

                            inputLineStr = split[0];

                            for (int i = 1; i < (split.Length - 1); i++)
                            {
                                inputLineStr += "\t" + split[i];
                            }
                            InitialString[NumberofPoints] = inputLineStr;

                            ColorValue[NumberofPoints] = (int) Convert.ToDouble(split[split.Length - 1]);
                            ++NumberofPoints;
                        }
                        catch (Exception e)
                        {
                            SALSAUtility.SALSAError("Failed to load data array " + inputLineStr + " "
                                                    + " "
                                                    + NumberofPoints.ToString() + " " + newlabelnumber.ToString() + " " +
                                                    e);

                            throw (e);
                        }
                    }

                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed to read data from " + MDSClusterFileName + " " + e);

                throw (e);
            }

            if (NumberofPoints != SALSAUtility.PointCount_Global)
            {
                Exception e = SALSAUtility.SALSAError("Incorrect #Points in File "
                                                      + NumberofPoints + " Expected " + SALSAUtility.PointCount_Global);

                throw (e);
            }

            return true;
        }

        // End ReadMDSClusterFile

        public static bool ReadClusterLabel_File(string LabelFileName, String[] LabelArray, ref int NumberofLabels)
        {
            NumberofLabels = 0;
            int startnumber = 1;

            try
            {
                // Check if file exists
                if (!File.Exists(LabelFileName))
                {
                    Exception e = SALSAUtility.SALSAError("File " + LabelFileName + " does not exists.");

                    throw (e);
                }

                // Create a new stream to read from a file
                using (StreamReader sr = File.OpenText(LabelFileName))
                {
                    // Read contents of a file, line by line, into a string
                    String inputLineStr;

                    while ((inputLineStr = sr.ReadLine()) != null)
                    {
                        if (inputLineStr.Length < 2)
                            continue; //replace empty line

                        inputLineStr = inputLineStr.Trim();

                        try
                        {
                            // Parse each record string
                            string[] split = inputLineStr.Split(new[] {' ', '\t'}, 2,
                                                                StringSplitOptions.RemoveEmptyEntries);
                            int newlabelnumber = Convert.ToInt32(split[0]);

                            if (NumberofLabels == 0)
                            {
                                startnumber = newlabelnumber;

                                if ((startnumber < 0) || (startnumber > 1))
                                {
                                    Exception e =
                                        SALSAUtility.SALSAError("Unexpected Start Number " + startnumber.ToString());

                                    throw (e);
                                }
                            }

                            if (NumberofLabels != (newlabelnumber + startnumber))
                            {
                                Exception e = SALSAUtility.SALSAError("Unexpected Label Number "
                                                                      + newlabelnumber.ToString() + " Expected " +
                                                                      NumberofLabels.ToString() + " + " +
                                                                      startnumber.ToString());

                                throw (e);
                            }

                            if (split[1].Length <= 0)
                            {
                                Exception e =
                                    SALSAUtility.SALSAError("Zero length label for point " + NumberofLabels.ToString());

                                throw (e);
                            }

                            LabelArray[NumberofLabels] = split[1];
                            ++NumberofLabels;
                        }
                        catch (Exception e)
                        {
                            SALSAUtility.SALSAError("Failed to load data from " + LabelFileName + " " + e);

                            throw (e);
                        }
                    }

                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed to read data from " + LabelFileName + " " + e);

                throw (e);
            }

            return true;
        }

        // End ReadClusterLabel_File


        // Write label-cluster results into a file
        public static void WriteColor_Cluster(string fname, string[] labels, int[] ColorValues, int dataPoints,
                                              bool append)
        {
            try
            {
                StreamWriter sw = null;

                if (!string.IsNullOrEmpty(fname))
                {
                    sw = new StreamWriter(fname, append, Encoding.UTF8);
                }

                if (sw != null)
                {
                    for (int i = 0; i < dataPoints; i++)
                    {
                        string stripped = labels[i].Trim(new[] {' ', '\t'});
                        sw.WriteLine(String.Format(stripped + "\t" + ColorValues[i].ToString() + ".0000000000"));
                    }
                }

                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed writing data on " + fname + " " + e);

                throw (e);
            }
        }

        // End WriteColor_Cluster

        // Write Point-cluster results into a file
        public static void WritePointCluster(string fname, int[] ColorValues, int dataPoints)
        {
            try
            {
                StreamWriter sw = null;

                if (!string.IsNullOrEmpty(fname))
                {
                    sw = new StreamWriter(fname, false, Encoding.UTF8);
                }

                if (sw != null)
                {
                    for (int i = 0; i < dataPoints; i++)
                    {
                        sw.WriteLine(String.Format((i + 1).ToString() + " " + ColorValues[i].ToString()));
                    }
                }
                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed writing data on " + fname + " " + e);

                throw (e);
            }
        }

        // End WritePointCluster

        public static bool ReadString_File(string StringFileName, String[] LabelArray, ref int NumberofLines)
        {
            NumberofLines = 0;

            try
            {
                // Check if file exists
                if (!File.Exists(StringFileName))
                {
                    Exception e = SALSAUtility.SALSAError("File " + StringFileName + " does not exists.");

                    throw (e);
                }

                // Create a new stream to read from a file
                using (StreamReader sr = File.OpenText(StringFileName))
                {
                    // Read contents of a file, line by line, into a string
                    String inputLineStr;

                    while ((inputLineStr = sr.ReadLine()) != null)
                    {
                        if (inputLineStr.Length < 2)
                            continue; //replace empty line

                        LabelArray[NumberofLines] = inputLineStr.Trim();
                        ++NumberofLines;
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed to read data from " + StringFileName + " " + e);

                throw (e);
            }

            return true;
        }

        // End ReadString_File
    }
}