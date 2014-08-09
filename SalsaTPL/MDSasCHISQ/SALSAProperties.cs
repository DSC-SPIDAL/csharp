using System;
using System.IO;
using System.Text;

//    OLD Comments on input Data
//    public string DataFileName = " "; // Input data file name
//    public string DataLabelsFileName = "";   // File with application specific labels for this data
//    public string ReducedVectorOutputFileName = "";  // File in ResultDirectoryName holding ReducedVector Output
//    public string ClusterDirectory = "";    // Directory to find Cluster Output
//    public int DataPoints = 0;  // Number of Original Data Points
//    public int ProcessingOption = 0; // Control Requested Processing
//    public int DistanceProcessingOption = 0; // Control Input data Processing ( = 2 SQUARE Input Values)
//    public int InitializationOption = 0;    // Control Initialization =0 None =1 Old Model =2 New File Type
//    public string InitializationFileName = "";  // File
//    public bool write2Das3D = true; // for 2D MDS write output for 3D viewer
//    public int LocalVectorDimension = 3;    // Final MDS dimension (2 or 3 expected)
//    public string selectedvariedpoints = "";  // List of varied points specified by point number or cluster number 
//    public string VariedPointCriterion = "all"; // If "all" or "rest" all the non-fixed points varied; if "allinfile" choose all in file and ignore selectedvariedpoints property
//    // if "originalpoint" select by Original Point; = "family1", "family2", "cluster", "group" select by this division
//    public string selectedvariedpointfile = "";  // File containing mapping of groups to original values 
//    public string selectedfixedpoints = "";  // List of fixed points specified by point number or group number 
//    public string FixedPointCriterion = "none"; // If "none" there are no fixed points; if "allinfile" choose all in file and ignore selectedfixedpoints property
//    // if "originalpoint" select by Original Point; = "family1", "family2", "cluster", "group" select by this division
//    public string selectedfixedpointfile = "";  // File containing coordinates of fixed points

//    public int RotationOption = 0;  // Control Rotation 
//    public string RotationLabelsFileName = "";   // File with application specific labels for data to be rotated

//    public int Chisqnorm = 0;   // Method of Normalizing Chisq
//    public int DistanceFormula = 0;   // Method of Calculating Distance
//    public int FullSecondDerivativeOption = 0; // Use Full second Derivative Calculation if > 0
//    public double MinimumDistance = -0.01;   // Minimum Distance for Chisq Normalization
//    public double FunctionErrorCalcMultiplier = 10.0; // Multiplicative Constant for Function (square this for Chisq) Error Calculation
//    public double ChisqPrintConstant = 1.0; // Multiplicative Constant for Chisq Printing
//    public int maxit = 60;  // Maximum Number of Manxcat Iterations
//    public int nbadgo = 6;  // Limit on Iterations with no progress
//    public double ChisqChangePerPoint = 0.001;  // Stop when this is Chisq change per point
//    public double FletcherRho = 0.25;   // Fletcher's parameter rho
//    public double FletcherSigma = 0.75; // Fletcher's Parameter sigma
//    public double Omega = 1.25; //  Parameter Omega deciding if to do Linear Search if OmegaOption = 1
//    public int OmegaOption = 0; // If nonzero do linear searches if promising
//    public double QHighInitialFactor = 0.01; // Initial call uses this fraction of Qhigh for Q (added 2009)
//    public double QgoodReductionFactor = 0.5;    // Set to zero to make Q = 0. for good solutions
//    public double addonforQcomputation = 2.0;   // Add this to make Power Method Safer
//    public double extraprecision = 0.05; //  Extra precision for Qlow
//    public int QLimitscalculationInterval = 1; //  recalculate Qhigh Qlow Trace and Norm after this number of iterations
//    public int InitialSteepestDescents = 0;  // Start minimization with this number of Steepest Descent Iterations
//    public double TimeCutmillisec = -1.0;  // if positive, Stop when time used in milliseconds greater than timmax
//    public bool derivtest = false; // If true, do a derivative 

//    public double CGResidualLimit = 0.00001;    // Stop when Norm of Residual is less than this
//    public int PowerIterationLimit = 200;    // Limit on Power Iterations
//    public double eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
//    public double eigenvectorchange = 0.001; // convergence test in power method for eigenvector

//    public string TimingOutputFileName = " "; // Timing Output file name
//    public string BaseResultDirectoryName = "";  // Base Directory holding all Runs
//    public string ResultDirectoryExtension = ""; // 
//    public string ControlDirectoryName = ""; // Directory holding Control Information for this RunSet
//    public string RunSetLabel = "";  //  Label for RunSet
//    public int RunNumber = 0;    // Unique number of run

//    public string pattern = ""; // Parallel Pattern
//    public int ThreadCount = 1;    // Number of Threads
//    public int NodeCount = 1;    // Number of Nodes
//    public int MPIperNodeCount = 1;    // Number of MPI Processes per Node
//    public int MPIIOStrategy = 0;   // Strategy for I/O in MPI

//    public int HistogramBinCount = 100; // Bin Count for Histograms
//    public string Comment = ""; // User Comment
//    public string Extradata1 = "";  // Use in special ways
//    public string Extradata2 = "";  // Use in special ways
//    public string Extradata3 = "";  // Use in special ways
//    public string Extradata4 = "";  // Use in special ways
//    public int ExtraOption1 = 0;    // use in various ways
//    public int DebugPrintOption = 1; // Control Debug Printing (= 0 None, = 1 Full, ==2 Summary)
//    public bool ConsoleDebugOutput = false; // If true send debug output to Console

namespace SALSALibrary
{
    [Serializable]
    public class SALSAFileProperties
    {
        public string ClusterName = ""; // Name of Cluster -- typically Run Label
        public int ClusterStartIndex = 0; // Cluster Indices start at 0 or 1
        public string Comment = ""; // General Comment

        public string FamilyName1 = ""; // Name of Family 1
        public string FamilyName2 = ""; // Name of Family 2

        public int FileGenerationType = 0;
                   // 0 Clustering; 1 Grouping; 2 MDS Original; 3 MDS Incremental; 4 MDS Fixed; 5 MDS Rotation

        public string GroupName = "Undefined Group"; // Name of Group
        public int LocalPointStartIndex = 0; // File Point Indices start at 0 or 1
        public int LocalVectorDimension = 3; // Vector dimension of Mapped Points
        public int NumberOriginalPoints = 0; // Number of Original Points
        public int NumberPointsinFile = 0; // Number of Points Written
        public int NumberRotationParameters = 0; // Number of Rotation Parameters
        public int OriginalPointStartIndex = 0; // Original Point Indices start at 0 or 1
        public string RotationParameters = ""; // List of Rotation parameters separated by commas

        public SALSAFileProperties ShallowCopy()
        {
            return (SALSAFileProperties) MemberwiseClone();
        }
    }

    // End SALSAFileProperties

    [Serializable]
    public class SALSADataPointProperties
    {
        public int FixedorVaried = 0; // = 0 Ignored = 1 Varied = 2 Fixed
        public int LocalPointNumber = -1; // Point Number in this File
        public int OriginalPointNumber = -1; // Point Number in Original analysis
        public int PointType = 0; // = 0 Normal, > 0 Center
        public int cluster = -1; // Index into Cluster
        public string clusterlabel = ""; // Label of this cluster member
        public bool errorsset = false; // If true x y z errors set in colon section
        public int family1 = -1; // Index into Family 1 -- all these indices UNDEFINED if negative
        public int family2 = -1; // Index into Family 2
        public string familylabel1 = ""; // Label of this family member
        public string familylabel2 = ""; // Label of this family member
        public int group = -1; // Index into Undefined Group
        public string grouplabel = ""; // Label of this group member
        public string pointlabel = ""; // Label of this point
        public string source = ""; // Source run producing  position values
        public bool valuesset = false; // If true x y z set in colon section
        public double x = 0.0; // x value
        public double xerr = 0.0; // x error value
        public double y = 0.0; // y value
        public double yerr = 0.0; // y error value
        public double z = 0.0; // z value
        public double zerr = 0.0; // z error value

        public SALSADataPointProperties ShallowCopy()
        {
            return (SALSADataPointProperties) MemberwiseClone();
        }
    }

    // End SALSADataPointProperties

    public class SALSA_Properties
    {
        public static bool ReadDataPointFile(string MDSClusterFileName, ref int FileType,
                                             SALSAFileProperties FileProperties,
                                             ref SALSADataPointProperties[] DataPoints, ref int NumberofPoints)
        {
            // FileType = 0 simple integers and positions:  PointNumber x y z LABEL or Group Number (LABEL if file name contains label)
            // FileType = 1 Just integers -- Point Number and Group Number: PointNumber  Group Number 
            // File Type = 2 Integer and Label: Point LABEL (LABEL if file name contains label)
            // File Type = 3 Pure colon style -- with  positions
            // File Type = 4 Pure Colon Style -- no    positions
            // File Type = 5 Pure Colon Style --    Labels
            // File Type = 6 Hybrid Style --     with positions
            // File Type = 7 Hybrid Style --     no positions
            // File Type = 8 Hybrid Style --     Labels


            bool positionsset = false;
            bool colonsinput = false;
            bool NonColonsInput = false;
            bool labelfile = false;
            string LowerFileName = MDSClusterFileName.ToLower();
            if (LowerFileName.Contains("label"))
                labelfile = true;

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
                    // Read contents of a file
                    String inputLineStr;
                    int newlabelnumber = -1;
                    int LocalStart = FileProperties.LocalPointStartIndex;
                    int OriginalStart = FileProperties.OriginalPointStartIndex;
                    while ((inputLineStr = sr.ReadLine()) != null)
                    {
                        if (inputLineStr.Length < 2)
                        {
                            continue; //replace empty line
                        }

                        inputLineStr = inputLineStr.Trim(new[] {' ', '\t'});
                        try
                        {
                            // Parse each record string
                            inputLineStr = inputLineStr.Replace("\t\t", "\t");

                            if (inputLineStr.Contains("FileProperties"))
                            {
                                // Not a Data Point
                                string EndofInput = inputLineStr.Replace("FileProperties", "");
                                ReadFileProperties(FileProperties, EndofInput);
                                colonsinput = true;
                                LocalStart = FileProperties.LocalPointStartIndex;
                                OriginalStart = FileProperties.OriginalPointStartIndex;
                            }
                            else
                            {
                                // Must be a Data Point
                                ArrayInitializer(ref DataPoints, SALSAUtility.NumberOriginalPoints,
                                                 FileProperties.NumberPointsinFile);
                                DataPoints[NumberofPoints] = new SALSADataPointProperties();
                                bool incrementcount = false;
                                int PointDataStarts = inputLineStr.IndexOf("PointProperties");
                                if (PointDataStarts >= 0)
                                {
                                    // Some Colon Information
                                    string EndofInput = inputLineStr.Substring(PointDataStarts);
                                    EndofInput = EndofInput.Replace("PointProperties", "");
                                    ReadDataPointProperties(DataPoints[NumberofPoints], EndofInput);
                                    colonsinput = true;
                                    if (DataPoints[NumberofPoints].valuesset)
                                        positionsset = true;
                                    incrementcount = true;
                                    DataPoints[NumberofPoints].LocalPointNumber = NumberofPoints +
                                                                                  FileProperties.LocalPointStartIndex;
                                } //  End Processing Colon Point Information

                                if (PointDataStarts < 0)
                                {
                                    // traditional bare line
                                    PointDataStarts = inputLineStr.Length;
                                }
                                if (PointDataStarts > 0)
                                {
                                    // Process number information
                                    string StartofString = inputLineStr.Substring(0, PointDataStarts);
                                    StartofString = StartofString.Trim(new[] {' ', '\t'});

                                    if (StartofString.Length > 0)
                                    {
                                        // You come here only for traditional bare line of x,y,z coordinates.

                                        string[] split = StartofString.Split(new[] {' ', '\t'},
                                                                             StringSplitOptions.RemoveEmptyEntries);

                                        if ((split.Length != 5) && (split.Length != 2) && (split.Length != 4))
                                        {
                                            Exception e =
                                                SALSAUtility.SALSAError(" Bad Line " + split.Length.ToString() + " "
                                                                        + NumberofPoints.ToString() + " " + inputLineStr);
                                            throw (e);
                                        }
                                        newlabelnumber = Convert.ToInt32(split[0]);

                                        if ((NumberofPoints + LocalStart) != newlabelnumber)
                                        {
                                            Exception e =
                                                SALSAUtility.SALSAError("Unexpected Label Number " + newlabelnumber +
                                                                        " Expected "
                                                                        + NumberofPoints.ToString() + " + " +
                                                                        LocalStart.ToString());
                                            throw (e);
                                        }
                                        if (DataPoints[NumberofPoints].LocalPointNumber >= 0)
                                        {
                                            if ((DataPoints[NumberofPoints].LocalPointNumber - LocalStart) !=
                                                NumberofPoints)
                                            {
                                                Exception e =
                                                    SALSAUtility.SALSAError("Unexpected Local Number " +
                                                                            DataPoints[NumberofPoints].LocalPointNumber +
                                                                            " Expected "
                                                                            + NumberofPoints.ToString() + " + " +
                                                                            LocalStart.ToString());
                                                throw (e);
                                            }
                                        }
                                        DataPoints[NumberofPoints].LocalPointNumber = NumberofPoints;
                                        if (DataPoints[NumberofPoints].OriginalPointNumber >= 0)
                                        {
                                            if ((DataPoints[NumberofPoints].OriginalPointNumber - OriginalStart) < 0)
                                            {
                                                Exception e =
                                                    SALSAUtility.SALSAError("Unexpected Original Number " +
                                                                            DataPoints[NumberofPoints].
                                                                                OriginalPointNumber + " Local Point "
                                                                            + NumberofPoints.ToString() +
                                                                            " Original Increment " +
                                                                            OriginalStart.ToString());
                                                throw (e);
                                            }
                                            DataPoints[NumberofPoints].OriginalPointNumber -= OriginalStart;
                                        }
                                        else
                                        {
                                            DataPoints[NumberofPoints].OriginalPointNumber = newlabelnumber;
                                        }

                                        if (labelfile)
                                            DataPoints[NumberofPoints].pointlabel = split[split.Length - 1];
                                        else
                                            DataPoints[NumberofPoints].group = Convert.ToInt32(split[split.Length - 1]);

                                        if (split.Length >= 4)
                                        {
                                            DataPoints[NumberofPoints].valuesset = true;
                                            DataPoints[NumberofPoints].x = Convert.ToDouble(split[1]);
                                            DataPoints[NumberofPoints].y = Convert.ToDouble(split[2]);
                                            positionsset = true;
                                        }
                                        if (split.Length == 5)
                                        {
                                            DataPoints[NumberofPoints].z = Convert.ToDouble(split[3]);
                                        }
                                        incrementcount = true;
                                        NonColonsInput = true;
                                    } // End Processing non colon Point information
                                }
                                if (incrementcount)
                                    ++NumberofPoints;
                            }
                        }
                        catch (Exception e)
                        {
                            SALSAUtility.SALSAError("Failed to load data array " + inputLineStr + " "
                                                    + " " + NumberofPoints.ToString() + " " + newlabelnumber.ToString() +
                                                    " " + e);
                            throw (e);
                        }
                    }

                    FileType = 1;
                    if (positionsset)
                        FileType = 0;
                    if (labelfile)
                        FileType = 2;
                    if (colonsinput)
                    {
                        if (NonColonsInput)
                            FileType += 6;
                        else
                            FileType += 3;
                    }
                    if (FileProperties.NumberOriginalPoints == 0)
                    {
                        FileProperties.NumberOriginalPoints = NumberofPoints;
                    }

                    if (FileProperties.NumberPointsinFile == 0)
                    {
                        FileProperties.NumberPointsinFile = NumberofPoints;
                    }

                    if (FileProperties.NumberPointsinFile != NumberofPoints)
                    {
                        Exception e =
                            SALSAUtility.SALSAError("Unexpected Number of Points in File " + NumberofPoints.ToString() +
                                                    " Read but Expected "
                                                    + FileProperties.NumberPointsinFile.ToString());
                        throw (e);
                    }
                    sr.Close();
                }
            }
            catch (Exception e)
            {
                SALSAUtility.SALSAError("Failed to read data from " + MDSClusterFileName + " " + e);
                throw (e);
            }
            return true;
        }

        // End ReadDataPointFile

        public static void ArrayInitializer(ref SALSADataPointProperties[] DataArray, int sizemax, int sizereadin)
        {
            if (DataArray != null)
            {
                if (DataArray.Length < sizereadin)
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Data Array too small for file Length " + DataArray.Length.ToString() +
                                                " Needs " + sizereadin.ToString());
                    throw (e);
                }
                return;
            }
            int size = sizereadin;
            if (size == 0)
                size = sizemax;
            DataArray = new SALSADataPointProperties[size];
            return;
        }

        // End ArrayInitializer

        // Write label-cluster results into a file
        public static void WriteDataPointFile(string CoreOutputFileName, bool write2Das3D, string OutputStyles,
                                              SALSAFileProperties FileProperties, SALSADataPointProperties[] DataPoints,
                                              int NumberofDataPoints)
        {
            var DothisOutputStyle = new bool[5];
            for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
            {
                DothisOutputStyle[StyleIndex] = false;
                if (OutputStyles.Contains("all"))
                    DothisOutputStyle[StyleIndex] = true;
            }
            if (OutputStyles.Contains("colon"))
                DothisOutputStyle[0] = true;
            if (OutputStyles.Contains("family1"))
                DothisOutputStyle[1] = true;
            if (OutputStyles.Contains("family2"))
                DothisOutputStyle[2] = true;
            if (OutputStyles.Contains("cluster"))
                DothisOutputStyle[3] = true;
            if (OutputStyles.Contains("group"))
                DothisOutputStyle[4] = true;

            bool setgroup = false;
            bool OutputValues = false;

            SALSADataPointProperties FirstOne = DataPoints[0];
            if (FirstOne.family1 < 0)
                DothisOutputStyle[1] = false;
            if (FirstOne.family2 < 0)
                DothisOutputStyle[2] = false;
            if (FirstOne.cluster < 0)
                DothisOutputStyle[3] = false;
            if (FirstOne.group < 0)
                DothisOutputStyle[4] = false;
            if ((FirstOne.family1 < 0) && (FirstOne.family2 < 0) && (FirstOne.cluster < 0) && (FirstOne.group < 0))
            {
                DothisOutputStyle[4] = true;
                setgroup = true;
            }
            if (FirstOne.valuesset)
                OutputValues = true;
            int LocalVectorDimension = FileProperties.LocalVectorDimension;
            int LocalPointIncrement = FileProperties.LocalPointStartIndex;

            int numberoffiles = 0;
            for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
            {
                if (DothisOutputStyle[StyleIndex])
                    ++numberoffiles;
            }
            if (numberoffiles == 0)
            {
                SALSAUtility.SALSAError("No files output for core name " + CoreOutputFileName);
                return;
            }
            if (numberoffiles > 1)
            {
                if (OutputStyles.Contains("SameFileName"))
                {
                    Exception e =
                        SALSAUtility.SALSAError("Attempt to generate multiple outputs to same file " +
                                                CoreOutputFileName);
                    throw (e);
                }
            }
            for (int StyleIndex = 0; StyleIndex < 5; StyleIndex++)
            {
                if (!DothisOutputStyle[StyleIndex])
                    continue;
                string OutputFileName = "";
                if (!OutputStyles.Contains("SameFileName"))
                {
                    if (StyleIndex == 0)
                        OutputFileName = CoreOutputFileName.Replace(".txt", "Colon.txt");
                    if (StyleIndex == 1)
                        OutputFileName = CoreOutputFileName.Replace(".txt", "Family1.txt");
                    if (StyleIndex == 2)
                        OutputFileName = CoreOutputFileName.Replace(".txt", "Family2.txt");
                    if (StyleIndex == 3)
                        OutputFileName = CoreOutputFileName.Replace(".txt", "Cluster.txt");
                    if (StyleIndex == 4)
                        OutputFileName = CoreOutputFileName.Replace(".txt", "Group.txt");
                }
                else
                {
                    OutputFileName = CoreOutputFileName;
                }

                try
                {
                    StreamWriter sw = null;
                    if (!string.IsNullOrEmpty(OutputFileName))
                    {
                        sw = new StreamWriter(OutputFileName, false, Encoding.UTF8);
                    }
                    if (sw != null)
                    {
                        if (StyleIndex == 0)
                            WriteFileProperties(FileProperties, sw); // Write Header of a Colon File

                        for (int GlobalDataPoint = 0; GlobalDataPoint < NumberofDataPoints; GlobalDataPoint++)
                        {
                            SALSADataPointProperties ThesePointProperties = DataPoints[GlobalDataPoint];
                            if ((LocalVectorDimension == 2) && write2Das3D && OutputValues)
                            {
                                ThesePointProperties.z = 0.0;
                                ThesePointProperties.zerr = 0.0;
                            }
                            if (StyleIndex == 0)
                            {
                                string OutputLine = "";
                                AppendDataPointProperties(ThesePointProperties, ref OutputLine);
                                sw.WriteLine(OutputLine);
                                continue;
                            }
                            int IntegerIndex = 0;
                            if (StyleIndex == 1)
                                IntegerIndex = ThesePointProperties.family1;
                            if (StyleIndex == 2)
                                IntegerIndex = ThesePointProperties.family2;
                            if (StyleIndex == 3)
                                IntegerIndex = ThesePointProperties.cluster;
                            if (StyleIndex == 4)
                            {
                                IntegerIndex = ThesePointProperties.group;
                                if (setgroup)
                                    IntegerIndex = 1;
                            }
                            string Coordinates = "";
                            if (OutputValues)
                            {
                                Coordinates = ThesePointProperties.x.ToString("E4") + "\t" +
                                              ThesePointProperties.y.ToString("E4") + "\t";
                                if ((LocalVectorDimension == 3) || write2Das3D)
                                    Coordinates += ThesePointProperties.z.ToString("E4") + "\t";
                            }
                            sw.WriteLine(
                                String.Format((GlobalDataPoint + LocalPointIncrement).ToString() + "\t" + Coordinates +
                                              IntegerIndex.ToString()));
                        }

                        sw.Flush();
                        sw.Close();
                    }
                }
                catch
                {
                    Exception e = SALSAUtility.SALSAError(" Failed to Write Properties File " + CoreOutputFileName);
                    throw (e);
                }
            } // End Loop over File Types
            return;
        }

        // End WriteDataPointFile


        public static void ReadFileProperties(SALSAFileProperties FileProps, string InputLine)
        {
            InputLine = InputLine.Trim(new[] {' ', '\t'});
            string[] colonsplit = InputLine.Split(new[] {':'}, 2);
            if (colonsplit.Length < 2)
            {
                Console.WriteLine("No deliminator in Line " + InputLine);
                return;
            }
            colonsplit[0] = colonsplit[0].Trim();
            colonsplit[1] = colonsplit[1].Trim();
            if (colonsplit[0].Equals("LocalVectorDimension"))
                FileProps.LocalVectorDimension = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("ClusterStartIndex"))
                FileProps.ClusterStartIndex = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("OriginalPointStartIndex"))
                FileProps.OriginalPointStartIndex = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("LocalPointStartIndex"))
                FileProps.LocalPointStartIndex = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("FileGenerationType"))
                FileProps.FileGenerationType = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("FamilyName1"))
                FileProps.FamilyName1 = colonsplit[1];
            if (colonsplit[0].Equals("FamilyName2"))
                FileProps.FamilyName2 = colonsplit[1];
            if (colonsplit[0].Equals("GroupName"))
                FileProps.GroupName = colonsplit[1];
            if (colonsplit[0].Equals("ClusterName"))
                FileProps.ClusterName = colonsplit[1];
            if (colonsplit[0].Equals("NumberOriginalPoints"))
                FileProps.NumberOriginalPoints = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("NumberPointsinFile"))
                FileProps.NumberPointsinFile = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("RotationParameters"))
                FileProps.RotationParameters = colonsplit[1];
            if (colonsplit[0].Equals("NumberRotationParameters"))
                FileProps.NumberRotationParameters = Convert.ToInt32(colonsplit[1]);
            if (colonsplit[0].Equals("Comment"))
            {
                if (FileProps.Comment != "")
                    FileProps.Comment += "\n";
                FileProps.Comment += colonsplit[1];
            }
            return;
        }

        // End ReadFileProperties

        public static void WriteFileProperties(SALSAFileProperties FileProps, StreamWriter sw)
        {
            sw.WriteLine(String.Format("FileProperties\tLocalVectorDimension:{0}", FileProps.LocalVectorDimension));
            sw.WriteLine(String.Format("FileProperties\tClusterStartIndex:{0}", FileProps.ClusterStartIndex));
            sw.WriteLine(String.Format("FileProperties\tOriginalPointStartIndex:{0}", FileProps.OriginalPointStartIndex));
            sw.WriteLine(String.Format("FileProperties\tLocalPointStartIndex:{0}", FileProps.LocalPointStartIndex));
            sw.WriteLine(String.Format("FileProperties\tFileGenerationType:{0}", FileProps.FileGenerationType));
            sw.WriteLine(String.Format("FileProperties\tFamilyName1:{0}", FileProps.FamilyName1));
            sw.WriteLine(String.Format("FileProperties\tFamilyName2:{0}", FileProps.FamilyName2));
            sw.WriteLine(String.Format("FileProperties\tGroupName:{0}", FileProps.GroupName));
            sw.WriteLine(String.Format("FileProperties\tClusterName:{0}", FileProps.ClusterName));
            sw.WriteLine(String.Format("FileProperties\tNumberOriginalPoints:{0}", FileProps.NumberOriginalPoints));
            sw.WriteLine(String.Format("FileProperties\tNumberPointsinFile:{0}", FileProps.NumberPointsinFile));
            sw.WriteLine(String.Format("FileProperties\tRotationParameters:{0}", FileProps.RotationParameters));
            sw.WriteLine(String.Format("FileProperties\tNumberRotationParameters:{0}",
                                       FileProps.NumberRotationParameters));

            string FileComments = FileProps.Comment;
            if (FileComments != "")
            {
                string[] split = FileComments.Split(new[] {'\n'});
                for (int runthroughcomments = 0; runthroughcomments < split.Length; runthroughcomments++)
                {
                    if (split[runthroughcomments] != "")
                        sw.WriteLine(String.Format("FileProperties\tComment:{0}", split[runthroughcomments]));
                }
            }
            return;
        }

        // End WriteFileProperties

        public static void ReadDataPointProperties(SALSADataPointProperties DataPointProps, string InputLine)
        {
            InputLine = InputLine.Trim(new[] {' ', '\t'});
            string[] split = InputLine.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
            DataPointProps.valuesset = false;
            DataPointProps.errorsset = false;

            for (int itembackwards = 0; itembackwards < split.Length; itembackwards++)
            {
                int item = split.Length - itembackwards - 1;
                if (split[item].Equals("PointProperties"))
                    continue;
                string[] colonsplit = split[item].Split(new[] {':'}, 2);
                if (colonsplit[0].Equals("PointProperties"))
                    continue;
                if (colonsplit.Length < 2)
                {
                    if (item == 0)
                        Console.WriteLine("No deliminator in item " + split[item] + " in Line " + InputLine);
                    else
                    {
                        split[item - 1] += " " + split[item];
                        continue;
                    }
                }
                colonsplit[0] = colonsplit[0].Trim();
                colonsplit[1] = colonsplit[1].Trim();
                if (colonsplit[0].Equals("x"))
                {
                    DataPointProps.x = Convert.ToDouble(colonsplit[1]);
                    DataPointProps.valuesset = true;
                }
                if (colonsplit[0].Equals("y"))
                    DataPointProps.y = Convert.ToDouble(colonsplit[1]);
                if (colonsplit[0].Equals("z"))
                    DataPointProps.z = Convert.ToDouble(colonsplit[1]);
                if (colonsplit[0].Equals("xerr"))
                {
                    DataPointProps.xerr = Convert.ToDouble(colonsplit[1]);
                    DataPointProps.errorsset = true;
                }
                if (colonsplit[0].Equals("source"))
                    DataPointProps.source = colonsplit[1];
                if (colonsplit[0].Equals("yerr"))
                    DataPointProps.yerr = Convert.ToDouble(colonsplit[1]);
                if (colonsplit[0].Equals("zerr"))
                    DataPointProps.zerr = Convert.ToDouble(colonsplit[1]);
                if (colonsplit[0].Equals("family1"))
                    DataPointProps.family1 = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("familylabel1"))
                    DataPointProps.familylabel1 = colonsplit[1];
                if (colonsplit[0].Equals("family2"))
                    DataPointProps.family2 = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("familylabel2"))
                    DataPointProps.familylabel2 = colonsplit[1];
                if (colonsplit[0].Equals("cluster"))
                    DataPointProps.cluster = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("clusterlabel"))
                    DataPointProps.clusterlabel = colonsplit[1];
                if (colonsplit[0].Equals("group"))
                    DataPointProps.group = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("grouplabel"))
                    DataPointProps.grouplabel = colonsplit[1];
                if (colonsplit[0].Equals("pointlabel"))
                    DataPointProps.pointlabel = colonsplit[1];
                if (colonsplit[0].Equals("FixedorVaried"))
                    DataPointProps.FixedorVaried = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("PointType"))
                    DataPointProps.PointType = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("LocalPointNumber"))
                    DataPointProps.LocalPointNumber = Convert.ToInt32(colonsplit[1]);
                if (colonsplit[0].Equals("OriginalPointNumber"))
                    DataPointProps.OriginalPointNumber = Convert.ToInt32(colonsplit[1]);
            }
            return;
        }

        // End ReadDataPointProperties

        public static void AppendDataPointProperties(SALSADataPointProperties DataPointProps, ref string InputLine)
        {
            //Calling Routine must supply any needed delimiter before this is appended
            InputLine += "PointProperties";
            if (DataPointProps.valuesset)
            {
                // x y z are set

                InputLine += "\tvaluesset:true\tsource:" + DataPointProps.source + "\tx:" +
                             DataPointProps.x.ToString("E4");
                InputLine += "\ty:" + DataPointProps.y.ToString("E4");
                InputLine += "\tz:" + DataPointProps.z.ToString("E4");
            }
            if (DataPointProps.errorsset)
            {
                // x y z errors are set
                InputLine += "\terrorsset:true\txerr:" + DataPointProps.xerr.ToString("E4");
                InputLine += "\tyerr:" + DataPointProps.yerr.ToString("E4");
                InputLine += "\tzerr:" + DataPointProps.zerr.ToString("E4");
            }
            if (DataPointProps.family1 >= 0)
            {
                InputLine += "\tfamily1:" + DataPointProps.family1.ToString();
                if (DataPointProps.familylabel1.Length > 0)
                    InputLine += "\tfamilylabel1:" + DataPointProps.familylabel1;
            }
            if (DataPointProps.family2 >= 0)
            {
                InputLine += "\tfamily2:" + DataPointProps.family2.ToString();
                if (DataPointProps.familylabel2.Length > 0)
                    InputLine += "\tfamilylabel2:" + DataPointProps.familylabel2;
            }
            if (DataPointProps.cluster >= 0)
            {
                InputLine += "\tcluster:" + DataPointProps.cluster.ToString();
                if (DataPointProps.clusterlabel.Length > 0)
                    InputLine += "\tclusterlabel:" + DataPointProps.clusterlabel;
            }
            if (DataPointProps.group >= 0)
            {
                InputLine += "\tgroup:" + DataPointProps.group.ToString();
                if (DataPointProps.grouplabel.Length > 0)
                    InputLine += "\tgrouplabel:" + DataPointProps.grouplabel;
            }
            if (DataPointProps.pointlabel.Length > 0)
                InputLine += "\tpointlabel:" + DataPointProps.pointlabel;
            InputLine += "\tFixedorVaried:" + DataPointProps.FixedorVaried.ToString();
            InputLine += "\tPointType:" + DataPointProps.PointType.ToString();
            if (DataPointProps.LocalPointNumber >= 0)
                InputLine += "\tLocalPointNumber:" + DataPointProps.LocalPointNumber.ToString();
            if (DataPointProps.OriginalPointNumber >= 0)
                InputLine += "\tOriginalPointNumber:" + DataPointProps.OriginalPointNumber.ToString();
        }

        // End AppendDataPointProperties
    }

    // End SALSA_Properties
}