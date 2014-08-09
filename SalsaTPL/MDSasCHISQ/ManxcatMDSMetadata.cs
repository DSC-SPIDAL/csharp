using System;
using System.Collections;
using System.IO;
using System.Text;

namespace Manxcat
{
    public class MetaDataIO
    {
        #region Deprecated Code

        //public class RunMetadata
        //{
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

        //    public static void Membercopy(ref RunMetadata One, ref RunMetadata Two)
        //    {
        //        Two = (Manxcat.RunMetadata)One.MemberwiseClone();
        //        return;
        //    }
        //}

        // Write Control data into a file
        //public static void WriteControl_Cluster(string ControlFileName, ref RunMetadata DataforthisRun)
        //{
        //    try
        //    {
        //        StreamWriter FileToWrite = null;
        //        if (!string.IsNullOrEmpty(ControlFileName))
        //        {
        //            FileToWrite = new StreamWriter(ControlFileName, false, Encoding.UTF8); // Overwrite
        //        }

        //        if (FileToWrite != null)
        //        {
        //            FileToWrite.WriteLine(String.Format("DataFileName:{0}", DataforthisRun.DataFileName));
        //            FileToWrite.WriteLine(String.Format("DataLabelsFileName:{0}", DataforthisRun.DataLabelsFileName));
        //            FileToWrite.WriteLine(String.Format("ReducedVectorOutputFileName:{0}", DataforthisRun.ReducedVectorOutputFileName));
        //            FileToWrite.WriteLine(String.Format("ClusterDirectory:{0}", DataforthisRun.ClusterDirectory));
        //            FileToWrite.WriteLine(String.Format("DataPoints:{0}", DataforthisRun.DataPoints));
        //            FileToWrite.WriteLine(String.Format("ProcessingOption:{0}", DataforthisRun.ProcessingOption));
        //            FileToWrite.WriteLine(String.Format("DistanceProcessingOption:{0}", DataforthisRun.DistanceProcessingOption));
        //            FileToWrite.WriteLine(String.Format("InitializationOption:{0}", DataforthisRun.InitializationOption));
        //            FileToWrite.WriteLine(String.Format("InitializationFileName:{0}", DataforthisRun.InitializationFileName));
        //            FileToWrite.WriteLine(String.Format("WeightingOption:{0}", DataforthisRun.WeightingOption));
        //            FileToWrite.WriteLine(String.Format("WeightingFileName:{0}", DataforthisRun.WeightingFileName));
        //            FileToWrite.WriteLine(String.Format("write2Das3D:{0}", DataforthisRun.write2Das3D));
        //            FileToWrite.WriteLine(String.Format("LocalVectorDimension:{0}", DataforthisRun.LocalVectorDimension));
        //            FileToWrite.WriteLine(String.Format("selectedvariedpoints:{0}", DataforthisRun.selectedvariedpoints));
        //            FileToWrite.WriteLine(String.Format("VariedPointCriterion:{0}", DataforthisRun.VariedPointCriterion));
        //            FileToWrite.WriteLine(String.Format("selectedvariedpointfile:{0}", DataforthisRun.selectedvariedpointfile));
        //            FileToWrite.WriteLine(String.Format("selectedfixedpoints:{0}", DataforthisRun.selectedfixedpoints));
        //            FileToWrite.WriteLine(String.Format("FixedPointCriterion:{0}", DataforthisRun.FixedPointCriterion));
        //            FileToWrite.WriteLine(String.Format("selectedfixedpointfile:{0}", DataforthisRun.selectedfixedpointfile));
        //            FileToWrite.WriteLine(String.Format("ConversionOption:{0}", DataforthisRun.ConversionOption));
        //            FileToWrite.WriteLine(String.Format("ConversionInformation:{0}", DataforthisRun.ConversionInformation));

        //            FileToWrite.WriteLine(String.Format("RotationOption:{0}", DataforthisRun.RotationOption));
        //            FileToWrite.WriteLine(String.Format("RotationLabelsFileName:{0}", DataforthisRun.RotationLabelsFileName));
        //            FileToWrite.WriteLine(String.Format("InitializationLoops:{0}", DataforthisRun.InitializationLoops));

        //            FileToWrite.WriteLine(String.Format("Chisqnorm:{0}", DataforthisRun.Chisqnorm));
        //            FileToWrite.WriteLine(String.Format("DistanceFormula:{0}", DataforthisRun.DistanceFormula));
        //            FileToWrite.WriteLine(String.Format("FullSecondDerivativeOption:{0}", DataforthisRun.FullSecondDerivativeOption));
        //            FileToWrite.WriteLine(String.Format("MinimumDistance:{0}", DataforthisRun.MinimumDistance));
        //            FileToWrite.WriteLine(String.Format("FunctionErrorCalcMultiplier:{0}", DataforthisRun.FunctionErrorCalcMultiplier));
        //            FileToWrite.WriteLine(String.Format("ChisqPrintConstant:{0}", DataforthisRun.ChisqPrintConstant));
        //            FileToWrite.WriteLine(String.Format("maxit:{0}", DataforthisRun.maxit));
        //            FileToWrite.WriteLine(String.Format("nbadgo:{0}", DataforthisRun.nbadgo));
        //            FileToWrite.WriteLine(String.Format("ChisqChangePerPoint:{0}", DataforthisRun.ChisqChangePerPoint));
        //            FileToWrite.WriteLine(String.Format("FletcherRho:{0}", DataforthisRun.FletcherRho));
        //            FileToWrite.WriteLine(String.Format("FletcherSigma:{0}", DataforthisRun.FletcherSigma));
        //            FileToWrite.WriteLine(String.Format("Omega:{0}", DataforthisRun.Omega));
        //            FileToWrite.WriteLine(String.Format("OmegaOption:{0}", DataforthisRun.OmegaOption));
        //            FileToWrite.WriteLine(String.Format("QHighInitialFactor:{0}", DataforthisRun.QHighInitialFactor));
        //            FileToWrite.WriteLine(String.Format("QgoodReductionFactor:{0}", DataforthisRun.QgoodReductionFactor));
        //            FileToWrite.WriteLine(String.Format("QLimitscalculationInterval:{0}", DataforthisRun.QLimitscalculationInterval));
        //            FileToWrite.WriteLine(String.Format("extraprecision:{0}", DataforthisRun.extraprecision));
        //            FileToWrite.WriteLine(String.Format("addonforQcomputation:{0}", DataforthisRun.addonforQcomputation));
        //            FileToWrite.WriteLine(String.Format("InitialSteepestDescents:{0}", DataforthisRun.InitialSteepestDescents));
        //            FileToWrite.WriteLine(String.Format("TimeCutmillisec:{0}", DataforthisRun.TimeCutmillisec));
        //            FileToWrite.WriteLine(String.Format("CGResidualLimit:{0}", DataforthisRun.CGResidualLimit));
        //            FileToWrite.WriteLine(String.Format("PowerIterationLimit:{0}", DataforthisRun.PowerIterationLimit));
        //            FileToWrite.WriteLine(String.Format("eigenvaluechange:{0}", DataforthisRun.eigenvaluechange));
        //            FileToWrite.WriteLine(String.Format("eigenvectorchange:{0}", DataforthisRun.eigenvectorchange));
        //            FileToWrite.WriteLine(String.Format("derivtest:{0}", DataforthisRun.derivtest));

        //            FileToWrite.WriteLine(String.Format("RunNumber:{0}", DataforthisRun.RunNumber));
        //            FileToWrite.WriteLine(String.Format("BaseResultDirectoryName:{0}", DataforthisRun.BaseResultDirectoryName));
        //            FileToWrite.WriteLine(String.Format("ResultDirectoryExtension:{0}", DataforthisRun.ResultDirectoryExtension));
        //            FileToWrite.WriteLine(String.Format("TimingOutputFileName:{0}", DataforthisRun.TimingOutputFileName));
        //            FileToWrite.WriteLine(String.Format("RunSetLabel:{0}", DataforthisRun.RunSetLabel));
        //            FileToWrite.WriteLine(String.Format("ControlDirectoryName:{0}", DataforthisRun.ControlDirectoryName));


        //            FileToWrite.WriteLine(String.Format("pattern:{0}", DataforthisRun.pattern));
        //            FileToWrite.WriteLine(String.Format("ThreadCount:{0}", DataforthisRun.ThreadCount));
        //            FileToWrite.WriteLine(String.Format("NodeCount:{0}", DataforthisRun.NodeCount));
        //            FileToWrite.WriteLine(String.Format("MPIperNodeCount:{0}", DataforthisRun.MPIperNodeCount));
        //            FileToWrite.WriteLine(String.Format("MPIIOStrategy:{0}", DataforthisRun.MPIIOStrategy));

        //            FileToWrite.WriteLine(String.Format("HistogramBinCount:{0}", DataforthisRun.HistogramBinCount));
        //            FileToWrite.WriteLine(String.Format("Extradata1:{0}", DataforthisRun.Extradata1));
        //            FileToWrite.WriteLine(String.Format("Extradata2:{0}", DataforthisRun.Extradata2));
        //            FileToWrite.WriteLine(String.Format("Extradata3:{0}", DataforthisRun.Extradata3));
        //            FileToWrite.WriteLine(String.Format("Extradata4:{0}", DataforthisRun.Extradata4));
        //            FileToWrite.WriteLine(String.Format("ExtraOption1:{0}", DataforthisRun.ExtraOption1));
        //            FileToWrite.WriteLine(String.Format("DebugPrintOption:{0}", DataforthisRun.DebugPrintOption));
        //            FileToWrite.WriteLine(String.Format("ConsoleDebugOutput:{0}", DataforthisRun.ConsoleDebugOutput));

        //            string usercomments = ManxcatCentral.Configuration.Comment; // todo: saliya - handle comments afterwards
        //            if (usercomments != "")
        //            {
        //                string[] split = usercomments.Split(new char[] { '\n' });
        //                for (int runthroughcomments = 0; runthroughcomments < split.Length; runthroughcomments++)
        //                {
        //                    if (split[runthroughcomments] != "")
        //                        FileToWrite.WriteLine(String.Format("Comment:{0}", split[runthroughcomments]));
        //                }
        //            }
        //        }

        //        FileToWrite.Flush();
        //        FileToWrite.Close();
        //    }
        //    catch (Exception e)
        //    {
        //        Console.WriteLine("Failed writing data" + e);
        //    }
        //}   // End WriteControl_Cluster

        // Write Results summary data into a file

        //  Read Control Data
        //public static void ReadControl_Cluster(string ControlFileName, ref RunMetadata DataforthisRun)
        //{
        //    DataforthisRun.DataFileName = "";
        //    DataforthisRun.DataLabelsFileName = "";
        //    DataforthisRun.ReducedVectorOutputFileName = "";
        //    DataforthisRun.ClusterDirectory = "";
        //    DataforthisRun.DataPoints = 0;
        //    DataforthisRun.ProcessingOption = 0;
        //    DataforthisRun.DistanceProcessingOption = 0;
        //    DataforthisRun.InitializationOption = 0;
        //    DataforthisRun.InitializationFileName = "";
        //    DataforthisRun.write2Das3D = true; // for 2D MDS write output for 3D viewer
        //    DataforthisRun.LocalVectorDimension = 3;    // Final MDS dimension (2 or 3 expected)
        //    DataforthisRun.selectedvariedpoints = "";  // List of varied points specified by point number or cluster number 
        //    DataforthisRun.VariedPointCriterion = "all"; // If 0 all the non-fixed points varied; if > 0 choose a subset of points to be varied; 1 =1 select by Original Point Number; =2 Select by group (family, cluster)
        //    DataforthisRun.selectedvariedpointfile = "";  // File containing mapping of groups to original values 
        //    DataforthisRun.selectedfixedpoints = "";  // List of fixed points specified by point number or group number 
        //    DataforthisRun.FixedPointCriterion = "none"; // If 0 there are no fixed points; if > 0 choose a subset of points to be fixed; 1 =1 select by Original Point Number; =2 Select by group (family, cluster)
        //    DataforthisRun.selectedfixedpointfile = "";  // File containing distances of fixed points

        //    DataforthisRun.RotationOption = 0;  // Control Rotation 
        //    DataforthisRun.RotationLabelsFileName = "";   // File with application specific labels for data to be rotated

        //    DataforthisRun.Chisqnorm = 0;
        //    DataforthisRun.DistanceFormula = 0;
        //    DataforthisRun.FullSecondDerivativeOption = 0;
        //    DataforthisRun.MinimumDistance = -0.001;   // Minimum Distance for Chisq Normalization
        //    DataforthisRun.FunctionErrorCalcMultiplier = 10.0; // Multiplicative Constant for Function (square this for Chisq) Error Calculation
        //    DataforthisRun.ChisqPrintConstant = 1.0; // Multiplicative Constant for Chisq Printing
        //    DataforthisRun.maxit = 20;  // Maximum Number of Manxcat Iterations
        //    DataforthisRun.nbadgo = 6;  // Limit on Iterations with no progress
        //    DataforthisRun.ChisqChangePerPoint = 0.001;  // Stop when this is Chisq change per point
        //    DataforthisRun.FletcherRho = 0.25;   // Fletcher's parameter rho
        //    DataforthisRun.FletcherSigma = 0.75; // Fletcher's Parameter sigma
        //    DataforthisRun.Omega = 1.25;
        //    DataforthisRun.OmegaOption = 0;
        //    DataforthisRun.QHighInitialFactor = 0.01; // Initial call uses this fraction of Qhigh for Q (added 2009)
        //    DataforthisRun.QgoodReductionFactor = 0.5;    // Set to zero to make Q = 0. for good solutions
        //    DataforthisRun.QLimitscalculationInterval = 1; //  recalculate Qhigh Qlow Trace and Norm after this number of iterations
        //    DataforthisRun.addonforQcomputation = 2.0;   // Add this to make Power Method Safer
        //    DataforthisRun.extraprecision = 0.05;
        //    DataforthisRun.InitialSteepestDescents = 0;  // Start minimization with this number of Steepest Descent Iterations
        //    DataforthisRun.TimeCutmillisec = -1.0;  // if positive, Stop when time used in milliseconds greater than timmax
        //    DataforthisRun.CGResidualLimit = 0.00001;    // Stop when Norm of Residual is less than this
        //    DataforthisRun.PowerIterationLimit = 200;    // Limit on Power Iterations
        //    DataforthisRun.eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
        //    DataforthisRun.eigenvectorchange = 0.001; // convergence test in power method for eigenvector
        //    DataforthisRun.derivtest = false;


        //    DataforthisRun.TimingOutputFileName = "";
        //    DataforthisRun.BaseResultDirectoryName = "";
        //    DataforthisRun.ResultDirectoryExtension = "";
        //    DataforthisRun.ControlDirectoryName = "";
        //    DataforthisRun.RunSetLabel = "";
        //    DataforthisRun.RunNumber = 0;

        //    DataforthisRun.pattern = "";
        //    DataforthisRun.ThreadCount = 1;    // Number of Threads
        //    DataforthisRun.NodeCount = 1;    // Number of Nodes
        //    DataforthisRun.MPIperNodeCount = 1; // Number of MPI processes per Node
        //    DataforthisRun.MPIIOStrategy = 0;   // Controls strategy of file handling with MPI =0 is ONE FILE

        //    DataforthisRun.HistogramBinCount = 100;
        //    DataforthisRun.Comment = "";
        //    DataforthisRun.Extradata1 = "";
        //    DataforthisRun.Extradata2 = "";
        //    DataforthisRun.Extradata3 = "";
        //    DataforthisRun.Extradata4 = "";
        //    DataforthisRun.ExtraOption1 = 0;
        //    DataforthisRun.DebugPrintOption = 1;
        //    DataforthisRun.ConsoleDebugOutput = false;

        //    try
        //    {
        //        // Check if file exists
        //        if (!File.Exists(ControlFileName))
        //        {
        //            Console.WriteLine("File {0} does not exist", ControlFileName);
        //            return;
        //        }
        //        // Create a new stream to read from a file
        //        using (StreamReader sr = File.OpenText(ControlFileName))
        //        {
        //            // Read contents of a file, line by line, into a string
        //            String inputLineStr;
        //            while ((inputLineStr = sr.ReadLine()) != null)
        //            {
        //                if (inputLineStr.Length < 2) continue; // Ignore empty line

        //                inputLineStr.Trim();    // Remove all beginning and Trailing white spaces
        //                try
        //                {
        //                    // Parse each record string
        //                    string[] split = inputLineStr.Split(new char[] { ':' }, 2);
        //                    if (split.Length < 2)
        //                    {
        //                        Console.WriteLine("No deliminator in " + inputLineStr);
        //                    }
        //                    split[0] = split[0].Trim();
        //                    split[1] = split[1].Trim();
        //                    if (split[0].Equals("DataFileName"))
        //                        DataforthisRun.DataFileName = split[1];
        //                    if (split[0].Equals("DataLabelsFileName"))
        //                        DataforthisRun.DataLabelsFileName = split[1];
        //                    if (split[0].Equals("ReducedVectorOutputFileName"))
        //                        DataforthisRun.ReducedVectorOutputFileName = split[1];
        //                    if (split[0].Equals("ClusterDirectory"))
        //                        DataforthisRun.ClusterDirectory = split[1];
        //                    if (split[0].Equals("DataPoints"))
        //                        DataforthisRun.DataPoints = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("ProcessingOption"))
        //                        DataforthisRun.ProcessingOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("DistanceProcessingOption"))
        //                        DataforthisRun.DistanceProcessingOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("InitializationOption"))
        //                        DataforthisRun.InitializationOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("InitializationFileName"))
        //                        DataforthisRun.InitializationFileName = split[1];
        //                    if (split[0].Equals("write2Das3D"))
        //                        DataforthisRun.write2Das3D = Convert.ToBoolean(split[1]);
        //                    if (split[0].Equals("LocalVectorDimension"))
        //                        DataforthisRun.LocalVectorDimension = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("selectedvariedpoints"))
        //                        DataforthisRun.selectedvariedpoints = split[1];
        //                    if (split[0].Equals("VariedPointCriterion"))
        //                        DataforthisRun.VariedPointCriterion = split[1].ToLower();
        //                    if (split[0].Equals("selectedvariedpointfile"))
        //                        DataforthisRun.selectedvariedpointfile = split[1];
        //                    if (split[0].Equals("selectedfixedpoints"))
        //                        DataforthisRun.selectedfixedpoints = split[1];
        //                    if (split[0].Equals("FixedPointCriterion"))
        //                        DataforthisRun.FixedPointCriterion = split[1].ToLower();
        //                    if (split[0].Equals("selectedfixedpointfile"))
        //                        DataforthisRun.selectedfixedpointfile = split[1];

        //                    if (split[0].Equals("RotationOption"))
        //                        DataforthisRun.RotationOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("RotationLabelsFileName"))
        //                        DataforthisRun.RotationLabelsFileName = split[1];

        //                    if (split[0].Equals("Chisqnorm"))
        //                        DataforthisRun.Chisqnorm = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("DistanceFormula"))
        //                        DataforthisRun.DistanceFormula = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("FullSecondDerivativeOption"))
        //                        DataforthisRun.FullSecondDerivativeOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("MinimumDistance"))
        //                        DataforthisRun.MinimumDistance = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("FunctionErrorCalcMultiplier"))
        //                        DataforthisRun.FunctionErrorCalcMultiplier = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("ChisqPrintConstant"))
        //                        DataforthisRun.ChisqPrintConstant = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("maxit"))
        //                        DataforthisRun.maxit = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("nbadgo"))
        //                        DataforthisRun.nbadgo = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("ChisqChangePerPoint"))
        //                        DataforthisRun.ChisqChangePerPoint = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("FletcherRho"))
        //                        DataforthisRun.FletcherRho = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("FletcherSigma"))
        //                        DataforthisRun.FletcherSigma = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("Omega"))
        //                        DataforthisRun.Omega = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("OmegaOption"))
        //                        DataforthisRun.OmegaOption = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("QHighInitialFactor"))
        //                        DataforthisRun.QHighInitialFactor = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("QgoodReductionFactor"))
        //                        DataforthisRun.QgoodReductionFactor = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("QLimitscalculationInterval"))
        //                        DataforthisRun.QLimitscalculationInterval = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("addonforQcomputation"))
        //                        DataforthisRun.addonforQcomputation = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("extraprecision"))
        //                        DataforthisRun.extraprecision = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("InitialSteepestDescents"))
        //                        DataforthisRun.InitialSteepestDescents = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("TimeCutmillisec"))
        //                        DataforthisRun.TimeCutmillisec = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("CGResidualLimit"))
        //                        DataforthisRun.CGResidualLimit = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("PowerIterationLimit"))
        //                        DataforthisRun.PowerIterationLimit = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("eigenvaluechange"))
        //                        DataforthisRun.eigenvaluechange = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("eigenvectorchange"))
        //                        DataforthisRun.eigenvectorchange = Convert.ToDouble(split[1]);
        //                    if (split[0].Equals("derivtest"))
        //                        DataforthisRun.derivtest = Convert.ToBoolean(split[1]);


        //                    if (split[0].Equals("TimingOutputFileName"))
        //                        DataforthisRun.TimingOutputFileName = split[1];
        //                    if (split[0].Equals("BaseResultDirectoryName"))
        //                        DataforthisRun.BaseResultDirectoryName = split[1];
        //                    if (split[0].Equals("ResultDirectoryExtension"))
        //                        DataforthisRun.ResultDirectoryExtension = split[1];
        //                    if (split[0].Equals("ControlDirectoryName"))
        //                        DataforthisRun.ControlDirectoryName = split[1];
        //                    if (split[0].Equals("RunSetLabel"))
        //                        DataforthisRun.RunSetLabel = split[1];
        //                    if (split[0].Equals("RunNumber"))
        //                        DataforthisRun.RunNumber = Convert.ToInt32(split[1]);

        //                    if (split[0].Equals("pattern"))
        //                        DataforthisRun.pattern = split[1];
        //                    if (split[0].Equals("ThreadCount"))
        //                        DataforthisRun.ThreadCount = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("NodeCount"))
        //                        DataforthisRun.NodeCount = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("MPIperNodeCount"))
        //                        DataforthisRun.MPIperNodeCount = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("MPIIOStrategy"))
        //                        DataforthisRun.MPIIOStrategy = Convert.ToInt32(split[1]);

        //                    if (split[0].Equals("HistogramBinCount"))
        //                        DataforthisRun.HistogramBinCount = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("Comment"))
        //                    {
        //                        if (DataforthisRun.Comment != "")
        //                            DataforthisRun.Comment += "\n";
        //                        DataforthisRun.Comment += split[1];
        //                    }
        //                    if (split[0].Equals("Extradata1"))
        //                        DataforthisRun.Extradata1 = split[1];
        //                    if (split[0].Equals("Extradata2"))
        //                        DataforthisRun.Extradata2 = split[1];
        //                    if (split[0].Equals("Extradata3"))
        //                        DataforthisRun.Extradata3 = split[1];
        //                    if (split[0].Equals("Extradata4"))
        //                        DataforthisRun.Extradata4 = split[1];
        //                    if (split[0].Equals("ExtraOption1"))
        //                        DataforthisRun.ExtraOption1 = Convert.ToInt32(split[1]);
        //                    if (split[0].Equals("DebugPrintOption"))
        //                        DataforthisRun.DebugPrintOption = Convert.ToInt32(split[1]); ;
        //                    if (split[0].Equals("ConsoleDebugOutput"))
        //                        DataforthisRun.ConsoleDebugOutput = Convert.ToBoolean(split[1]);

        //                }
        //                catch
        //                {
        //                    Console.WriteLine("Failed to load Control File");
        //                }
        //            }   // End while reading lines
        //            sr.Close();
        //        }
        //    }
        //    catch (Exception e)
        //    {
        //        Console.WriteLine("Failed to read Control File " + e);
        //    }
        //}   // End ReadControl_Cluster

        #endregion

        public static void WriteResults_Cluster(string ControlFileName, ArrayList linestooutput)
        {
            try
            {
                StreamWriter FileToWrite = null;

                if (!string.IsNullOrEmpty(ControlFileName))
                {
                    FileToWrite = new StreamWriter(ControlFileName, true, Encoding.UTF8); // Overwrite
                }


                if (FileToWrite != null)
                {
                    FileToWrite.WriteLine();

                    int LineLength = linestooutput.Count;

                    for (int countlines = 0; countlines < LineLength; countlines++)
                    {
                        FileToWrite.WriteLine(linestooutput[countlines]);
                    }
                }
                FileToWrite.Flush();
                FileToWrite.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing data {0} " + e, ControlFileName);
                throw (e);
            }
        }

        // End WriteResults_Cluster
    }

    // End class MetaDataIO
}

// End Namespace MDS