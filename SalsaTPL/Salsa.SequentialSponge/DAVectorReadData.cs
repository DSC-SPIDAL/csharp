using System;
using System.IO;
using MPI;
using Salsa.Core.Blas;
using Salsa.Core;

using Salsa.DAVectorSponge;

namespace SALSALibrary
{
    public class DAVectorReadData
    {
        public static void ReadLabelsFromFile(string fname)
        {
            char[] _sep = new[] { ' ', '\t' };

            int MinSplitSize = 8;
            int SplitPosition = 3;
            int MedeaPosition = 7;
            int GoldenIDPosition = 5;
            int MclustPosition = 6;
            int GoldenLabelPosition = 4;

            bool success = false;
            int count = 0;

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                    sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    string line;
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            if (splits.Length < MinSplitSize)
                            {
                                Exception e = DAVectorUtility.SALSAError("Count " + count.ToString() + "Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                throw (e);
                            }
                            int Charge = -101;
                            if (!Int32.TryParse(splits[SplitPosition], out Charge))
                                continue;
                            if (Charge != Program.SelectedInputLabel)
                                continue;
                            int OKinteger = -1;
                            if (!Int32.TryParse(splits[GoldenIDPosition], out OKinteger))
                                OKinteger = -1;
                            Program.GoldenPeaks.PointstoClusterIDs[count] = OKinteger;
                            GoldenExamination.GoldenID[count] = OKinteger;
                            GoldenExamination.GoldenLabel[count] = splits[GoldenLabelPosition];
                            GoldenExamination.PeakPosition[count][0] = double.Parse(splits[1]);
                            GoldenExamination.PeakPosition[count][1] = double.Parse(splits[2]);

                            OKinteger = 0;
                            if (!Int32.TryParse(splits[MclustPosition], out OKinteger))
                            {
                                double tryagain = 0.0;
                                if (!Double.TryParse(splits[MclustPosition], out tryagain))
                                {
                                    Exception e = DAVectorUtility.SALSAError("Count " + count.ToString() + "Illegal Mclust value on Point file " + splits.Length.ToString()
                                         + " " + line);
                                    throw (e);
                                }
                                OKinteger = (int) Math.Floor(tryagain + 0.001);
                            }
                            Program.MclustClusters.PointstoClusterIDs[count] = OKinteger;

                            OKinteger = 0;
                            if (!Int32.TryParse(splits[MedeaPosition], out OKinteger))
                            {
                                double tryagain = 0.0;
                                if (!Double.TryParse(splits[MedeaPosition], out tryagain))
                                {
                                    Exception e = DAVectorUtility.SALSAError("Count " + count.ToString() + "Illegal Medea value on Point file " + splits.Length.ToString()
                                         + " " + line);
                                    throw (e);
                                }
                                OKinteger = (int)Math.Floor(tryagain + 0.001);
                            }
                            Program.MedeaClusters.PointstoClusterIDs[count] = OKinteger;
                            count++;
                        }
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Points data" + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File Analyze error " + fname);
                throw (e);
            }
            Program.GoldenPeaks.setup();
            Program.MclustClusters.setup();
            Program.MedeaClusters.setup();

        }   // End ReadLabelsFromFile

        // read point data from file to memory
        // OutputFileType = 0 Harvard Format Cosmic Index mz RT Charge Peptide Peptide.id Mclust Medea
        // OutputFileType = 1 Ingest output file:   Index mz RT  0.0 Cluster#
        public static void AnalyzeDataFromFile(string fname)
        {
            char[] _sep = new[] { ' ', '\t' };
            int Maxcounts = 200;
            int[] CountLabels = new int[Maxcounts];
            for (int ChargeIndex = 0; ChargeIndex < Maxcounts; ChargeIndex++)
                CountLabels[ChargeIndex] = 0;
            int MinSplitSize = Program.ParameterVectorDimension + 3;
            int SplitPosition = 1 + Program.ParameterVectorDimension;
            if (Program.InputFileType == 1)
            {
                MinSplitSize = 5;
                SplitPosition = 4;
            }
            bool success = false;
            int count = 0;

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                        sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    string line;
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            if (splits.Length < MinSplitSize)
                            {
                                Exception e = DAVectorUtility.SALSAError("Count " + count.ToString() + "Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                throw (e);
                            }
                            int Charge = -101;
                            if (!Int32.TryParse(splits[SplitPosition], out Charge))
                                continue;
                            if ((Program.SelectedInputLabel < 0 ) && (Charge == -Program.SelectedInputLabel))
                                continue;
                            int position = 1 + Charge;
                            position = Math.Max(0, position);
                            position = Math.Min(Maxcounts-1, position);
                            ++CountLabels[position];
                            count++;
                        }
                    }
                    success = true;
                    if (Program.SelectedInputLabel >= 0)
                        DAVectorUtility.PointCount_Global = CountLabels[1 + Program.SelectedInputLabel];
                    else
                    {
                        DAVectorUtility.PointCount_Global = 0;
                        for (int selection = 0; selection < Maxcounts; selection++)
                            DAVectorUtility.PointCount_Global += CountLabels[selection];
                    }
                    DAVectorUtility.SALSAPrint(1, "File Analyzed " + fname + " Total " + count.ToString());
                    string label = "Charge";
                    if (Program.InputFileType == 1)
                        label = "Cluster";
                    if( CountLabels[0] > 0 )
                        DAVectorUtility.SALSAPrint(1, "Negative " + label + "s " + CountLabels[0].ToString());
                    for (int LabelIndex = 1; LabelIndex < (Maxcounts - 1); LabelIndex++)
                    {
                        if( CountLabels[LabelIndex] > 0 )
                            DAVectorUtility.SALSAPrint(1, label + " " + (LabelIndex - 1).ToString() + " Count " + CountLabels[LabelIndex].ToString());
                    }
                    if( CountLabels[Maxcounts-1] > 0 )
                        DAVectorUtility.SALSAPrint(1, "Overflow " + label + "s " + CountLabels[Maxcounts-1].ToString());
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Points data" + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File Analyze error " + fname);
                throw (e);
            }
        }   // End AnalyzeDataFromFile

        //  Readvectors = 0 Normal
        //  Readvectors = 1 Read points changing cluster labels
        //  Readvectors = 2 Read Cluster Centers
        public static void ReadDataFromFile(string fname, int ReadVectorsOption)
        {
            char[] _sep = new[] { ' ', '\t' };

            int FirstPointPosition = 0;
            int TotalNumberPointstoRead = 0;
            if (ReadVectorsOption <= 1)
            {
                FirstPointPosition = DAVectorUtility.PointStart_Process;
                TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;
            }
            int BeginLabel = 0;
            if (ReadVectorsOption == 1)
                BeginLabel = ParallelClustering.RunningSolution.Ncent_Global - 1;

            int MinSplitSize = Program.ParameterVectorDimension + 3;
            int SplitPosition = 1 + Program.ParameterVectorDimension;
            int LabelPosition = 4 + Program.ParameterVectorDimension;

            if (Program.InputFileType == 1)
            {
                MinSplitSize = 5;
                SplitPosition = 4;
                LabelPosition = 4;
            }

            bool success = false;
            string line = " Unset";
            int CountLinesinFile = 0;
            int CurrentLocalPointtoUse = 0;

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                    sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            if (splits.Length < MinSplitSize)
                            {
                                Exception e = DAVectorUtility.SALSAError("Count " + CountLinesinFile.ToString() + "Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                throw (e);
                            }
                            int Charge = -101;
                            if (!Int32.TryParse(splits[SplitPosition], out Charge))
                                continue;
                            if (Program.SelectedInputLabel >= 0)
                            {
                                if (Charge != Program.SelectedInputLabel)
                                    continue;
                            }
                            else
                            {
                                if (ReadVectorsOption>=2)
                                {
                                    if (Charge != -Program.SelectedInputLabel)
                                        continue;
                                }
                                else
                                {
                                    if (Charge == -Program.SelectedInputLabel)
                                        continue;
                                }
                            }
                            if ((ReadVectorsOption<=0) && (CountLinesinFile < FirstPointPosition))
                            {
                                CountLinesinFile += Program.Replicate;
                                continue;
                            }

                            int ActualPointPosition = 0;
                            if( ReadVectorsOption==0 )
                                ActualPointPosition = CountLinesinFile - FirstPointPosition;

                            else if (ReadVectorsOption == 1)
                            {
                                double CurrentLineY0 = double.Parse(splits[1]);
                                bool enditall = false;
                                bool endinputline = false;
                                while (true)
                                {
                                    if (CurrentLocalPointtoUse >= TotalNumberPointstoRead)
                                    {
                                        enditall = true;
                                        break;
                                    }
                                    double CurrentPointY0 = Program.PointPosition[CurrentLocalPointtoUse][0];
                                    if (CurrentLineY0 < (CurrentPointY0 - 0.000001))
                                    {
                                        endinputline = true;
                                        break;
                                    }
                                    if (CurrentLineY0 > (CurrentPointY0 + 0.000001))
                                    {
                                        ++CurrentLocalPointtoUse;
                                        continue;
                                    }
                                    else
                                    {
                                        if (Math.Abs(Program.PointPosition[CurrentLocalPointtoUse][1] - double.Parse(splits[2])) > 0.00001)
                                        {
                                            ++CurrentLocalPointtoUse;
                                            continue;
                                        } 
                                        if (Program.PointLabel[CurrentLocalPointtoUse] != 0)
                                        {
                                            Exception e = DAVectorUtility.SALSAError(fname + " Position " + (CurrentLocalPointtoUse + FirstPointPosition).ToString() 
                                                + " Inconsistent Cluster Number " + Program.PointLabel[CurrentLocalPointtoUse].ToString()
                                                + " Y0 " + Program.PointPosition[CurrentLocalPointtoUse][0].ToString("F5") + " " + line);
                                            throw (e);
                                        }
                                        ActualPointPosition = CurrentLocalPointtoUse;
                                        ++CurrentLocalPointtoUse;
                                        break;
                                    }
                                }
                                if (endinputline)
                                    continue;
                                if (enditall)
                                    break;
                            }
                            for (int CountReplicas = 0; CountReplicas < Program.Replicate; CountReplicas++)
                            {
                                double label;
                                if (ReadVectorsOption>=2)
                                {
                                    ParallelClustering.RunningSolution.Y_k_i_[ParallelClustering.RunningSolution.Ncent_Global][0] = double.Parse(splits[1]);
                                    ParallelClustering.RunningSolution.Y_k_i_[ParallelClustering.RunningSolution.Ncent_Global][1] = double.Parse(splits[2]);
                                    if (Program.ParameterVectorDimension > 2)
                                    {
                                        for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                                            ParallelClustering.RunningSolution.Y_k_i_[ParallelClustering.RunningSolution.Ncent_Global][VectorIndex] = double.Parse(splits[VectorIndex + 1]);
                                    }
                                    ++ParallelClustering.RunningSolution.Ncent_Global;
                                }
                                else
                                {
                                    if (ReadVectorsOption == 0)
                                    {
                                        Program.PointPosition[ActualPointPosition][0] = double.Parse(splits[1]);
                                        Program.PointPosition[ActualPointPosition][1] = double.Parse(splits[2]);
                                        Program.PointOriginalIndex[ActualPointPosition] = int.Parse(splits[0]);
                                        if (Program.ParameterVectorDimension > 2)
                                        {
                                            for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                                                Program.PointPosition[ActualPointPosition][VectorIndex] = double.Parse(splits[VectorIndex + 1]);
                                        }
                                    }
                                    if (!Double.TryParse(splits[LabelPosition], out label))
                                        label = 0.0;
                                    Program.PointLabel[ActualPointPosition] = (int)(label + 0.001);
                                    if (Program.PointLabel[ActualPointPosition] != 0)
                                        Program.PointLabel[ActualPointPosition] += BeginLabel;
                                }
                                ++ActualPointPosition;
                            }
                            CountLinesinFile += Program.Replicate;
                            if ((ReadVectorsOption<=0) && (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead)) )
                                break;
                            
                        }
                    }
                    if ( (ReadVectorsOption<=0) && (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead)) )
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank.ToString() 
                            + " Lines in File " + CountLinesinFile.ToString() + " Number to Read " + TotalNumberPointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Points data " + DAVectorUtility.MPI_Rank.ToString() + " " + CountLinesinFile.ToString() + " Start "
                    + FirstPointPosition.ToString() + " Number " + TotalNumberPointstoRead.ToString() + " " + line + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File read error " + fname);
                throw (e);
            }
        }   // End ReadDataFromFile

        public static void ReadDataFromFile(string fname, int ClusterPosition, int[] InitialPointAssignment, int StartPointPosition)
        {
            char[] _sep = new[] { ' ', ',', '\t' };

            int FirstPointPosition = 0;
            int TotalNumberPointstoRead = 0;
            FirstPointPosition = DAVectorUtility.PointStart_Process;
            TotalNumberPointstoRead = DAVectorUtility.PointCount_Process;

            int MinSplitSize = ClusterPosition + 1;
            if (StartPointPosition >= 0)
                MinSplitSize = Math.Max(MinSplitSize, StartPointPosition + Program.ParameterVectorDimension);

            bool success = false;
            string line = " Unset";
            int CountLinesinFile = 0;
            string stringtest = "";

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                    sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep, StringSplitOptions.RemoveEmptyEntries);
                            if (splits.Length < MinSplitSize)
                            {
                                DAVectorUtility.SALSAPrint(0, "Count " + CountLinesinFile.ToString() + " Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                continue;
                            }   // Skip header lines
                            double junk;
                            if (!Double.TryParse(splits[StartPointPosition], out junk))
                                continue;   // Skip header lines

                            if (CountLinesinFile < FirstPointPosition)
                            {
                                CountLinesinFile += 1;
                                continue;
                            }

                            int ActualPointPosition = CountLinesinFile - FirstPointPosition;
                            int label = 0;
                            if (StartPointPosition >= 0)
                            {
                                stringtest = "0 *" + splits[StartPointPosition];
                                Program.PointPosition[ActualPointPosition][0] = double.Parse(splits[StartPointPosition]);
                                stringtest = "1 *" + splits[StartPointPosition + 1];
                                Program.PointPosition[ActualPointPosition][1] = double.Parse(splits[StartPointPosition + 1]);
                                if (Program.ParameterVectorDimension > 2)
                                {
                                    for (int VectorIndex = 2; VectorIndex < Program.ParameterVectorDimension; VectorIndex++)
                                    {
                                        stringtest = VectorIndex.ToString() + " *" + splits[StartPointPosition + VectorIndex];
                                        Program.PointPosition[ActualPointPosition][VectorIndex] = double.Parse(splits[VectorIndex + StartPointPosition]);
                                    }
                                }
                            }
                            if (ClusterPosition >= 0)
                            {
                                if (!Int32.TryParse(splits[ClusterPosition], out label))
                                    label = 1;
                                InitialPointAssignment[ActualPointPosition] = label - 1;
                            }

                            ++ActualPointPosition;
                            ++CountLinesinFile;
                            if (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead))
                                break;

                        }
                    }
                    if (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead))
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank.ToString()
                            + " Lines in File " + CountLinesinFile.ToString() + " Number to Read " + TotalNumberPointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine(stringtest + "* Failed reading Points data " + DAVectorUtility.MPI_Rank.ToString() + " " + CountLinesinFile.ToString() + " Start "
                    + FirstPointPosition.ToString() + " Number " + TotalNumberPointstoRead.ToString() + " " + line + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File read error " + fname);
                throw (e);
            }
        }   // End ReadDataFromFile  

        public static void Read3DDataFromFile(string fname, int VectorSize, int StartPointPosition)
        {
            char[] _sep = new[] { ' ', ',', '\t' };

            int FirstPointPosition = 0;
            int TotalNumberPointstoRead = 0;
            TotalNumberPointstoRead = DAVectorUtility.PointCount_Global;

            int MinSplitSize = StartPointPosition + VectorSize;

            bool success = false;
            string line = " Unset";
            int CountLinesinFile = 0;

            try
            {
                StreamReader sr = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    Stream stream = File.Open(fname, FileMode.Open, FileAccess.Read, FileShare.Read);
                    sr = new StreamReader(stream);
                }
                if (sr != null)
                {
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            if (splits.Length < MinSplitSize)
                            {
                                DAVectorUtility.SALSAPrint(0, "Count " + CountLinesinFile.ToString() + " Illegal data length on Point file " + splits.Length.ToString()
                                    + " " + MinSplitSize.ToString() + " " + line);
                                continue;
                            }   // Skip header lines
                            double junk;
                            if (!Double.TryParse(splits[StartPointPosition], out junk))
                                continue;   // Skip header lines

                            if (CountLinesinFile < FirstPointPosition)
                            {
                                CountLinesinFile += 1;
                                continue;
                            }

                            int ActualPointPosition = CountLinesinFile - FirstPointPosition;
                            if (StartPointPosition >= 0)
                            {
                                Program.FullPoint3DPosition[ActualPointPosition][0] = double.Parse(splits[StartPointPosition]);
                                Program.FullPoint3DPosition[ActualPointPosition][1] = double.Parse(splits[StartPointPosition + 1]);
                                if (VectorSize > 2)
                                {
                                    for (int VectorIndex = 2; VectorIndex < VectorSize; VectorIndex++)
                                        Program.FullPoint3DPosition[ActualPointPosition][VectorIndex] = double.Parse(splits[VectorIndex + StartPointPosition]);
                                }
                            }

                            ++ActualPointPosition;
                            ++CountLinesinFile;
                            if (CountLinesinFile >= (FirstPointPosition + TotalNumberPointstoRead))
                                break;

                        }
                    }
                    if (CountLinesinFile != (FirstPointPosition + TotalNumberPointstoRead))
                    {
                        Exception e = DAVectorUtility.SALSAError("Illegal count on Points file " + fname + " Rank " + DAVectorUtility.MPI_Rank.ToString()
                            + " Lines in File " + CountLinesinFile.ToString() + " Number to Read " + TotalNumberPointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Points data " + DAVectorUtility.MPI_Rank.ToString() + " " + CountLinesinFile.ToString() + " Start "
                    + FirstPointPosition.ToString() + " Number " + TotalNumberPointstoRead.ToString() + " " + line + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = DAVectorUtility.SALSAError("DA Vector File read error " + fname);
                throw (e);
            }
        }   // End Read3DDataFromFile  

    }   // End class DAVectorReadData
}
