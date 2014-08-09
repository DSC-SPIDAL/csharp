using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MPI;

using SALSALibrary;

namespace Salsa.PairwiseClusteringTPL
{
    public class FindCenters
    {
        // Sequence Space Measures
        //  SmallestDistanceMeans Find Sequence that has Minimum Mean Distance to all other sequences in cluster using Smith Waterman Gotoh distance
        //  Bucket Method defaults f= 0.15,0.4,0.75 for buckets 0 1 2
        //  Find value V(f) such that fraction f of all pairwise distances have size < V(f)
        //  Find Center with maximum number of points satisfying distance < V(f) for distances to center

        //  MDS Space
        //  SmallestMDSDistanceMeans Find Sequence that has Minimum Mean Distance to all other sequences in cluster using distances from mapping to 3D
        //  SmallestMDSCoG Take cluster points in 3D. Find geometric center with position gotten by averaging all positions in cluster. Find sequence nearest the "center of gravity"
        
        
        //  GroupIndex[GlobalPointIndex] is Cluster (group) number with label starting at StartPosition
        //  NumberofGroups is Number of clusters (families)
        //  PWCUtility.NumberofCenters is Number of Centers of each type to be found
        //  StartPosition is lowest cluster number 0 or 1
        //  CenterOutputFileName is file name of output file containing centers
        //  MDSInputFileName is file containg MDS  data
        //  Ignore group if number of points in group <= PWCUtility.NumberofCenters

        public static double[][] MDSvalues; // Array to store MDS positions if available

        public static void FindGroupCenters(int[] GroupIndex, int NumberofGroups, int StartPosition, string CenterOutputFileName, string MDSInputFileName, string LabelInputFileName)
        {
            if (PWCUtility.NumberofCenters <= 0)
                return;

            // for Bucket Method
            int NumberofBins = 1000;

            // For total evaluation
            int numberofcategories = 3 + PWCUtility.NumberofBuckets;
            int[] CombinedListIndices = new int[numberofcategories * PWCUtility.NumberofCenters];
            double[] CombinedListRatings = new double[numberofcategories * PWCUtility.NumberofCenters];
            double[] CombinedListRatingsCount = new double[numberofcategories * PWCUtility.NumberofCenters];
            bool[] CombinedListSourceSeq = new bool[numberofcategories * PWCUtility.NumberofCenters];
            bool[] CombinedListSourceMDS = new bool[numberofcategories * PWCUtility.NumberofCenters];
            double[] TopMeanDistance = new double[PWCUtility.NumberofCenters];
            double[] TopMDSMeanDistance = new double[PWCUtility.NumberofCenters];
            double[] TopMDSCoGDistance = new double[PWCUtility.NumberofCenters];
            double[] TopBucketDistance = new double[PWCUtility.NumberofCenters];

            // For Print Out
            int[] UtilityIndex = new int[numberofcategories * PWCUtility.NumberofCenters];
            string[] PointProperties = new string[numberofcategories * PWCUtility.NumberofCenters];

            double[][] GroupSmallestfromMinMeans = new double[NumberofGroups][];   // List means from of centers from min means
            int[][] GroupIndexfromMinMeans = new int[NumberofGroups][];            //  List of Global Point indices for centers for min means
            double[] GlobalGroupMean = new double[NumberofGroups];          // Mean distance in groups
            double[] GlobalGroupMax = new double[NumberofGroups];           // Maximum distances in groups
            int[] GroupCount = new int[NumberofGroups];                     // Count of points in a group

            // For MDS need to read all points into all processes
            MDSvalues = new double[PWCUtility.PointCount_Global][];
            if (PWCUtility.addMDS > 0)
            {
                for (int GlobalPointIndex = 0; GlobalPointIndex < PWCUtility.PointCount_Global; GlobalPointIndex++)
                {
                    MDSvalues[GlobalPointIndex] = new double[3];
                }
                ReadMDS(MDSInputFileName, MDSvalues, 0, PWCUtility.PointCount_Global);
            }

            // Read labels if they exist
            bool LabelsAvailable = false;
            string[] SequenceLabels = new string[PWCUtility.PointCount_Global];
            int[] SequenceLengths = new int[PWCUtility.PointCount_Global];
            if (LabelInputFileName.Length > 0)
            {
                LabelsAvailable = true;
                if (PWCUtility.MPI_Rank == 0)
                {
                    ReadLabels(LabelInputFileName, SequenceLabels, SequenceLengths, 0, PWCUtility.PointCount_Global);
                    int maxlen = 0;
                    for (int looplabel = 0; looplabel < PWCUtility.PointCount_Global; looplabel++)
                        maxlen = Math.Max(maxlen, SequenceLabels[looplabel].Length);
                    for (int looplabel = 0; looplabel < PWCUtility.PointCount_Global; looplabel++)
                    {
                        for (int addblanks = SequenceLabels[looplabel].Length; addblanks < maxlen; addblanks++)
                            SequenceLabels[looplabel] += " ";
                    }
                }
                PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                PWCUtility.MPI_communicator.Broadcast<int>(ref SequenceLengths, 0); ;
                PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
            }
            else
            {
                PWCUtility.LengthCut1 = -1;
                PWCUtility.LengthCut2 = -1;
            }

            double[][] GroupMDSCoG = new double[NumberofGroups][];    // Mean Position(CoG) of Groups in MDS space
            double[] GlobalGroupMDSMean = new double[NumberofGroups];   // Mean of Group in MDS space

            int[][] GroupIndexfromMDSMinMeans = new int[NumberofGroups][];     // Indices from min mean in MDS space
            double[][] GroupSmallestfromMDSMinMeans = new double[NumberofGroups][];    // min means from min mean in MDS space
            
            int[][] GroupIndexfromMDSCoG = new int[NumberofGroups][];   // Indices nearest Center of Mass in MDS space
            double[][] GroupSmallestfromMDSCoG = new double[NumberofGroups][];  // Distances from Center of Mass in MDS space

            for (int group = 0; group < NumberofGroups; group++)
            {
                GroupCount[group] = 0;
                GlobalGroupMean[group] = 0.0;
                GlobalGroupMax[group] = 0.0;
                GroupSmallestfromMinMeans[group] = new double[PWCUtility.NumberofCenters];
                GroupIndexfromMinMeans[group] = new int[PWCUtility.NumberofCenters];

                GlobalGroupMDSMean[group] = 0.0;
                GroupMDSCoG[group] = new double[3];
                GroupIndexfromMDSMinMeans[group] = new int[PWCUtility.NumberofCenters];
                GroupSmallestfromMDSMinMeans[group] = new double[PWCUtility.NumberofCenters];
                GroupIndexfromMDSCoG[group] = new int[PWCUtility.NumberofCenters];
                GroupSmallestfromMDSCoG[group] = new double[PWCUtility.NumberofCenters];
            }

            for (int GlobalPointindex = 0; GlobalPointindex < PWCUtility.PointCount_Global; GlobalPointindex++)
            {
                int group = GroupIndex[GlobalPointindex] - StartPosition;
                ++GroupCount[group];
            }
            for (int group = 0; group < NumberofGroups; group++)
            {
                if (GroupCount[group] <= PWCUtility.NumberofCenters)
                    GroupCount[group] = 0;
            }

            // Initialize accumulation data structures for each group
            GlobalReductions.FindDoubleMax[] FindGroupMax = new GlobalReductions.FindDoubleMax[NumberofGroups];
            GlobalReductions.FindMeanSigma[] FindGroupMeansigma = new GlobalReductions.FindMeanSigma[NumberofGroups];
            GlobalReductions.FindMeanSigma[] FindGroupMDSMeansigma = new GlobalReductions.FindMeanSigma[NumberofGroups];
            GlobalReductions.FindMeanSigma[][] FindGroupMDSCoG = new GlobalReductions.FindMeanSigma[NumberofGroups][];

            GlobalReductions.FindManyMinValuewithIndex[] FindGroupOriginalDataCenters_mean = new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];
            GlobalReductions.FindManyMinValuewithIndex[] FindGroupMDSCenters_mean = new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];
            GlobalReductions.FindManyMinValuewithIndex[] FindGroupMDSCenters_CoG = new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];

            for (int group = 0; group < NumberofGroups; group++)
            {
                if (GroupCount[group] == 0)
                    continue;

                int Countlinks = 0;
                if (LabelsAvailable)
                {
                    for (int GlobalPointIndex = 0; GlobalPointIndex < PWCUtility.PointCount_Global; GlobalPointIndex++)
                    {
                        int group_Point = GroupIndex[GlobalPointIndex] - StartPosition;
                        if (group != group_Point)
                            continue;
                        if ((SequenceLengths[GlobalPointIndex] > PWCUtility.LengthCut1) && (SequenceLengths[GlobalPointIndex] > PWCUtility.LengthCut2))
                            ++Countlinks;

                    }
                }
                else Countlinks = GroupCount[group];
                if (Countlinks < PWCUtility.LinkCountinCenterFinding)
                {
                    GroupCount[group] = 0;
                    continue;
                }

                FindGroupMax[group] = new GlobalReductions.FindDoubleMax(PWCUtility.ThreadCount);
                FindGroupMeansigma[group] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
                FindGroupMDSMeansigma[group] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
                FindGroupMDSCoG[group] = new GlobalReductions.FindMeanSigma[3];
                for (int MDSIndex = 0; MDSIndex < 3; MDSIndex++)
                {
                    FindGroupMDSCoG[group][MDSIndex] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
                }

                FindGroupOriginalDataCenters_mean[group] = new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
                FindGroupMDSCenters_mean[group] = new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
                FindGroupMDSCenters_CoG[group] = new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
            }



            //  Loop over points GlobalPointIndex1 to find MDS means which are needed for next step
            if (PWCUtility.addMDS > 0)
            {
                Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                {
                    int indexlen = PWCUtility.PointsperThread[ThreadNo];
                    int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;

                    for (int index = beginpoint; index < indexlen + beginpoint; index++)
                    {
                        int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;

                        int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                        if ((group1 < 0) || (group1 >= NumberofGroups))
                        {
                            Exception e = PWCUtility.SALSAError(" Illegal group number " + group1.ToString() + " Point " + GlobalPointIndex1.ToString());
                            throw (e);
                        }
                        if (GroupCount[group1] <= 0)
                            continue;

                        for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global; GlobalPointIndex2++)
                        {
                            if (GlobalPointIndex1 == GlobalPointIndex2)
                                continue;
                            int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                            if ((group2 < 0) || (group2 >= NumberofGroups))
                            {
                                Exception e = PWCUtility.SALSAError(" Illegal group number " + group2.ToString() + " Point " + GlobalPointIndex2.ToString());
                                throw (e);
                            }
                            if (group1 != group2)
                                continue;

                            double tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                            FindGroupMDSMeansigma[group1].addapoint(ThreadNo, tmp);
                        }   // End Loop over GlobalPointIndex2

                        for (int MDSindex = 0; MDSindex < 3; MDSindex++)
                            FindGroupMDSCoG[group1][MDSindex].addapoint(ThreadNo, MDSvalues[GlobalPointIndex1][MDSindex]);
                    }
                });  // End Parallel Section defining MDS average values

                for (int group = 0; group < NumberofGroups; group++)
                {
                    if (GroupCount[group] == 0)
                        continue;
                    FindGroupMDSMeansigma[group].sumoverthreadsandmpi();
                    GlobalGroupMDSMean[group] = FindGroupMDSMeansigma[group].Totalmean;
                    for (int MDSindex = 0; MDSindex < 3; MDSindex++)
                    {
                        FindGroupMDSCoG[group][MDSindex].sumoverthreadsandmpi();
                        GroupMDSCoG[group][MDSindex] = FindGroupMDSCoG[group][MDSindex].Totalmean;
                    }

                }
            }   // End case where MDS values exist

            //  Loop over points GlobalPointIndex1
            Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
            {
                int indexlen = PWCUtility.PointsperThread[ThreadNo];
                int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;

                for (int index = beginpoint; index < indexlen + beginpoint; index++)
                {
                    int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;

                    int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                    if ((group1 < 0) || (group1 >= NumberofGroups))
                    {
                        Exception e = PWCUtility.SALSAError(" Illegal group number " + group1.ToString() + " Point " + GlobalPointIndex1.ToString());
                        throw(e);
                    }
                    if (GroupCount[group1] <= 0)
                        continue;
                    int LinkCutUsed = PWCUtility.LinkCountinCenterFinding;
                    LinkCutUsed = Math.Min(LinkCutUsed, GroupCount[group1] - 20);
                    LinkCutUsed = Math.Max(LinkCutUsed, 1);

                    double thispointmean = 0.0;
                    double thispointMDSmean = 0.0;
                    int Countlinks = 0;
                    for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global; GlobalPointIndex2++)
                    {
                        if (GlobalPointIndex1 == GlobalPointIndex2)
                            continue;
                        int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                        if ((group2 < 0) || (group2 >= NumberofGroups))
                        {
                            Exception e = PWCUtility.SALSAError(" Illegal group number " + group2.ToString() + " Point " + GlobalPointIndex2.ToString());
                            throw (e);
                        }
                        if (group1 != group2)
                            continue;

                        double tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                        if ( tmp > PWCUtility.MinimumDistanceCut )
                        {
                            if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2) )
                            {
                                ++Countlinks;
                                thispointmean += tmp;
                                FindGroupMax[group1].addapoint(ThreadNo, tmp);
                                FindGroupMeansigma[group1].addapoint(ThreadNo, tmp);
                            }
                        }
                        if (PWCUtility.addMDS > 0)
                        {
                            tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                            thispointMDSmean += tmp;
                        }

                    }   // End Loop over GlobalPointIndex2
                    if (Countlinks >= LinkCutUsed) 
                        FindGroupOriginalDataCenters_mean[group1].addapoint(ThreadNo, GlobalPointIndex1, thispointmean/Countlinks);

                    if (PWCUtility.addMDS > 0)
                    {
                        double MDSmeandistance = getMDSDistancefromPoint(GlobalPointIndex1, GroupMDSCoG[group1]);
                        FindGroupMDSCenters_CoG[group1].addapoint(ThreadNo, GlobalPointIndex1, MDSmeandistance);
                        if (Countlinks > 0)
                            FindGroupMDSCenters_mean[group1].addapoint(ThreadNo, GlobalPointIndex1, thispointMDSmean / Countlinks); 
                    }
                }
            });  // End Parallel Section defining GlobalPointIndex1

            // Finishup finding centers from sequence distance means and MDS CoG and means
            for (int group = 0; group < NumberofGroups; group++)
            {
                if (GroupCount[group] <= 0)
                    continue;
                FindGroupMax[group].sumoverthreadsandmpi();
                FindGroupMeansigma[group].sumoverthreadsandmpi();
                GlobalGroupMean[group] = FindGroupMeansigma[group].Totalmean;
                GlobalGroupMax[group] = FindGroupMax[group].TotalMax;

                FindGroupOriginalDataCenters_mean[group].sumoverthreadsandmpi();
                int TotalFound = (int) (FindGroupOriginalDataCenters_mean[group].TotalNumberofPoints + 0.01);
                PWCUtility.SALSAPrint(0, "Group " + group.ToString() + " Number in Mean " + GroupCount[group].ToString() + " After Cuts " + TotalFound.ToString());

                if (PWCUtility.addMDS > 0)
                {
                    FindGroupMDSCenters_CoG[group].sumoverthreadsandmpi();
                    FindGroupMDSCenters_mean[group].sumoverthreadsandmpi();
                }

                for (int CenterIndex = 0; CenterIndex < PWCUtility.NumberofCenters; CenterIndex++)
                {
                    GroupSmallestfromMinMeans[group][CenterIndex] = FindGroupOriginalDataCenters_mean[group].OrderedMinValue[CenterIndex];
                    GroupIndexfromMinMeans[group][CenterIndex] = FindGroupOriginalDataCenters_mean[group].OrderedIndexValue[CenterIndex];

                    if (PWCUtility.addMDS > 0)
                    {
                        GroupIndexfromMDSMinMeans[group][CenterIndex] = FindGroupMDSCenters_mean[group].OrderedIndexValue[CenterIndex];     // Indices from min mean in MDS space
                        GroupSmallestfromMDSMinMeans[group][CenterIndex] = FindGroupMDSCenters_mean[group].OrderedMinValue[CenterIndex];    // min means from min mean in MDS space

                        GroupIndexfromMDSCoG[group][CenterIndex] = FindGroupMDSCenters_CoG[group].OrderedIndexValue[CenterIndex];   // Indices nearest Center of Mass in MDS space
                        GroupSmallestfromMDSCoG[group][CenterIndex] = FindGroupMDSCenters_CoG[group].OrderedMinValue[CenterIndex];  // Distances from group CoG
                    }
                }
            }

            // Loop over groups -- output Min means MDS and bucket solutions -- also find bucket answers
            for (int group = 0; group < NumberofGroups; group++)
            {
                if (GroupCount[group] <= 0)
                    continue;

                //  Overall header
                PWCUtility.SALSAPrint(1, "\n++++++++++++++++++\nFind Centers Group=" + (group + StartPosition).ToString() + " Count " + GroupCount[group].ToString()
                    + " Mean " + GlobalGroupMean[group].ToString("F4") + " MDS Mean " + GlobalGroupMDSMean[group].ToString("F4") + " Max " + GlobalGroupMax[group].ToString("F4") );
                PWCUtility.SALSAPrint(1, "MDS Center of Gravity " + GroupMDSCoG[group][0].ToString("F4") + " " + GroupMDSCoG[group][1].ToString("F4") + " " + GroupMDSCoG[group][2].ToString("F4"));

                // Output conventional min Means with weight of 1
                int countindicesfound = 0;
                TopMeanDistance[0] = 0.0;

                //  Output min means from sequence space
                int NumMinMeansFound = 0;
                for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++)
                {
                    int GlobalPointIndex1 = GroupIndexfromMinMeans[group][CenterIndex1];
                    if (GlobalPointIndex1 < 0)
                    {
                        PWCUtility.SALSAPrint(0,"Group " + group.ToString() + " Center Index for Sequence Min-Mean " + CenterIndex1.ToString() + " Undefined");
                        continue;
                    }
                    ++NumMinMeansFound;

                    // Add into Global Rating List
                    int useposition = countindicesfound;
                    if (countindicesfound > 0)
                    {
                        for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++)
                        {
                            if (CombinedListIndices[loopindexlist] != GlobalPointIndex1)
                                continue;
                            useposition = loopindexlist;
                            break;
                        }
                    }
                    if (useposition < countindicesfound) {
                        CombinedListRatingsCount[useposition] += 1.0;
                        CombinedListRatings[useposition] +=  CenterIndex1;
                        CombinedListSourceSeq[useposition] = true;
                    }
                    else {
                        CombinedListIndices[useposition] = GlobalPointIndex1;
                        CombinedListRatingsCount[useposition] = 1.0;
                        CombinedListRatings[useposition] = CenterIndex1;
                        CombinedListSourceSeq[useposition] = true;
                        CombinedListSourceMDS[useposition] = false;
                        ++countindicesfound;
                    }

                    //  Calculate Averages and list each choice
                    double avgmean = 0.0;
                    double Avg_CoGcenters = 0.0;
                    int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                    if (owner == PWCUtility.MPI_Rank)
                    {
                        for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++)
                        {
                            int GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                            if (GlobalPointIndex2 < 0)
                            {
                                Exception e = PWCUtility.SALSAError("Means Group-2 " + group.ToString() + " Center " + CenterIndex2.ToString() + " Undefined");
                                throw (e);
                            }
                            double tmp;
                            if (CenterIndex1 != CenterIndex2)
                            {
                                tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (CenterIndex1 == 0)
                                    TopMeanDistance[CenterIndex2] = tmp;
                                avgmean = avgmean + tmp;
                            }
                            if (PWCUtility.addMDS > 0)
                            {
                                GlobalPointIndex2 = GroupIndexfromMDSCoG[group][CenterIndex2];
                                tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (CenterIndex2 == 0)
                                    TopMDSCoGDistance[CenterIndex1] = tmp;
                                Avg_CoGcenters += tmp;
                            }
                        }
                        avgmean = avgmean / (PWCUtility.NumberofCenters - 1);
                        if (PWCUtility.addMDS > 0)
                            Avg_CoGcenters = Avg_CoGcenters / PWCUtility.NumberofCenters;
                    }
                    if (PWCUtility.MPI_Size > 1)
                    {
                        PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                        PWCUtility.MPI_communicator.Broadcast<double>(ref avgmean, owner);
                        if (PWCUtility.addMDS > 0)
                        {
                            PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_CoGcenters, owner);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref TopMDSCoGDistance[CenterIndex1], owner);
                        }
                        if (CenterIndex1 == 0)
                            PWCUtility.MPI_communicator.Broadcast<double>(ref TopMeanDistance, owner);
                        PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                    }

                    if (PWCUtility.MPI_Rank == 0)
                    {
                        string Seqlabel = "";
                        if (LabelsAvailable)
                        {
                            Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" + SequenceLengths[GlobalPointIndex1].ToString();
                        }
                        PointProperties[CenterIndex1] = "Measure=" + GroupSmallestfromMinMeans[group][CenterIndex1].ToString("F4")
                            + " Method=SmallestDistanceMeans Group=" + (group + StartPosition).ToString() + Seqlabel
                             + " Date=\"" + DateTime.Now.ToLocalTime() + "\"";
                        string MDSstring = "";
                        if (PWCUtility.addMDS > 0)
                            MDSstring = " Avg MDS Distce CoG's " + Avg_CoGcenters.ToString("F4") + " to top CoG " + TopMDSCoGDistance[CenterIndex1].ToString("F4");
                        PWCUtility.SALSAPrint(1, "Mean " + GroupSmallestfromMinMeans[group][CenterIndex1].ToString("F4")
                            + " Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Avg Distce-means " + avgmean.ToString("E3")
                            + " Distce to top mean " + TopMeanDistance[CenterIndex1].ToString("E3") + MDSstring + Seqlabel);
                    }
                }

                if (PWCUtility.MPI_Rank == 0)
                {
                    if(NumMinMeansFound > 0)
                    {
                        for (int CenterIndex = 0; CenterIndex < NumMinMeansFound; CenterIndex++)
                        {
                            UtilityIndex[CenterIndex] = GroupIndexfromMinMeans[group][CenterIndex];
                        }
                        WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, NumMinMeansFound, true);
                    }
                }

                //  Output MDS results -- min mean method
                if (PWCUtility.addMDS > 0)
                {
                    double scaleMDSratings = 1.0;   // Use to scale ratings
                    TopMDSMeanDistance[0] = 0.0;
                    int NumMDSMeansFound = 0;
                    for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++)
                    {
                        int GlobalPointIndex1 = GroupIndexfromMDSMinMeans[group][CenterIndex1];
                        if (GlobalPointIndex1 < 0)
                        {
                            PWCUtility.SALSAPrint(0, "MDS Means Group " + group.ToString() + " Center " + CenterIndex1.ToString() + " Undefined");
                            continue;
                        }
                        NumMDSMeansFound++;
                        // Add into Global Rating List
                        int useposition = countindicesfound;
                        if (countindicesfound > 0)
                        {
                            for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++)
                            {
                                if (CombinedListIndices[loopindexlist] != GlobalPointIndex1)
                                    continue;
                                useposition = loopindexlist;
                                break;
                            }
                        }
                        if (useposition < countindicesfound)
                        {
                            CombinedListRatingsCount[useposition] += scaleMDSratings;
                            CombinedListRatings[useposition] += scaleMDSratings * CenterIndex1;
                            CombinedListSourceMDS[useposition] = true;
                        }
                        else
                        {
                            CombinedListIndices[useposition] = GlobalPointIndex1;
                            CombinedListRatingsCount[useposition] = scaleMDSratings;
                            CombinedListRatings[useposition] = scaleMDSratings * CenterIndex1;
                            CombinedListSourceSeq[useposition] = false;
                            CombinedListSourceMDS[useposition] = true;
                            ++countindicesfound;
                        }

                        double Avg_MeanMethod = 0.0;
                        double Avg_Internal = 0.0;
                        int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                        if (owner == PWCUtility.MPI_Rank)
                        {

                            for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++)
                            {
                                int GlobalPointIndex2 = GroupIndexfromMDSMinMeans[group][CenterIndex2];
                                if (GlobalPointIndex2 < 0)
                                {
                                    Exception e = PWCUtility.SALSAError("MDS Means Group-2 " + group.ToString() + " Center " + CenterIndex2.ToString() + " Undefined");
                                    throw (e);
                                }
                                double tmp;
                                if (CenterIndex1 != CenterIndex2)
                                {
                                    tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                    if (CenterIndex1 == 0)
                                        TopMDSMeanDistance[CenterIndex2] = tmp;
                                    Avg_Internal = Avg_Internal + tmp;
                                }
                                GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                                tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                Avg_MeanMethod = Avg_MeanMethod + tmp;
                                if (CenterIndex2 == 0)
                                    TopMeanDistance[CenterIndex1] = tmp;
                            }

                        }

                        if (PWCUtility.MPI_Size > 1)
                        {
                            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_Internal, owner);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_MeanMethod, owner);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref TopMeanDistance[CenterIndex1], owner);
                            if (CenterIndex1 == 0)
                                PWCUtility.MPI_communicator.Broadcast<double>(ref TopMDSMeanDistance, owner);
                            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                        }
                        Avg_MeanMethod = Avg_MeanMethod / PWCUtility.NumberofCenters;
                        Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);

                        if (PWCUtility.MPI_Rank == 0)
                        {
                            string Seqlabel = "";
                            if (LabelsAvailable)
                            {
                                Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" + SequenceLengths[GlobalPointIndex1].ToString();
                            }
                            PointProperties[CenterIndex1] = "Measure=" + GroupSmallestfromMDSMinMeans[group][CenterIndex1].ToString("F4")
                                + " Method=SmallestMDSDistanceMeans Group=" + (group + StartPosition).ToString() + Seqlabel
                                 + " Date=\"" + DateTime.Now.ToLocalTime() + "\"";
                            PWCUtility.SALSAPrint(1, "MDS Mean " + GroupSmallestfromMDSMinMeans[group][CenterIndex1].ToString("F4")
                                + " Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-MDSmeans " + Avg_Internal.ToString("E3")
                                + " Distce to MDS Mean top " + TopMDSMeanDistance[CenterIndex1].ToString("E3")
                                + " Distce-OriginalMeans " + Avg_MeanMethod.ToString("E3") + " Distce Original Mean Top " + TopMeanDistance[CenterIndex1].ToString("E3") + Seqlabel);
                        }
                    }

                    if (PWCUtility.MPI_Rank == 0)
                    {
                        if (NumMDSMeansFound > 0)
                        {
                            for (int CenterIndex = 0; CenterIndex < NumMDSMeansFound; CenterIndex++)
                            {
                                UtilityIndex[CenterIndex] = GroupIndexfromMDSMinMeans[group][CenterIndex];
                            }
                            WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, NumMDSMeansFound, true);
                        }
                    }

                    //  Output MDS results -- CoG method
                    TopMDSCoGDistance[0] = 0.0;
                    for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++)
                    {
                        int GlobalPointIndex1 = GroupIndexfromMDSCoG[group][CenterIndex1];
                        if (GlobalPointIndex1 < 0)
                        {
                            Exception e = PWCUtility.SALSAError("MDS CoG Group-1 " + group.ToString() + " Center " + CenterIndex1.ToString() + " Undefined");
                            throw (e);
                        }

                        // Add into Global Rating List
                        int useposition = countindicesfound;
                        if (countindicesfound > 0)
                        {
                            for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++)
                            {
                                if (CombinedListIndices[loopindexlist] != GlobalPointIndex1)
                                    continue;
                                useposition = loopindexlist;
                                break;
                            }
                        }
                        if (useposition < countindicesfound)
                        {
                            CombinedListRatingsCount[useposition] += scaleMDSratings;
                            CombinedListRatings[useposition] += scaleMDSratings * CenterIndex1;
                            CombinedListSourceMDS[useposition] = true;
                        }
                        else
                        {
                            CombinedListIndices[useposition] = GlobalPointIndex1;
                            CombinedListRatingsCount[useposition] = scaleMDSratings;
                            CombinedListRatings[useposition] = scaleMDSratings * CenterIndex1;
                            CombinedListSourceSeq[useposition] = false;
                            CombinedListSourceMDS[useposition] = true;
                            ++countindicesfound;
                        }

                        double Avg_MeanMethod = 0.0;
                        int Num_MeanMethod = 0;
                        double Avg_Internal = 0.0;
                        int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                        if (owner == PWCUtility.MPI_Rank)
                        {

                            for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++)
                            {
                                int GlobalPointIndex2 = GroupIndexfromMDSCoG[group][CenterIndex2];
                                if (GlobalPointIndex2 < 0)
                                {
                                    Exception e = PWCUtility.SALSAError("MDS CoG Group-2 " + group.ToString() + " Center " + CenterIndex2.ToString() + " Undefined");
                                    throw (e);
                                }
                                double tmp;
                                if (CenterIndex1 != CenterIndex2)
                                {
                                    tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                    if (CenterIndex1 == 0)
                                        TopMDSCoGDistance[CenterIndex2] = tmp;
                                    Avg_Internal = Avg_Internal + tmp;
                                }
                                GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                                if (GlobalPointIndex2 < 0)
                                {
                                    if (CenterIndex2 == 0)
                                        TopMeanDistance[CenterIndex1] = 0.0;
                                }
                                else
                                {
                                    tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                    Avg_MeanMethod = Avg_MeanMethod + tmp;
                                    ++Num_MeanMethod;
                                    if (CenterIndex2 == 0)
                                        TopMeanDistance[CenterIndex1] = tmp;
                                }
                            }

                        }

                        if (PWCUtility.MPI_Size > 1)
                        {
                            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_Internal, owner);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_MeanMethod, owner);
                            PWCUtility.MPI_communicator.Broadcast<int>(ref Num_MeanMethod, owner);
                            PWCUtility.MPI_communicator.Broadcast<double>(ref TopMeanDistance[CenterIndex1], owner);
                            if (CenterIndex1 == 0)
                                PWCUtility.MPI_communicator.Broadcast<double>(ref TopMDSCoGDistance, owner);
                            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                        } 
                        if(Num_MeanMethod > 0)
                            Avg_MeanMethod = Avg_MeanMethod / Num_MeanMethod;
                        Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);

                        if (PWCUtility.MPI_Rank == 0)
                        {
                            string Seqlabel = "";
                            if (LabelsAvailable)
                            {
                                Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" + SequenceLengths[GlobalPointIndex1].ToString();
                            }
                            PointProperties[CenterIndex1] = "Measure=" + GroupSmallestfromMDSCoG[group][CenterIndex1].ToString("F4")
                                + " Method=SmallestMDSCoG Group=" + (group + StartPosition).ToString() + Seqlabel
                                 + " Date=\"" + DateTime.Now.ToLocalTime() + "\"";
                            PWCUtility.SALSAPrint(1, "MDS CoG  " + GroupSmallestfromMDSCoG[group][CenterIndex1].ToString("F4")
                                + " Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-MDS CoG " + Avg_Internal.ToString("E3")
                                + " Distce to MDS CoG top " + TopMDSCoGDistance[CenterIndex1].ToString("E3")
                                + " Distce-OriginalMeans " + Avg_MeanMethod.ToString("E3") + " Distce Original Mean Top " + TopMeanDistance[CenterIndex1].ToString("E3") + Seqlabel);
                        }
                    }

                    if (PWCUtility.MPI_Rank == 0)
                    {
                        for (int CenterIndex = 0; CenterIndex < PWCUtility.NumberofCenters; CenterIndex++)
                        {
                            UtilityIndex[CenterIndex] = GroupIndexfromMDSCoG[group][CenterIndex];
                        }
                        WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, PWCUtility.NumberofCenters, true);
                    }

                }   // End Output of MDS

                // Do Buckets *******************************************************************************************
                if (PWCUtility.NumberofBuckets > 0)
                {
                    double[] BucketRadii = new double[PWCUtility.NumberofBuckets];

                    //  First Histogram distance values so we can convert bucket fractions into radii
                    double fudge = (double)NumberofBins / GlobalGroupMax[group];
                    GlobalReductions.FindDoubleArraySum DistanceHistogramBinCounts = new GlobalReductions.FindDoubleArraySum(PWCUtility.ThreadCount, NumberofBins);

                    //  Loop over points selecting those in this group
                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;
                        DistanceHistogramBinCounts.startthread(ThreadNo);

                        for (int index = beginpoint; index < indexlen + beginpoint; index++)
                        {
                            int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;

                            int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                            if (group1 != group)
                                continue;
                            for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global; GlobalPointIndex2++)
                            {
                                if (GlobalPointIndex1 == GlobalPointIndex2)
                                    continue;
                                int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                                if (group2 != group)
                                    continue;
                                double tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (tmp > PWCUtility.MinimumDistanceCut)
                                {
                                    if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2))
                                    {
                                        int itmp = (int)Math.Floor(tmp * fudge);
                                        if (itmp >= NumberofBins)
                                            itmp = NumberofBins - 1;
                                        DistanceHistogramBinCounts.addapoint(ThreadNo, itmp);
                                    }
                                }

                            }   // End Loop over GlobalPointIndex2
                        }   // End Loop over GlobalPointIndex1
                    });  // End Parallel Section defining GlobalPointIndex1

                    DistanceHistogramBinCounts.sumoverthreadsandmpi();
                    double NumberinHistogram = DistanceHistogramBinCounts.TotalNumberofPoints;

                    // Find Bucket Distance Cuts from histogram
                    for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                    {
                        double target = (PWCUtility.BucketFractions[BucketIndex] * NumberinHistogram);
                        double TotalbinCounts = 0;
                        for (int BinIndex = 0; BinIndex < NumberofBins; BinIndex++)
                        {
                            if (TotalbinCounts >= target)
                            {
                                BucketRadii[BucketIndex] = (BinIndex + 0.5) / fudge;
                                break;
                            }
                            TotalbinCounts += DistanceHistogramBinCounts.TotalSum[BinIndex];
                        }
                    }

                    GlobalReductions.FindManyMinValuewithIndex[] FindCentersbybuckets = new GlobalReductions.FindManyMinValuewithIndex[PWCUtility.NumberofBuckets];
                    for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                    {
                        FindCentersbybuckets[BucketIndex] = new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
                    }

                    //  Loop over points
                    Parallel.For(0, Program.ParallelOptions.MaxDegreeOfParallelism, Program.ParallelOptions, (ThreadNo) =>
                    {
                        int indexlen = PWCUtility.PointsperThread[ThreadNo];
                        int beginpoint = PWCUtility.StartPointperThread[ThreadNo] - PWCUtility.PointStart_Process;

                        for (int index = beginpoint; index < indexlen + beginpoint; index++)
                        {
                            int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;
                            int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                            if (group != group1)
                                continue;
                            double[] BucketCounts = new double[PWCUtility.NumberofBuckets];
                            for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                                BucketCounts[BucketIndex] = 0.0;

                            for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global; GlobalPointIndex2++)
                            {
                                if (GlobalPointIndex1 == GlobalPointIndex2)
                                    continue;
                                int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                                if (group != group2)
                                    continue;
                                double tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (tmp > PWCUtility.MinimumDistanceCut)
                                {
                                    if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2))
                                    {
                                        for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                                        {
                                            if (tmp <= BucketRadii[BucketIndex])
                                                BucketCounts[BucketIndex] += 1.0;
                                        }
                                    }
                                }

                            }   // End Loop over GlobalPointIndex2

                            for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                            {
                                if( BucketCounts[BucketIndex] < 0.5)
                                    continue;
                                double tmp = (double)GroupCount[group] - BucketCounts[BucketIndex];
                                FindCentersbybuckets[BucketIndex].addapoint(ThreadNo, GlobalPointIndex1, tmp);
                            }
                        }   // End Loop over GlobalPointIndex1
                    });  // End Parallel Section defining GlobalPointIndex1

                    for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                    {
                        FindCentersbybuckets[BucketIndex].sumoverthreadsandmpi();
                    }

                    // Output Bucket results
                    for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++)
                    {
                        PWCUtility.SALSAPrint(1, "\n----------------------\nGroup " + (group + StartPosition).ToString() + " Max "
                            + GlobalGroupMax[group].ToString("F4") + " Count " + GroupCount[group].ToString() + " Find Buckets CutOff Fraction "
                            + PWCUtility.BucketFractions[BucketIndex].ToString("F3") + " Avg Bucket Distce " + BucketRadii[BucketIndex].ToString("F4"));

                        TopBucketDistance[0] = 0.0;
                        int NumBucketEntriesFound = 0;
                        for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++)
                        {
                            int GlobalPointIndex1 = FindCentersbybuckets[BucketIndex].OrderedIndexValue[CenterIndex1];
                            if (GlobalPointIndex1 < 0)
                            {
                                PWCUtility.SALSAPrint(0, "Group " + group.ToString() + " Bucket " + BucketIndex.ToString() + " Center Index " + CenterIndex1.ToString() + " Undefined");
                                continue;
                            }
                            ++NumBucketEntriesFound;

                            // Add into Global Rating List
                            int useposition = countindicesfound;
                            if (countindicesfound > 0)
                            {
                                for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++)
                                {
                                    if (CombinedListIndices[loopindexlist] != GlobalPointIndex1)
                                        continue;
                                    useposition = loopindexlist;
                                    break;
                                }
                            }
                            if (useposition < countindicesfound)
                            {
                                CombinedListRatingsCount[useposition] += 1.0/PWCUtility.NumberofBuckets;
                                CombinedListRatings[useposition] += CenterIndex1/PWCUtility.NumberofBuckets;
                                CombinedListSourceSeq[useposition] = true;
                            }
                            else
                            {
                                CombinedListIndices[useposition] = GlobalPointIndex1;
                                CombinedListRatingsCount[useposition] = 1.0/PWCUtility.NumberofBuckets;
                                CombinedListRatings[useposition] = CenterIndex1/PWCUtility.NumberofBuckets;
                                CombinedListSourceSeq[useposition] = true;
                                CombinedListSourceMDS[useposition] = false;
                                ++countindicesfound;
                            }
                            double Avg_Internal = 0.0;
                            double Avg_MeanMethod = 0.0;
                            int Num_MeanMethod = 0;
                            int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                            if (owner == PWCUtility.MPI_Rank)
                            {
                                for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++)
                                {
                                    int GlobalPointIndex2 = FindCentersbybuckets[BucketIndex].OrderedIndexValue[CenterIndex2];
                                    if (GlobalPointIndex2 < 0)
                                    {
                                        Exception e = PWCUtility.SALSAError("Means Group-2 " + group.ToString() + " Center " + CenterIndex2.ToString() + " Undefined");
                                        throw (e);
                                    }
                                    double tmp;
                                    if (CenterIndex1 != CenterIndex2)
                                    {
                                        tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                        if (CenterIndex1 == 0)
                                            TopBucketDistance[CenterIndex2] = tmp;
                                        Avg_Internal = Avg_Internal + tmp;
                                    }

                                    GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                                    if (GlobalPointIndex2 < 0)
                                    {
                                        if (CenterIndex2 == 0)
                                            TopMeanDistance[CenterIndex1] = 0.0;
                                    }
                                    else
                                    {
                                        tmp = PWCParallelism.getDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                        Avg_MeanMethod = Avg_MeanMethod + tmp;
                                        ++Num_MeanMethod;
                                        if (CenterIndex2 == 0)
                                            TopMeanDistance[CenterIndex1] = tmp;
                                    }
                                }
                            }
                            if (PWCUtility.MPI_Size > 1)
                            {
                                PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                                PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_Internal, owner);
                                PWCUtility.MPI_communicator.Broadcast<double>(ref Avg_MeanMethod, owner);
                                PWCUtility.MPI_communicator.Broadcast<int>(ref Num_MeanMethod, owner);
                                if (CenterIndex1 == 0)
                                    PWCUtility.MPI_communicator.Broadcast<double>(ref TopBucketDistance, owner);
                                PWCUtility.MPI_communicator.Broadcast<double>(ref TopMeanDistance[CenterIndex1], owner);
                                PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                            }
                            Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);
                            if(Num_MeanMethod > 0)
                                Avg_MeanMethod = Avg_MeanMethod / Num_MeanMethod;
                            if (PWCUtility.MPI_Rank == 0)
                            {
                                string Seqlabel = "";
                                if (LabelsAvailable)
                                {
                                    Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" + SequenceLengths[GlobalPointIndex1].ToString();
                                }
                                int NumPointsinBucket = (int)Math.Floor((double)GroupCount[group] - FindCentersbybuckets[BucketIndex].OrderedMinValue[CenterIndex1] + 0.5);
                                PointProperties[CenterIndex1] = "Measure=" + NumPointsinBucket.ToString() + " Method=LargestCountsBucket-"
                                    + BucketIndex.ToString() + " Group=" + (group + StartPosition).ToString() + Seqlabel
                                    + " Date=\"" + DateTime.Now.ToLocalTime() + "\"";
                                PWCUtility.SALSAPrint(1, "Count " + NumPointsinBucket.ToString("F2")
                                    + " Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-Internal " + Avg_Internal.ToString("E3")
                                    + " Distce-means " + Avg_MeanMethod.ToString("E3") + " Distce to top mean "
                                    + TopMeanDistance[CenterIndex1].ToString("E3") + " To Top Bucket " + TopBucketDistance[CenterIndex1].ToString("E3") + Seqlabel);
                            }
                        }

                        // Output
                        if (PWCUtility.MPI_Rank == 0)
                        {
                            if (NumBucketEntriesFound > 0)
                            {
                                for (int CenterIndex = 0; CenterIndex < NumBucketEntriesFound; CenterIndex++)
                                {
                                    UtilityIndex[CenterIndex] = FindCentersbybuckets[BucketIndex].OrderedIndexValue[CenterIndex];
                                }
                                WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, NumBucketEntriesFound, true);
                            }
                        }

                    }   // End BucketIndex Loop
                }

                //  Global Rating
                if (PWCUtility.MPI_Rank == 0)
                {
                    PWCUtility.SALSAPrint(1, "\n---------------- Best Rated Points");
                    double bestrating = -1.0;
                    int bestposition = -1;
                    for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++)
                        CombinedListRatings[loopindexlist] = CombinedListRatings[loopindexlist] / CombinedListRatingsCount[loopindexlist];

                    int bestratedpoints = 0;
                    for (int loopindexlist1 = 0; loopindexlist1 < countindicesfound; loopindexlist1++)
                    {
                        bestrating = -1.0;
                        bestposition = -1;
                        for (int loopindexlist2 = 0; loopindexlist2 < countindicesfound; loopindexlist2++)
                        {
                            if (CombinedListRatingsCount[loopindexlist2] < -0.5)
                                continue;
                            if ((bestposition >= 0) && (CombinedListRatings[loopindexlist2] >= bestrating))
                                continue;
                            bestposition = loopindexlist2;
                            bestrating = CombinedListRatings[loopindexlist2];
                        }
                        if (CombinedListRatingsCount[bestposition] > 0.6)
                        {
                            int GlobalPointIndex = CombinedListIndices[bestposition];
                            if (GlobalPointIndex < 0)
                                continue;

                            string Seqlabel = "";
                            if (LabelsAvailable)
                            {
                                Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex] + " Length=" + SequenceLengths[GlobalPointIndex].ToString();
                            }
                            string source = "";
                            if (CombinedListSourceMDS[bestposition])
                            {
                                if (CombinedListSourceSeq[bestposition])
                                    source = " Both";
                                else
                                    source = " MDS ";
                            }
                            else
                                source = " Seq ";

                            UtilityIndex[bestratedpoints] = GlobalPointIndex;
                            PointProperties[bestratedpoints] = "Measure=" +
                                                               CombinedListRatings[bestposition].ToString("F2") +
                                                               " Count=" +
                                                               CombinedListRatingsCount[bestposition].ToString("F1")
                                                               + " Source=" + source.Trim() + " Method=OverallBest" + " Group=" + (group + StartPosition).ToString() + Seqlabel + " Date=\"" + DateTime.Now.ToLocalTime() + "\"";

                            ++bestratedpoints;
                            PWCUtility.SALSAPrint(1, "Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex, 6) + " Rating " + CombinedListRatings[bestposition].ToString("F2")
                                 + " Count " + CombinedListRatingsCount[bestposition].ToString("F1") + source + Seqlabel);
                        }
                        CombinedListRatingsCount[bestposition] = -1.0;
                    }
                    if (bestratedpoints > 0)
                        WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, bestratedpoints, true);
                }   // End Output in rank 0 only of best rated points

            }   // End loop over groups
            return;
        }

        public static double getMDSDistanceValue(int GlobalPointIndex1, int GlobalPointIndex2)
        {
            if (GlobalPointIndex1 == GlobalPointIndex2)
                return 0.0;
            double tmp0 = MDSvalues[GlobalPointIndex1][0] - MDSvalues[GlobalPointIndex2][0];
            double tmp1 = MDSvalues[GlobalPointIndex1][1] - MDSvalues[GlobalPointIndex2][1];
            double tmp2 = MDSvalues[GlobalPointIndex1][2] - MDSvalues[GlobalPointIndex2][2];
            return Math.Sqrt(tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2);
        }

        public static double getMDSDistancefromPoint(int GlobalPointIndex, double[] MDSPosition)
        {
            double tmp0 = MDSvalues[GlobalPointIndex][0] - MDSPosition[0];
            double tmp1 = MDSvalues[GlobalPointIndex][1] - MDSPosition[1];
            double tmp2 = MDSvalues[GlobalPointIndex][2] - MDSPosition[2];
            return Math.Sqrt(tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2);
        }

        // Read MDS values
        public static void ReadMDS(string fname, double[][] MDSvaluestoberead, int BeginPoint, int PointstoRead)
        {
            char[] _sep = new[] { ' ', '\t' };

            bool success = false; 
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
                    int count = 0;
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            int index = int.Parse(splits[0]);
                            if (index < BeginPoint)
                                continue;
                            if (index >= BeginPoint + PointstoRead)
                                break;
                            MDSvaluestoberead[index - BeginPoint][0] = double.Parse(splits[1]);
                            MDSvaluestoberead[index - BeginPoint][1] = double.Parse(splits[2]);
                            MDSvaluestoberead[index - BeginPoint][2] = double.Parse(splits[3]);
                            count++;
                        }
                    }
                    if (count != PointstoRead) 
                    {
                        Exception e = PWCUtility.SALSAError("Illegal count on MDS file " + count.ToString() + " " + PointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading MDS data" + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = PWCUtility.SALSAError("MDS File read error " + fname);
                throw (e);
            }

        }   // End ReadMDS

        // Read Point Labels
        public static void ReadLabels(string fname, string[] SequenceLabels, int[] SequenceLengths, int BeginPoint, int PointstoRead)
        {
            char[] _sep = new[] { ' ', '\t' };

            bool success = false;
            string line ="NotSet";
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
                    while (!sr.EndOfStream)
                    {
                        line = sr.ReadLine();
                        if (!string.IsNullOrEmpty(line))
                        {
                            string[] splits = line.Trim().Split(_sep);
                            int index = int.Parse(splits[0]);
                            if (index < BeginPoint)
                                continue;
                            if (index >= BeginPoint + PointstoRead)
                                break;
                            if (splits.Length == 3)
                            {
                                SequenceLabels[index - BeginPoint] = splits[1];
                                SequenceLengths[index - BeginPoint] = int.Parse(splits[2]);
                            }
                            else if (splits.Length > 3)
                            {
                                // Assume splits[1] to splits[splits.Length - 2] contribute to the name
                                SequenceLabels[index - BeginPoint] =  string.Join(" ", splits, 1, (splits.Length - 2));
                                SequenceLengths[index - BeginPoint] = int.Parse(splits[splits.Length - 1]);
                            }
                            else
                            {
                                // Assume at least two splits where splits[0] is the index and splits[1] is the length
                                SequenceLabels[index - BeginPoint] = string.Empty;
                                SequenceLengths[index - BeginPoint] = int.Parse(splits[1]);
                            }
                            count++;
                        }
                    }
                    if (count != PointstoRead)
                    {
                        Exception e = PWCUtility.SALSAError("Illegal count on Label file " + count.ToString() + " " + PointstoRead.ToString());
                        throw (e);
                    }
                    success = true;
                }
                sr.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed reading Label data: count " + count.ToString() + " Line " + line + "\n" + e);
                throw (e);
            }
            if (!success)
            {
                Exception e = PWCUtility.SALSAError("Label File read error: count " + count.ToString() + " Line " + line + "\n" + fname);
                throw (e);
            }

        }   // End ReadSequenceLabels
 
        // Write General Point results into a file
        public static void WritePointProperties(string fname, int[] PointNumbers, string[] labels, int dataPoints, bool append)
        {
            if (PWCUtility.MPI_Rank != 0)
                return;
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
                        int index = PointNumbers[i];
                        string outputpoint = PWCUtility.LeftPrintFixedInteger(PointNumbers[i], 6);
                        sw.WriteLine("PointNumber=" + outputpoint + " " + labels[i]);
                    }
                }

                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing data" + e);
                throw (e);
            }
        }   // End WritePointProperties



    }   // End class FindCenters
}
