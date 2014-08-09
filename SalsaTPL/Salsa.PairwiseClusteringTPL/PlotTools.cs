using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml.Linq;
using SALSALibrary;

namespace Salsa.PairwiseClusteringTPL
{
    public class PlotTools
    {
        public static void CreatePlotWithCenters(string centerFile, string pointsFile, string clusterNumberFile, int numberOfCenterPointsToIncludeInEachCenterType, string centerPlotFile, string plotDescription)
        {
            /* Generate all types of center clusters per cluster 
             * 
             * Center clusters are,
             *  1. Original Min Mean
             *  2. MDS Min Mean
             *  3. MDS Center of Gravity (CoG)
             *  4. Overall Best
             *  5. Bucket Fraction 0
             *     Bucket Fraction 1 and so on
             *     
             * Number of center points to include in each center type = n 
             * n <= N, which is the number of center points found for each center type by PWC
             * N is specified through NumberOfCenters parameter in PWC
             * 
             * Assumes a center file from a PWC center finding run
             * Assumes a points file, which has each point mapped to its cluster in the format 
             *  PointNumber<TAB>Xcoord<TAB>Ycoord<TAB>Zcoord<TAB>ClusterNumber
             */

            
            /* Colors to use with PlotViz 
               reads color info from Matlab50.txt file */
            List<Color> matlab50Colors = GenerateMatlab50Colors();

            /* XML elements to hold points and clusters to be used in PlotViz file */
            XElement clustersElement = new XElement("clusters");
            XElement pointsElement = new XElement("points");

            /* Hashtable mapping point number to a PlotVizPoint data structure for the points in the given points file */
            Hashtable existingPointsTable = new Hashtable();

            /* Maximum number of points int the points file */
            int maxpnum;
            /* Maximum number of clusters that points are mapped to in the points file*/
            int maxcnum;

            ProcessPointsFile(pointsFile, clusterNumberFile, clustersElement, pointsElement, out maxpnum, out maxcnum, existingPointsTable, matlab50Colors);

            /* Table mapping each cluster (i.e. group) number to another table called method table
             * method table maps each method (e.g. smallest distance mean, smallest MDS distance mean, etc.) name to the list center points for that particular method
             * the order of points in the list is as same as in the given center file */
            Hashtable groupTable = ProcessCenterFile(centerFile);

            CreatePlotWithCentersInternal(centerPlotFile, plotDescription, clustersElement, pointsElement, maxpnum, existingPointsTable,
                                          maxcnum, matlab50Colors, groupTable,
                                          numberOfCenterPointsToIncludeInEachCenterType);

        }

        private static void CreatePlotWithCentersInternal(string centerPlotFile, string plotDescription, XElement clustersElement, XElement pointsElement, int maxpnum, Hashtable existingPointsTable, int maxcnum, List<Color> matlab50Colors, Hashtable groupTable, int numberOfCenterPointsToIncludeInEachCenterType)
        {
            ++maxcnum;
            foreach (DictionaryEntry groupToMethodTable in groupTable)
            {
                var group = (int)groupToMethodTable.Key; // group is the original cluster number
                var methodTable = (Hashtable)groupToMethodTable.Value;
                int methodCount = methodTable.Count;
                int tempCount = methodCount;
                foreach (DictionaryEntry methodToCenterPoints in methodTable)
                {
                    var method = (string)methodToCenterPoints.Key; // method is one of smallest distance mean, smallest MDS mean, etc.

                    // cluster number to be used in PlotViz for this center type 
                    int methodNumber = methodCount - tempCount--;
                    var clusterNumberForCenterType = group * methodCount + methodNumber + maxcnum;

                    // cluster name to be used in PlotViz for this center type
                    var centerTypeName = group + "." + method + ".centerpoints";

                    // add an XML element to represent this center type as a cluster in PlotViz
                    clustersElement.Add(CreateClusterElement(clusterNumberForCenterType, centerTypeName,
                                                             matlab50Colors[group % matlab50Colors.Count], false, 2.0,methodNumber));

                    var cps = (List<CenterInfo>)methodToCenterPoints.Value;
                    // Picking the topmost n point for each method
                    for (int i = 0; i < numberOfCenterPointsToIncludeInEachCenterType; i++)
                    {
                        CenterInfo cp = cps[i];
                        PlotVizPoint p = (PlotVizPoint)existingPointsTable[cp.Pnum];
                        pointsElement.Add(CreatePointElement(++maxpnum, clusterNumberForCenterType,
                                                             ("cluster:" + group + "-idx:" + p.Index + "method:" +
                                                              method),
                                                             p.X, p.Y, p.Z));
                    }
                }
            }

            XElement plotElement = CreatePlotElement(plotDescription, true);
            XElement plotvizElement = new XElement("plotviz");
            plotvizElement.Add(plotElement);
            plotvizElement.Add(clustersElement);
            plotvizElement.Add(pointsElement);
            plotvizElement.Save(centerPlotFile);

        }



        private static void ProcessPointsFile(string pointsFile, string clusterNumberFile, XElement clusters, XElement points, out int maxpnum, out int maxcnum, Hashtable pointsTable, List<Color> matlab50Colors)
        {
            using (StreamReader preader = new StreamReader(pointsFile), creader = new StreamReader(clusterNumberFile))
            {
                HashSet<int> clusterNumbers = new HashSet<int>();
                maxpnum = -1;
                while (!preader.EndOfStream)
                {
                    string pline = preader.ReadLine();
                    string cline = creader.ReadLine();
                    if (!string.IsNullOrEmpty(pline) && !string.IsNullOrEmpty(cline))
                    {
                        PlotVizPoint p = ReadPointLine(pline.Trim());
                        if (maxpnum < p.Index)
                        {
                            maxpnum = p.Index;
                        }
                        pointsTable.Add(p.Index, p);

                        int cnum = ReadCnum(cline);
                        p.Cluster = cnum;
                        if (!clusterNumbers.Contains(p.Cluster))
                        {
                            clusterNumbers.Add(p.Cluster);
                            clusters.Add(CreateClusterElement(p.Cluster,
                                                              p.Cluster.ToString(CultureInfo.InvariantCulture),
                                                              matlab50Colors[p.Cluster % matlab50Colors.Count], true, 0.1, Glyphs.Hexagon2D));
                        }
                        points.Add(CreatePointElement(p.Index, p.Cluster, string.Empty, p.X, p.Y, p.Z));
                    }
                }
                maxcnum = clusterNumbers.Max();
            }
        }

        private static int ReadCnum(string line)
        {
            char[] sep = new[] { ' ', '\t' };
            string[] splits = line.Split(sep, StringSplitOptions.RemoveEmptyEntries);
            return splits.Length == 2 ? int.Parse(splits[1]) : splits.Length == 5 ? int.Parse(splits[4]) : 0;
        }

        private static List<Color> GenerateMatlab50Colors()
        {
            using (Stream stream = Assembly.GetExecutingAssembly().GetManifestResourceStream("Salsa.PairwiseClusteringTPL.Matlab50.txt"))                                                                                              
            {
                if (stream != null)
                {
                    using (StreamReader reader = new StreamReader(stream))
                    {
                        List<Color> colors = new List<Color>();
                        char[] sep = new[] {' ', '\t'};
                        string[] splits;
                        string split;
                        int startIdx = 3;
                        int r, g, b, a;
                        while (!reader.EndOfStream)
                        {
                            string line = reader.ReadLine();
                            if (!string.IsNullOrEmpty(line))
                            {
                                splits = line.Trim().Split(sep);

                                split = splits[0];
                                r = int.Parse(split.Substring(startIdx, (split.Length - (startIdx + 1))));

                                split = splits[1];
                                g = int.Parse(split.Substring(startIdx, (split.Length - (startIdx + 1))));

                                split = splits[2];
                                b = int.Parse(split.Substring(startIdx, (split.Length - (startIdx + 1))));

                                split = splits[3];
                                a = int.Parse(split.Substring(startIdx, (split.Length - (startIdx + 1))));

                                colors.Add(Color.FromArgb(a, r, g, b));
                            }
                        }
                        return colors;
                    }
                }
                else
                {
                    throw new Exception("Unable to load embedded resource: Matlab50.txt");
                }
            }
        }

        private static PlotVizPoint ReadPointLine(string line)
        {
            char[] sep = new[] { ' ', '\t' };
            string[] splits = line.Split(sep, StringSplitOptions.RemoveEmptyEntries);
            PlotVizPoint p = new PlotVizPoint(double.Parse(splits[1]), double.Parse(splits[2]), double.Parse(splits[3]), int.Parse(splits[0]), int.Parse(splits[4]));
            return p;
        }

        private static CenterInfo ReadCenterLine(string line)
        {
            char[] sep = new[] { ' ', '\t' };
            char[] eqsep = new[] { '=' };
            string[] splits = line.Split(sep, StringSplitOptions.RemoveEmptyEntries);
            int pnum = int.Parse(splits[0].Split(eqsep)[1]);
            double measure = double.Parse(splits[1].Split(eqsep)[1]);
            int methodIdx = 2;
            string source = string.Empty;
            double count = 0.0;
            if (splits[2].StartsWith("Count"))
            {
                methodIdx = 4;
                count = double.Parse(splits[2].Split(eqsep)[1]);
                source = splits[3].Split(eqsep)[1];
            }
            string method = splits[methodIdx].Split(eqsep)[1];
            int group = int.Parse(splits[methodIdx + 1].Split(eqsep)[1]);
            string seqName = splits[methodIdx + 2].Split(eqsep)[1];
            for (int i = methodIdx + 3; i < splits.Length - 4; ++i)
            {
                seqName += (" " + splits[i]);
            }
            int seqLength = int.Parse(splits[splits.Length - 4].Split(eqsep)[1]);
            return new CenterInfo(pnum, measure, method, group, seqName, seqLength, source, count);
        }

        private static Hashtable ProcessCenterFile(string centerFile)
        {
            using (StreamReader reader = new StreamReader(centerFile))
            {
                Hashtable groupTable = new Hashtable();
                while (!reader.EndOfStream)
                {
                    CenterInfo cp = ReadCenterLine(reader.ReadLine());
                    AddToGroupTable(groupTable, cp);
                }
                return groupTable;
            }
        }

        private static void AddToGroupTable(Hashtable groupTable, CenterInfo cp)
        {
            if (groupTable.ContainsKey(cp.Cluster))
            {
                Hashtable methodTable = (Hashtable)groupTable[cp.Cluster];
                if (methodTable.ContainsKey(cp.Method))
                {
                    // Need a list to maintain the order of points
                    List<CenterInfo> cps = (List<CenterInfo>)methodTable[cp.Method];
                    cps.Add(cp);
                }
                else
                {
                    // Need a list to maintain the order of points
                    List<CenterInfo> cps = new List<CenterInfo> { cp };
                    methodTable[cp.Method] = cps;
                }
            }
            else
            {
                // Need a list to maintain the order of points
                List<CenterInfo> cps = new List<CenterInfo> { cp };
                Hashtable methodTable = new Hashtable();
                methodTable[cp.Method] = cps;
                groupTable[cp.Cluster] = methodTable;
            }
        }


        private static XElement CreatePlotElement(string name, bool glyphVisible)
        {
            XElement plot =
                new XElement("plot",
                             new XElement("title", name),
                             new XElement("pointsize", 1),
                             new XElement("glyph",
                                          new XElement("visible", glyphVisible ? 1 : 0),
                                          new XElement("scale", 1)),
                             new XElement("camera",
                                          new XElement("focumode", 0),
                                          new XElement("focus",
                                                       new XAttribute("x", 0),
                                                       new XAttribute("y", 0),
                                                       new XAttribute("z", 0))));
            return plot;
        }

        private static XElement CreateClusterElement(int key, string label, Color color, bool isDefault, double size, int shape)
        {
            XElement cluster =
                new XElement("cluster",
                             new XElement("key", key),
                             new XElement("label", label),
                             new XElement("visible", 1),
                             new XElement("default", isDefault ? 1 : 0),
                             new XElement("color",
                                          new XAttribute("r", color.R),
                                          new XAttribute("g", color.G),
                                          new XAttribute("b", color.B),
                                          new XAttribute("a", color.A)),
                             new XElement("size", size),
                             new XElement("shape", shape));
            return cluster;
        }

        private static XElement CreatePointElement(int key, int clusterKey, string label, double x, double y, double z)
        {
            XElement point =
                new XElement("point",
                             new XElement("key", key),
                             new XElement("clusterkey", clusterKey),
                             new XElement("label", label),
                             new XElement("location",
                                          new XAttribute("x", x),
                                          new XAttribute("y", y),
                                          new XAttribute("z", z)));
            return point;
        }

        struct Glyphs
        {
            public static int Triangle2D = 0;
            public static int Rectangle2D = 1;
            public static int Pentagon2D = 2;
            public static int Hexagon2D = 3;
            public static int Tetrahedron3D = 4;
            public static int Cube3D = 5;
            public static int Sphere3D = 6;
            public static int Cylinder3D = 7;
        }
    }
    
}
