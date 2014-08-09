using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace Salsa.PairwiseClusteringTPL
{
    public class PlotVizCluster
    {
        private int _clusterNumber;
        private string _label;

        public PlotVizCluster(int clusterNumber, string label)
        {
            _clusterNumber = clusterNumber;
            _label = label;
        }

        public PlotVizCluster(int clusterNumber)
        {
            _clusterNumber = clusterNumber;
        }

        public PlotVizCluster(string label)
        {
            _label = label;
        }

        public int ClusterNumber
        {
            get { return _clusterNumber; }
            set { _clusterNumber = value; }
        }

        public string Label
        {
            get { return _label; }
            set { _label = value; }
        }

        public XElement ToClusterElement(Color color, bool isDefault, double size)
        {
            return new XElement("cluster",
                             new XElement("key", _clusterNumber),
                             new XElement("label", _label),
                             new XElement("visible", 1),
                             new XElement("default", isDefault ? 1 : 0),
                             new XElement("color",
                                          new XAttribute("r", color.R),
                                          new XAttribute("g", color.G),
                                          new XAttribute("b", color.B),
                                          new XAttribute("a", color.A)),
                             new XElement("size", size));
        }
    }
}
