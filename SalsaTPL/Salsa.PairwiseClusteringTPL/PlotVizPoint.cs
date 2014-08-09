using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace Salsa.PairwiseClusteringTPL
{
    public class PlotVizPoint
    {
        private double _x;
        private double _y;
        private double _z;
        private int _index;
        private int _cluster;
        private string _label;

        public PlotVizPoint(double x, double y, double z, int index, int cluster)
        {
            _x = x;
            _y = y;
            _z = z;
            _index = index;
            _cluster = cluster;
        }

        public PlotVizPoint(double x, double y, double z, int index, int cluster, string label)
        {
            _x = x;
            _y = y;
            _z = z;
            _index = index;
            _cluster = cluster;
            _label = label;
        }

        public double X
        {
            get { return _x; }
            set { _x = value; }
        }

        public double Y
        {
            get { return _y; }
            set { _y = value; }
        }

        public double Z
        {
            get { return _z; }
            set { _z = value; }
        }

        public int Index
        {
            get { return _index; }
            set { _index = value; }
        }

        public int Cluster
        {
            get { return _cluster; }
            set { _cluster = value; }
        }

        public string Label
        {
            get { return _label; }
            set { _label = value; }
        }

        public XElement ToPvizPointElement()
        {
           return new XElement("point",
                             new XElement("key", _index),
                             new XElement("clusterkey", _cluster),
                             new XElement("label", _label),
                             new XElement("location",
                                          new XAttribute("x", _x),
                                          new XAttribute("y", _y),
                                          new XAttribute("z", _z)));
        }

        public PlotVizPoint Clone()
        {
            return new PlotVizPoint(_x,_y,_z,_index,_cluster,_label);
        }
    }
}
