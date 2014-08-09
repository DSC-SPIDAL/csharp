using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Salsa.PairwiseClusteringTPL
{
    public class CenterInfo
    {
        private int _pnum;
        private double _measure;
        private string _method;
        private int _cluster;
        private string _seqName;
        private int _seqLength;
        private string _source;
        private double _count;

        private PlotVizPoint _plotVizPoint;
        public CenterInfo(int pnum, double measure, string method, int cluster, string seqName, int seqLength)
        {
            _pnum = pnum;
            _measure = measure;
            _method = method;
            _cluster = cluster;
            _seqName = seqName;
            _seqLength = seqLength;
        }

        public CenterInfo(int pnum, double measure, string method, int cluster, string seqName, int seqLength, string source, double count)
        {
            _pnum = pnum;
            _measure = measure;
            _method = method;
            _cluster = cluster;
            _seqName = seqName;
            _seqLength = seqLength;
            _source = source;
            _count = count;
        }

        public int Pnum
        {
            get { return _pnum; }
            set { _pnum = value; }
        }

        public double Measure
        {
            get { return _measure; }
            set { _measure = value; }
        }

        public string Method
        {
            get { return _method; }
            set { _method = value; }
        }

        public int Cluster
        {
            get { return _cluster; }
            set { _cluster = value; }
        }

        public string SeqName
        {
            get { return _seqName; }
            set { _seqName = value; }
        }

        public int SeqLength
        {
            get { return _seqLength; }
            set { _seqLength = value; }
        }

        public string Source
        {
            get { return _source; }
            set { _source = value; }
        }

        public double Count
        {
            get { return _count; }
            set { _count = value; }
        }

        public PlotVizPoint PlotVizPoint
        {
            get { return _plotVizPoint; }
            set { _plotVizPoint = value; }
        }

        public new string ToString()
        {
            return "PointNumber=" + _pnum + "\tMeasure=" + _measure + "\tMethod=" + _method + "\tGroup=" + _cluster +
                   "\tSequence=" + _seqName + "\tLength=" + _seqLength;
        }
    }
}
