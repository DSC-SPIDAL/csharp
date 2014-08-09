using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Salsa.PairwiseClusteringTPL
{
    public class CenterFileReader : IDisposable
    {
        private readonly char[] _sep = new[] {' ', '\t'};
        private readonly char[] _eqsep = new[] {'='};
        private readonly StreamReader _reader;

        public CenterFileReader(string centerfile)
        {
            _reader = new StreamReader(centerfile);
        }

        public CenterInfo ReadCenterLine()
        {
            string line = _reader.ReadLine();
            
            string[] splits = line.Split(_sep, StringSplitOptions.RemoveEmptyEntries);
            int pnum = int.Parse(splits[0].Split(_eqsep)[1]);
            double measure = double.Parse(splits[1].Split(_eqsep)[1]);
            int methodIdx = 2;
            string source = string.Empty;
            double count = 0.0;
            if (splits[2].StartsWith("Count"))
            {
                methodIdx = 4;
                count = double.Parse(splits[2].Split(_eqsep)[1]);
                source = splits[3].Split(_eqsep)[1];
            }
            string method = splits[methodIdx].Split(_eqsep)[1];
            int group = int.Parse(splits[methodIdx + 1].Split(_eqsep)[1]);
            string seqName = splits[methodIdx + 2].Split(_eqsep)[1];
            for (int i = methodIdx + 3; i < splits.Length - 4; ++i)
            {
                seqName += (" " + splits[i]);
            }
            int seqLength = int.Parse(splits[splits.Length - 4].Split(_eqsep)[1]);
            return new CenterInfo(pnum, measure, method, group, seqName, seqLength, source, count);
        }

        public bool EndOfStream
        {
            get { return _reader.EndOfStream; }
        }


        public void Dispose()
        {
            _reader.Close();
            _reader.Dispose();
        }
    }
}
