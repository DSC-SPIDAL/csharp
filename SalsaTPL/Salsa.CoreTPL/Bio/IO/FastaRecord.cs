using System;
using System.Collections.Generic;
using Salsa.Core.Bio.Algorithms;
using Salsa.Core.Bio.IO.Internal;

namespace Salsa.Core.Bio.IO
{
    public interface IFastaRecord
    {
        string Label { get; set; }
        string Residues { get; set; }
    }

    public interface IFastaRecordFactory<TRecord> where TRecord : IFastaRecord
    {
        TRecord Create();
    }

    public class FastaParser<TRecord> where TRecord : IFastaRecord
    {
        private readonly IFastaRecordFactory<TRecord> _recordFactory;

        public FastaParser(IFastaRecordFactory<TRecord> recordFactory)
        {
            _recordFactory = recordFactory;
        }

        public IList<TRecord> Parse(string fileName)
        {
            var scanner = new Scanner(fileName);
            var parser = new Parser<TRecord>(scanner, _recordFactory);
            parser.Parse();
            return parser.Records;
        }
    }


    public class FastaSimpleRecord : IFastaRecord
    {
        #region Members

        private string _residues = string.Empty;

        #endregion

        public FastaSimpleRecord()
        {
        }

        public FastaSimpleRecord(string label, string sequence)
        {
            _residues = sequence;
            Label = label;
        }

        #region IFastaRecord Members

        public string Label { get; set; }

        public string Residues
        {
            get { return _residues; }
            set { _residues = value; }
        }

        #endregion

        public override string ToString()
        {
            return Label + Environment.NewLine + Residues;
        }
    }

    public class FastaSimpleRecordFactory : IFastaRecordFactory<FastaSimpleRecord>
    {
        #region IFastaRecordFactory<FastaSimpleRecord> Members

        FastaSimpleRecord IFastaRecordFactory<FastaSimpleRecord>.Create()
        {
            return Create();
        }

        #endregion

        public FastaSimpleRecord Create()
        {
            return new FastaSimpleRecord();
        }
    }

    public class FastaSimpleParser : FastaParser<FastaSimpleRecord>
    {
        public FastaSimpleParser()
            : base(new FastaSimpleRecordFactory())
        {
        }
    }


    public class FastaSequenceFactory : IFastaRecordFactory<Sequence>
    {
        #region IFastaRecordFactory<Sequence> Members

        public Sequence Create()
        {
            return new Sequence();
        }

        #endregion
    }

    public class FastaSequenceParser : FastaParser<Sequence>
    {
        public FastaSequenceParser()
            : base(new FastaSequenceFactory())
        {
        }
    }
}