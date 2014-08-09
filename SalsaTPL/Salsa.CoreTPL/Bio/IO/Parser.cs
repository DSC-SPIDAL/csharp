using System;
using System.Collections.Generic;
using System.IO;

namespace Salsa.Core.Bio.IO.Internal
{
    internal class Parser<TRecord> where TRecord : IFastaRecord
    {
        public const int _EOF = 0;
        public const int _string = 1;
        public const int maxT = 5;

        private const bool T = true;
        private const bool x = false;
        private const int minErrDist = 2;

        private static readonly bool[,] set =
            {
                {T, x, x, x, x, x, x}
            };

        private readonly IFastaRecordFactory<TRecord> _factory;
        private readonly List<TRecord> _records = new List<TRecord>();
        private int errDist = minErrDist;

        internal Errors errors;

        internal Token la; // lookahead token
        internal Scanner scanner;
        internal Token t; // last recognized token


        internal Parser(Scanner scanner, IFastaRecordFactory<TRecord> factory)
        {
            this.scanner = scanner;
            _factory = factory;
            errors = new Errors();
        }

        internal List<TRecord> Records
        {
            get { return _records; }
        }

        private void SynErr(int n)
        {
            if (errDist >= minErrDist) errors.SynErr(la.line, la.col, n);
            errDist = 0;
        }

        internal void SemErr(string msg)
        {
            if (errDist >= minErrDist) errors.SemErr(t.line, t.col, msg);
            errDist = 0;
        }

        private void Get()
        {
            for (;;)
            {
                t = la;
                la = scanner.Scan();
                if (la.kind <= maxT)
                {
                    ++errDist;
                    break;
                }

                la = t;
            }
        }

        private void Expect(int n)
        {
            if (la.kind == n) Get();
            else
            {
                SynErr(n);
            }
        }

        private bool StartOf(int s)
        {
            return set[s, la.kind];
        }

        private void ExpectWeak(int n, int follow)
        {
            if (la.kind == n) Get();
            else
            {
                SynErr(n);
                while (!StartOf(follow)) Get();
            }
        }


        private bool WeakSeparator(int n, int syFol, int repFol)
        {
            int kind = la.kind;
            if (kind == n)
            {
                Get();
                return true;
            }
            else if (StartOf(repFol))
            {
                return false;
            }
            else
            {
                SynErr(n);
                while (!(set[syFol, kind] || set[repFol, kind] || set[0, kind]))
                {
                    Get();
                    kind = la.kind;
                }
                return StartOf(syFol);
            }
        }


        private void FastaFile()
        {
            TRecord record = default(TRecord);
            while (la.kind == 2)
            {
                FastaEntry(ref record);
                _records.Add(record);
            }
        }

        private void FastaEntry(ref TRecord record)
        {
            record = _factory.Create();
            string label = string.Empty;
            string residues = string.Empty;
            LabelLine(out label);
            record.Label = label;
            while (la.kind == 1)
            {
                ResidueLine(out residues);
                record.Residues += residues.ToUpper();
            }
        }

        private void LabelLine(out string label)
        {
            Expect(2);
            Expect(1);
            label = t.val.Trim();
            while (la.kind == 3 || la.kind == 4)
            {
                if (la.kind == 3)
                {
                    Get();
                }
                else
                {
                    Get();
                }
            }
        }

        private void ResidueLine(out string residues)
        {
            Expect(1);
            residues = t.val.Trim();
            while (la.kind == 3 || la.kind == 4)
            {
                if (la.kind == 3)
                {
                    Get();
                }
                else
                {
                    Get();
                }
            }
        }


        internal void Parse()
        {
            la = new Token();
            la.val = "";
            Get();
            FastaFile();

            Expect(0);
        }
    }

    // end Parser

    internal class Errors
    {
        internal int count = 0; // number of errors detected
        internal string errMsgFormat = "-- line {0} col {1}: {2}"; // 0=line, 1=column, 2=text
        internal TextWriter errorStream = Console.Out; // error messages go to this stream

        internal void SynErr(int line, int col, int n)
        {
            string s;
            switch (n)
            {
                case 0:
                    s = "EOF expected";
                    break;
                case 1:
                    s = "string expected";
                    break;
                case 2:
                    s = "\">\" expected";
                    break;
                case 3:
                    s = "\"\\r\" expected";
                    break;
                case 4:
                    s = "\"\\n\" expected";
                    break;
                case 5:
                    s = "??? expected";
                    break;

                default:
                    s = "error " + n;
                    break;
            }
            errorStream.WriteLine(errMsgFormat, line, col, s);
            count++;
        }

        internal void SemErr(int line, int col, string s)
        {
            errorStream.WriteLine(errMsgFormat, line, col, s);
            count++;
        }

        internal void SemErr(string s)
        {
            errorStream.WriteLine(s);
            count++;
        }

        internal void Warning(int line, int col, string s)
        {
            errorStream.WriteLine(errMsgFormat, line, col, s);
        }

        internal void Warning(string s)
        {
            errorStream.WriteLine(s);
        }
    }

    // Errors

    internal class FatalError : Exception
    {
        internal FatalError(string m) : base(m)
        {
        }
    }
}