using System;
using System.IO;

namespace Salsa.Core.Blas
{
    public class PartialMatrixTextReader
    {
        private readonly int _colCount;
        private readonly char _delimiter;
        private readonly int _globalColStartIndex;
        private readonly int _globalRowStartIndex;
        private readonly int _rowCount;

        #region Constructors

        public PartialMatrixTextReader(Range globalRowRange, Range globalColumnRange)
            : this(
                globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                globalColumnRange.EndIndex, ',')
        {
        }

        public PartialMatrixTextReader(Range globalRowRange, int globalColumnStartIndex, int globalColumnEndIndex)
            : this(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnStartIndex, globalColumnEndIndex, ','
                )
        {
        }

        public PartialMatrixTextReader(int globalRowStartIndex, int globalRowEndIndex, Range globalColumnRange)
            : this(globalRowStartIndex, globalRowEndIndex, globalColumnRange.StartIndex, globalColumnRange.EndIndex, ','
                )
        {
        }

        public PartialMatrixTextReader(Range globalRowRange, Range globalColumnRange, char delimiter)
            : this(
                globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                globalColumnRange.EndIndex, delimiter)
        {
        }

        public PartialMatrixTextReader(Range globalRowRange, int globalColumnStartIndex, int globalColumnEndIndex,
                                       char delimiter)
            : this(
                globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnStartIndex, globalColumnEndIndex,
                delimiter)
        {
        }

        public PartialMatrixTextReader(int globalRowStartIndex, int globalRowEndIndex, Range globalColumnRange,
                                       char delimiter)
            : this(
                globalRowStartIndex, globalRowEndIndex, globalColumnRange.StartIndex, globalColumnRange.EndIndex,
                delimiter)
        {
        }

        public PartialMatrixTextReader(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                                       int globalColumnEndIndex, char delimiter)
        {
            _globalRowStartIndex = globalRowStartIndex;
            _globalColStartIndex = globalColumnStartIndex;
            _rowCount = globalRowEndIndex - globalRowStartIndex + 1;
            _colCount = globalColumnEndIndex - globalColumnStartIndex + 1;
            _delimiter = delimiter;
        }

        #endregion

        public PartialMatrix<short> ReadInt16(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadInt16(reader);
            }
        }

        public PartialMatrix<ushort> ReadUInt16(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadUInt16(reader);
            }
        }

        public PartialMatrix<int> ReadInt32(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadInt32(reader);
            }
        }

        public PartialMatrix<uint> ReadUInt32(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadUInt32(reader);
            }
        }

        public PartialMatrix<float> ReadSingle(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadSingle(reader);
            }
        }

        public PartialMatrix<double> ReadDouble(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadDouble(reader);
            }
        }

        public PartialMatrix<short> ReadInt16(StreamReader reader)
        {
            short[][] elements = PartialMatrix<short>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToInt16(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<short>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                            _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<ushort> ReadUInt16(StreamReader reader)
        {
            ushort[][] elements = PartialMatrix<ushort>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToUInt16(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<ushort>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                             _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<int> ReadInt32(StreamReader reader)
        {
            int[][] elements = PartialMatrix<int>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToInt32(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<int>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                          _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<uint> ReadUInt32(StreamReader reader)
        {
            uint[][] elements = PartialMatrix<uint>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToUInt32(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<uint>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                           _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<float> ReadSingle(StreamReader reader)
        {
            float[][] elements = PartialMatrix<float>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToSingle(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<float>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                            _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<double> ReadDouble(StreamReader reader)
        {
            double[][] elements = PartialMatrix<double>.CreateElements(_rowCount, _colCount);

            int row = 0;

            // Read the Body
            while (reader.EndOfStream == false)
            {
                string[] fields = reader.ReadLine().Split(new[] {_delimiter});

                if (fields.Length != _colCount)
                {
                    throw new InvalidDataException("Unexpected number of columns on line " + row);
                }

                for (int col = 0; col < fields.Length; col++)
                {
                    elements[row][col] = Convert.ToDouble(fields[col].Trim());
                }

                row++;
            }

            if (row != _rowCount)
            {
                throw new InvalidDataException("Unexpected number of rows in file");
            }

            return new PartialMatrix<double>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                             _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }
    }
}