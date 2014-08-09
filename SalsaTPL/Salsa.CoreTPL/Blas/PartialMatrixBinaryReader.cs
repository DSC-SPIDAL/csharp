using System.IO;

namespace Salsa.Core.Blas
{
    public class PartialMatrixBinaryReader
    {
        private readonly int _colCount;
        private readonly int _globalColStartIndex;
        private readonly int _globalRowStartIndex;
        private readonly int _rowCount;

        #region Constructors

        public PartialMatrixBinaryReader(Range globalRowRange, Range globalColumnRange)
            : this(
                globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                globalColumnRange.EndIndex)
        {
        }

        public PartialMatrixBinaryReader(Range globalRowRange, int globalColumnStartIndex, int globalColumnEndIndex)
            : this(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnStartIndex, globalColumnEndIndex)
        {
        }

        public PartialMatrixBinaryReader(int globalRowStartIndex, int globalRowEndIndex, Range globalColumnRange)
            : this(globalRowStartIndex, globalRowEndIndex, globalColumnRange.StartIndex, globalColumnRange.EndIndex)
        {
        }

        public PartialMatrixBinaryReader(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                                         int globalColumnEndIndex)
        {
            _globalRowStartIndex = globalRowStartIndex;
            _globalColStartIndex = globalColumnStartIndex;
            _rowCount = globalRowEndIndex - globalRowStartIndex + 1;
            _colCount = globalColumnEndIndex - globalColumnStartIndex + 1;
        }

        #endregion

        public PartialMatrix<short> ReadInt16(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt16(reader);
                }
            }
        }

        public PartialMatrix<ushort> ReadUInt16(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt16(reader);
                }
            }
        }

        public PartialMatrix<int> ReadInt32(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt32(reader);
                }
            }
        }

        public PartialMatrix<uint> ReadUInt32(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt32(reader);
                }
            }
        }

        public PartialMatrix<float> ReadSingle(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadSingle(reader);
                }
            }
        }

        public PartialMatrix<double> ReadDouble(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadDouble(reader);
                }
            }
        }

        public PartialMatrix<short> ReadInt16(BinaryReader reader)
        {
            short[][] elements = PartialMatrix<short>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt16();
                }
            }

            return new PartialMatrix<short>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                            _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<ushort> ReadUInt16(BinaryReader reader)
        {
            ushort[][] elements = PartialMatrix<ushort>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt16();
                }
            }

            return new PartialMatrix<ushort>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                             _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<int> ReadInt32(BinaryReader reader)
        {
            int[][] elements = PartialMatrix<int>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt32();
                }
            }

            return new PartialMatrix<int>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                          _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<uint> ReadUInt32(BinaryReader reader)
        {
            uint[][] elements = PartialMatrix<uint>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt32();
                }
            }

            return new PartialMatrix<uint>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                           _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<float> ReadSingle(BinaryReader reader)
        {
            float[][] elements = PartialMatrix<float>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadSingle();
                }
            }

            return new PartialMatrix<float>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                            _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }

        public PartialMatrix<double> ReadDouble(BinaryReader reader)
        {
            double[][] elements = PartialMatrix<double>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadDouble();
                }
            }

            return new PartialMatrix<double>(_globalRowStartIndex, _globalRowStartIndex + _rowCount - 1,
                                             _globalColStartIndex, _globalColStartIndex + _colCount - 1, elements);
        }
    }
}