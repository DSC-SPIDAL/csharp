using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;

namespace Salsa.Core.Blas
{
    public class MatrixBinaryReader
    {
        private readonly int _colCount;
        private readonly int _rowCount;

        public MatrixBinaryReader(int rowCount, int colCount)
        {
            _rowCount = rowCount;
            _colCount = colCount;
        }

        #region Read matrix through global row and column counts as set via constructor

        public Matrix<short> ReadInt16(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt16(reader);
                }
            }
        }

        public Matrix<ushort> ReadUInt16(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt16(reader);
                }
            }
        }

        public Matrix<int> ReadInt32(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt32(reader);
                }
            }
        }

        public Matrix<uint> ReadUInt32(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt32(reader);
                }
            }
        }

        public Matrix<float> ReadSingle(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadSingle(reader);
                }
            }
        }

        public Matrix<double> ReadDouble(string fileName)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadDouble(reader);
                }
            }
        }

        public Matrix<short> ReadInt16(BinaryReader reader)
        {
            short[][] elements = Matrix<short>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt16();
                }
            }

            return new Matrix<short>(elements);
        }

        public Matrix<ushort> ReadUInt16(BinaryReader reader)
        {
            ushort[][] elements = Matrix<ushort>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt16();
                }
            }

            return new Matrix<ushort>(elements);
        }

        public Matrix<int> ReadInt32(BinaryReader reader)
        {
            int[][] elements = Matrix<int>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt32();
                }
            }

            return new Matrix<int>(elements);
        }

        public Matrix<uint> ReadUInt32(BinaryReader reader)
        {
            uint[][] elements = Matrix<uint>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt32();
                }
            }

            return new Matrix<uint>(elements);
        }

        public Matrix<float> ReadSingle(BinaryReader reader)
        {
            float[][] elements = Matrix<float>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadSingle();
                }
            }

            return new Matrix<float>(elements);
        }

        public Matrix<double> ReadDouble(BinaryReader reader)
        {
            double[][] elements = Matrix<double>.CreateElements(_rowCount, _colCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadDouble();
                }
            }

            return new Matrix<double>(elements);
        }

        #endregion

        #region Read matrix from a set of files with consecutive numbering

        /// <summary>
        /// Reads a matrix of Int16 from the set of files in the given path with the given prefix and
        /// indices lying between startIndex and endIndex inclusively.
        /// </summary>
        /// <param name="path">The local path where partial matrix files are stored</param>
        /// <param name="prefix">The name prefix of the files</param>
        /// <param name="ext">File extension of the partial matrix files</param>
        /// <param name="startIndex">The start index of file names to read</param>
        /// <param name="endIndex">The end index (inclusively) of file names to read</param>
        /// <returns></returns>
        public Matrix<short> ReadInt16(string path, string prefix, string ext, int startIndex, int endIndex,
                                       Range[] ranges)
        {
            short[][] elements = Matrix<short>.CreateElements(_rowCount, _colCount);
            string fileName;
            int row = 0;
            for (int i = startIndex; i <= endIndex; i++)
            {
                fileName = string.Format("{0}_{1}{2}", prefix, i, ext);
                using (
                    Stream stream = File.Open(Path.Combine(path, fileName), FileMode.Open, FileAccess.Read,
                                              FileShare.Read))
                {
                    using (var reader = new BinaryReader(stream))
                    {
                        for (int j = ranges[i].StartIndex; j <= ranges[i].EndIndex; j++)
                        {
                            for (int k = 0; k < _colCount; k++)
                            {
                                elements[row][k] = reader.ReadInt16();
                            }
                            row++;
                        }
                    }
                }
            }
            return new Matrix<short>(elements);
        }

        #endregion

        #region Read matrix only for a given set of rows starting from a particular index

        public Matrix<short> ReadInt16(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt16(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<ushort> ReadUInt16(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt16(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<int> ReadInt32(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadInt32(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<uint> ReadUInt32(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadUInt32(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<float> ReadSingle(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadSingle(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<double> ReadDouble(string fileName, int startRowIndex, int rowCount)
        {
            using (Stream stream = File.Open(fileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (var reader = new BinaryReader(stream))
                {
                    return ReadDouble(reader, startRowIndex, rowCount);
                }
            }
        }

        public Matrix<short> ReadInt16(BinaryReader reader, int startRowIndex, int rowCount)
        {
            short[][] elements = Matrix<short>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (short);

            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt16();
                }
            }

            return new Matrix<short>(elements);
        }

        public Matrix<ushort> ReadUInt16(BinaryReader reader, int startRowIndex, int rowCount)
        {
            ushort[][] elements = Matrix<ushort>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (ushort);
            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt16();
                }
            }

            return new Matrix<ushort>(elements);
        }

        public Matrix<int> ReadInt32(BinaryReader reader, int startRowIndex, int rowCount)
        {
            int[][] elements = Matrix<int>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (int);
            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadInt32();
                }
            }

            return new Matrix<int>(elements);
        }

        public Matrix<uint> ReadUInt32(BinaryReader reader, int startRowIndex, int rowCount)
        {
            uint[][] elements = Matrix<uint>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (uint);
            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadUInt32();
                }
            }

            return new Matrix<uint>(elements);
        }

        public Matrix<float> ReadSingle(BinaryReader reader, int startRowIndex, int rowCount)
        {
            float[][] elements = Matrix<float>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (float);
            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadSingle();
                }
            }

            return new Matrix<float>(elements);
        }

        public Matrix<double> ReadDouble(BinaryReader reader, int startRowIndex, int rowCount)
        {
            double[][] elements = Matrix<double>.CreateElements(rowCount, _colCount);

            long offset = startRowIndex;
            offset *= _colCount;
            offset *= sizeof (double);
            reader.BaseStream.Seek(offset, SeekOrigin.Begin);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    elements[i][j] = reader.ReadDouble();
                }
            }

            return new Matrix<double>(elements);
        }

        #endregion
    }

    public class MatrixBinaryInt16RowReader : IDisposable, IEnumerable<MatrixRow<short>>
    {
        private readonly int _colCount;
        private readonly int _rowCount;
        private BinaryReader _reader;

        public MatrixBinaryInt16RowReader(string fileName, int rowCount, int colCount)
        {
            _rowCount = rowCount;
            _colCount = colCount;
            _reader = new BinaryReader(File.OpenRead(fileName));
        }

        #region IDisposable Members

        void IDisposable.Dispose()
        {
            Dispose(true);
        }

        #endregion

        #region IEnumerable<MatrixRow<short>> Members

        public IEnumerator<MatrixRow<short>> GetEnumerator()
        {
            for (int i = 0; i < _rowCount; i++)
            {
                var values = new short[_colCount];

                for (int j = 0; j < _colCount; j++)
                {
                    values[j] = _reader.ReadInt16();
                }

                yield return new MatrixRow<short>(i, values);
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        #endregion

        public virtual void Close()
        {
            Dispose(true);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposing)
            {
                if (_reader != null)
                {
                    _reader.Dispose();
                }
            }

            _reader = null;
        }
    }

    public class MatrixRow<T>
    {
        private readonly int _rowIndex;
        private readonly T[] _values;

        internal MatrixRow(int rowIndex, T[] values)
        {
            _rowIndex = rowIndex;
            _values = values;
        }

        public int RowIndex
        {
            get { return _rowIndex; }
        }

        public T this[int colIndex]
        {
            get { return _values[colIndex]; }
        }
    }
}