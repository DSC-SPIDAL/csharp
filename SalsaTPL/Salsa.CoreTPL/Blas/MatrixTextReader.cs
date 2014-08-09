using System;
using System.IO;

namespace Salsa.Core.Blas
{
    public class MatrixTextReader
    {
        private readonly int _colCount;
        private readonly char _delimiter;
        private readonly int _rowCount;

        public MatrixTextReader(int rowCount, int colCount)
            : this(rowCount, colCount, ',')
        {
        }

        public MatrixTextReader(int rowCount, int colCount, char delimiter)
        {
            _rowCount = rowCount;
            _colCount = colCount;
            _delimiter = delimiter;
        }

        public Matrix<short> ReadInt16(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadInt16(reader);
            }
        }

        public Matrix<ushort> ReadUInt16(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadUInt16(reader);
            }
        }

        public Matrix<int> ReadInt32(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadInt32(reader);
            }
        }

        public Matrix<uint> ReadUInt32(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadUInt32(reader);
            }
        }

        public Matrix<float> ReadSingle(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadSingle(reader);
            }
        }

        public Matrix<double> ReadDouble(string fileName)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadDouble(reader);
            }
        }

        public Matrix<short> ReadInt16(StreamReader reader)
        {
            short[][] elements = Matrix<short>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<short>(elements);
        }

        public Matrix<ushort> ReadUInt16(StreamReader reader)
        {
            ushort[][] elements = Matrix<ushort>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<ushort>(elements);
        }

        public Matrix<int> ReadInt32(StreamReader reader)
        {
            int[][] elements = Matrix<int>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<int>(elements);
        }

        public Matrix<uint> ReadUInt32(StreamReader reader)
        {
            uint[][] elements = Matrix<uint>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<uint>(elements);
        }

        public Matrix<float> ReadSingle(StreamReader reader)
        {
            float[][] elements = Matrix<float>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<float>(elements);
        }

        public Matrix<double> ReadDouble(StreamReader reader)
        {
            double[][] elements = Matrix<double>.CreateElements(_rowCount, _colCount);

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

            return new Matrix<double>(elements);
        }

        public Matrix<short> ReadInt16(string fileName, int rowStartIndex, int rowCount)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadInt16(reader, rowStartIndex, rowCount);
            }
        }

        public Matrix<ushort> ReadUInt16(string fileName, int rowStartIndex, int rowCount)
        {
            using (var reader = new StreamReader(fileName))
            {
                return ReadUInt16(reader, rowStartIndex, rowCount);
            }
        }

        public Matrix<short> ReadInt16(StreamReader reader, int rowStartIndex, int rowCount)
        {
            short[][] elements = Matrix<short>.CreateElements(_rowCount, _colCount);
            var separators = new[] {_delimiter};

            int row = 0;

            while (reader.EndOfStream == false)
            {
                string line = reader.ReadLine();

                if (row >= rowStartIndex + rowCount)
                {
                    break;
                }
                else if (row >= rowStartIndex)
                {
                    string[] fields = line.Split(separators);

                    if (fields.Length != _colCount)
                    {
                        throw new InvalidDataException("Unexpected number of columns on line " + row);
                    }

                    for (int col = 0; col < fields.Length; col++)
                    {
                        elements[row][col] = Convert.ToInt16(fields[col].Trim());
                    }
                }

                row++;
            }

            return new Matrix<short>(elements);
        }

        public Matrix<ushort> ReadUInt16(StreamReader reader, int rowStartIndex, int rowCount)
        {
            ushort[][] elements = Matrix<ushort>.CreateElements(_rowCount, _colCount);
            var separators = new[] {_delimiter};

            int row = 0;

            while (reader.EndOfStream == false)
            {
                string line = reader.ReadLine();

                if (row >= rowStartIndex + rowCount)
                {
                    break;
                }
                else if (row >= rowStartIndex)
                {
                    string[] fields = line.Split(separators);

                    if (fields.Length != _colCount)
                    {
                        throw new InvalidDataException("Unexpected number of columns on line " + row);
                    }

                    for (int col = 0; col < fields.Length; col++)
                    {
                        elements[row][col] = Convert.ToUInt16(fields[col].Trim());
                    }
                }

                row++;
            }

            return new Matrix<ushort>(elements);
        }
    }
}