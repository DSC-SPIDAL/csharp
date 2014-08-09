using System.IO;

namespace Salsa.Core.Blas
{
    public sealed unsafe class MatrixBinaryWriter
    {
        public void Write(Matrix<short> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<ushort> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<int> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<uint> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<float> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<double> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(Matrix<short> matrix, BinaryWriter writer)
        {
            short[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (short* dataPtr = element[i])
                {
                    short* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }

        public void Write(Matrix<ushort> matrix, BinaryWriter writer)
        {
            ushort[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (ushort* dataPtr = element[i])
                {
                    ushort* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }

        public void Write(Matrix<int> matrix, BinaryWriter writer)
        {
            int[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (int* dataPtr = element[i])
                {
                    int* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }

        public void Write(Matrix<uint> matrix, BinaryWriter writer)
        {
            uint[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (uint* dataPtr = element[i])
                {
                    uint* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }

        public void Write(Matrix<float> matrix, BinaryWriter writer)
        {
            float[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (float* dataPtr = element[i])
                {
                    float* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }

        public void Write(Matrix<double> matrix, BinaryWriter writer)
        {
            double[][] element = matrix.Elements;
            int rows = matrix.RowCount;
            int cols = matrix.ColumnCount;

            for (int i = 0; i < rows; i++)
            {
                fixed (double* dataPtr = element[i])
                {
                    double* ptr = dataPtr;

                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(*ptr);
                        ptr++;
                    }
                }
            }
        }
    }
}