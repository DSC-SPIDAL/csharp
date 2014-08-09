using System.IO;

namespace Salsa.Core.Blas
{
    public class MatrixTextWriter
    {
        private readonly char _delimiter;

        public MatrixTextWriter()
            : this(',')
        {
        }

        public MatrixTextWriter(char delimiter)
        {
            _delimiter = delimiter;
        }

        public void Write(Matrix<short> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<ushort> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<int> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<uint> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<float> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<double> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(Matrix<short> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }

        public void Write(Matrix<ushort> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }

        public void Write(Matrix<int> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }

        public void Write(Matrix<uint> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }

        public void Write(Matrix<float> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }

        public void Write(Matrix<double> matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    if (j > 0)
                    {
                        writer.Write(_delimiter);
                    }

                    writer.Write(matrix.Elements[i][j].ToString());
                }

                writer.WriteLine();
            }
        }
    }
}