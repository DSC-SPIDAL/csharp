using System.IO;

namespace Salsa.Core.Blas
{
    public class PartialMatrixTextWriter
    {
        private readonly char _delimiter;

        public PartialMatrixTextWriter() : this(',')
        {
        }

        public PartialMatrixTextWriter(char delimiter)
        {
            _delimiter = delimiter;
        }

        public void Write(PartialMatrix<short> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<ushort> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<int> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<uint> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<float> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<double> matrix, string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                Write(matrix, writer);
            }
        }

        public void Write(PartialMatrix<short> matrix, StreamWriter writer)
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

        public void Write(PartialMatrix<ushort> matrix, StreamWriter writer)
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

        public void Write(PartialMatrix<int> matrix, StreamWriter writer)
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

        public void Write(PartialMatrix<uint> matrix, StreamWriter writer)
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

        public void Write(PartialMatrix<float> matrix, StreamWriter writer)
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

        public void Write(PartialMatrix<double> matrix, StreamWriter writer)
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