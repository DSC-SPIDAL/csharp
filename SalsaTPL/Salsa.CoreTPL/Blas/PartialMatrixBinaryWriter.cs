using System.IO;

namespace Salsa.Core.Blas
{
    public class PartialMatrixBinaryWriter
    {
        public void Write(PartialMatrix<short> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<ushort> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<int> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<uint> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<float> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<double> matrix, string fileName)
        {
            using (Stream stream = File.Create(fileName))
            {
                using (var writer = new BinaryWriter(stream))
                {
                    Write(matrix, writer);
                }
            }
        }

        public void Write(PartialMatrix<short> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }

        public void Write(PartialMatrix<ushort> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }

        public void Write(PartialMatrix<int> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }

        public void Write(PartialMatrix<uint> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }

        public void Write(PartialMatrix<float> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }

        public void Write(PartialMatrix<double> matrix, BinaryWriter writer)
        {
            for (int i = 0; i < matrix.RowCount; i++)
            {
                for (int j = 0; j < matrix.ColumnCount; j++)
                {
                    writer.Write(matrix.Elements[i][j]);
                }
            }
        }
    }
}