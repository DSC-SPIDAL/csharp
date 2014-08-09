using System;
using System.Text;

namespace Salsa.Core.Blas
{
    [Serializable]
    public class Matrix<T>
    {
        private readonly int _colCount;
        private readonly T[][] _elements;
        private readonly int _rowCount;

        public Matrix(T[][] values)
        {
            _elements = values;
            GetRowColumnCount(_elements, out _rowCount, out _colCount);
        }

        public Matrix(int rowCount, int columnCount)
        {
            _rowCount = rowCount;
            _colCount = columnCount;
            _elements = CreateElements(_rowCount, _colCount);
        }

        public T this[int rowIndex, int columnIndex]
        {
            get { return _elements[rowIndex][columnIndex]; }
            set { _elements[rowIndex][columnIndex] = value; }
        }

        public int RowCount
        {
            get { return _rowCount; }
        }

        public int ColumnCount
        {
            get { return _colCount; }
        }

        public bool IsSquare
        {
            get { return (_rowCount == _colCount); }
        }

        internal T[][] Elements
        {
            get { return _elements; }
        }

        public Matrix<T> Transpose()
        {
            T[][] leftElements = CreateElements(_colCount, _rowCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    leftElements[j][i] = _elements[i][j];
                }
            }

            return new Matrix<T>(leftElements);
        }

        public void SetAllValues(T value)
        {
            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    _elements[i][j] = value;
                }
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            for (int i = 0; i < Math.Min(80, _rowCount); i++)
            {
                for (int j = 0; j < Math.Min(80, _colCount); j++)
                {
                    if (j > 0)
                    {
                        sb.Append(",");
                    }

                    sb.AppendFormat("{0}", _elements[i][j]);
                }

                sb.AppendLine();
            }

            return sb.ToString();
        }

        public static implicit operator T[][](Matrix<T> matrix)
        {
            return matrix.Elements;
        }

        internal static T[][] CreateElements(int rowCount, int columnCount)
        {
            var elements = new T[rowCount][];

            for (int i = 0; i < rowCount; i++)
            {
                elements[i] = new T[columnCount];
            }

            return elements;
        }

        internal static void GetRowColumnCount(T[][] data, out int rows, out int columns)
        {
            rows = data.Length;
            columns = (rows == 0) ? 0 : data[0].Length;
        }

        #region Row Operations

        public T[] GetRowValues(int rowIndex)
        {
            var values = new T[_colCount];
            Array.Copy(_elements[rowIndex], values, _colCount);
            return values;
        }

        public void SetRowValues(int rowIndex, T[] values)
        {
            Array.Copy(values, _elements[rowIndex], values.Length);
        }

        #endregion

        #region Column Operations

        public T[] GetColumnValues(int columnIndex)
        {
            var leftElements = new T[_rowCount];

            for (int i = 0; i < _rowCount; i++)
            {
                leftElements[i] = _elements[i][columnIndex];
            }

            return leftElements;
        }

        public void SetColumnValues(int columnIndex, T[] values)
        {
            for (int i = 0; i < values.Length; i++)
            {
                _elements[i][columnIndex] = values[i];
            }
        }

        #endregion

        #region Block Operations

        public T[][] GetBlockValues(Block block)
        {
            return GetBlockValues(block.RowRange, block.ColumnRange);
        }

        public T[][] GetBlockValues(Range rowRange, Range columnRange)
        {
            return GetBlockValues(rowRange.StartIndex, rowRange.EndIndex, columnRange.StartIndex, columnRange.EndIndex);
        }

        public T[][] GetBlockValues(int rowStartIndex, int rowEndIndex, int columnStartIndex, int columnEndIndex)
        {
            T[][] leftElements = CreateElements(rowEndIndex - rowStartIndex + 1, columnEndIndex - columnStartIndex + 1);

            for (int i = rowStartIndex; i <= rowEndIndex; i++)
            {
                for (int j = columnStartIndex; j <= columnEndIndex; j++)
                {
                    leftElements[i - rowStartIndex][j - columnStartIndex] = _elements[i][j];
                }
            }


            return leftElements;
        }

        public void SetBlockValues(Block block, T[][] values)
        {
            SetBlockValues(block.RowRange, block.ColumnRange, values);
        }

        public void SetBlockValues(Range rowRange, Range columnRange, T[][] values)
        {
            SetBlockValues(rowRange.StartIndex, rowRange.EndIndex, columnRange.StartIndex, columnRange.EndIndex, values);
        }

        public void SetBlockValues(int rowStartIndex, int rowEndIndex, int columnStartIndex, int columnEndIndex,
                                   T[][] values)
        {
            for (int i = rowStartIndex; i <= rowEndIndex; i++)
            {
                for (int j = columnStartIndex; j <= columnEndIndex; j++)
                {
                    _elements[i][j] = values[i - rowStartIndex][j - columnStartIndex];
                }
            }
        }

        #endregion
    }
}