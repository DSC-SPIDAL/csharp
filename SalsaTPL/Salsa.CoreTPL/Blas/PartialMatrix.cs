using System;
using System.Text;

namespace Salsa.Core.Blas
{
    [Serializable]
    public class PartialMatrix<T>
    {
        private readonly int _colCount;
        private readonly T[][] _elements;
        private readonly int _globalColStartIndex;
        private readonly int _globalRowStartIndex;
        private readonly int _rowCount;

        public PartialMatrix(Range globalRowRange, Range globalColumnRange)
            : this(
                globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                globalColumnRange.EndIndex)
        {
        }

        public PartialMatrix(Range globalRowRange, int globalColumnStartIndex, int globalColumnEndIndex)
            : this(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnStartIndex, globalColumnEndIndex)
        {
        }

        public PartialMatrix(int globalRowStartIndex, int globalRowEndIndex, Range globalColumnRange)
            : this(globalRowStartIndex, globalRowEndIndex, globalColumnRange.StartIndex, globalColumnRange.EndIndex)
        {
        }

        /// <summary>
        /// Creates a Partial Matrix
        /// </summary>
        /// <param name="rowStartIndex">Inclusive row starting index</param>
        /// <param name="rowEndIndex">Inclusive row ending index</param>
        /// <param name="columnStartIndex">Inclusive column starting index</param>
        /// <param name="columnEndIndex">Inclusive column ending index</param>
        public PartialMatrix(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                             int globalColumnEndIndex)
        {
            _globalRowStartIndex = globalRowStartIndex;
            _globalColStartIndex = globalColumnStartIndex;
            _rowCount = globalRowEndIndex - globalRowStartIndex + 1;
            _colCount = globalColumnEndIndex - globalColumnStartIndex + 1;
            _elements = CreateElements(_rowCount, _colCount);
        }

        internal PartialMatrix(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                               int globalColumnEndIndex, T[][] values)
        {
            _globalRowStartIndex = globalRowStartIndex;
            _globalColStartIndex = globalColumnStartIndex;
            _rowCount = globalRowEndIndex - globalRowStartIndex + 1;
            _colCount = globalColumnEndIndex - globalColumnStartIndex + 1;
            _elements = values;
        }

        public T this[int globalRowIndex, int globalColumnIndex]
        {
            get { return _elements[globalRowIndex - _globalRowStartIndex][globalColumnIndex - _globalColStartIndex]; }
            set { _elements[globalRowIndex - _globalRowStartIndex][globalColumnIndex - _globalColStartIndex] = value; }
        }

        public int RowCount
        {
            get { return _rowCount; }
        }

        public int ColumnCount
        {
            get { return _colCount; }
        }

        public int GlobalRowStartIndex
        {
            get { return _globalRowStartIndex; }
        }

        public int GlobalRowEndIndex
        {
            get { return _globalRowStartIndex + _rowCount - 1; }
        }

        public int GlobalColumnStartIndex
        {
            get { return _globalColStartIndex; }
        }

        public int GlobalColumnEndIndex
        {
            get { return _globalColStartIndex + _colCount - 1; }
        }

        public bool IsSquare
        {
            get { return (_rowCount == _colCount); }
        }

        internal T[][] Elements
        {
            get { return _elements; }
        }

        public PartialMatrix<T> Transpose()
        {
            T[][] leftElements = CreateElements(_colCount, _rowCount);

            for (int i = 0; i < _rowCount; i++)
            {
                for (int j = 0; j < _colCount; j++)
                {
                    leftElements[j][i] = _elements[i][j];
                }
            }

            return new PartialMatrix<T>(GlobalColumnStartIndex, GlobalColumnEndIndex, GlobalRowStartIndex,
                                        GlobalRowEndIndex, leftElements);
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

        public void CopyTo(Matrix<T> matrix)
        {
            matrix.SetBlockValues(GlobalRowStartIndex, GlobalRowEndIndex, GlobalColumnStartIndex, GlobalColumnEndIndex,
                                  _elements);
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            for (int i = 0; i < Math.Min(20, _rowCount); i++)
            {
                for (int j = 0; j < Math.Min(20, _colCount); j++)
                {
                    if (j > 0)
                    {
                        sb.Append(", ");
                    }

                    sb.AppendFormat("{0,4}", _elements[i][j]);

                    if (i == 20 - 1)
                    {
                        sb.Append("...");
                    }
                }

                sb.AppendLine();
            }

            return sb.ToString();
        }

        public static implicit operator T[][](PartialMatrix<T> matrix)
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

        public T[] GetRowValues(int globalRowIndex)
        {
            globalRowIndex = globalRowIndex - _globalRowStartIndex;
            var values = new T[_colCount];
            Array.Copy(_elements[globalRowIndex], values, _colCount);
            return values;
        }

        public void SetRowValues(int globalRowIndex, T[] values)
        {
            globalRowIndex = globalRowIndex - _globalRowStartIndex;
            Array.Copy(values, _elements[globalRowIndex], values.Length);
        }

        #endregion

        #region Column Operations

        public T[] GetColumnValues(int globalColumnIndex)
        {
            globalColumnIndex = globalColumnIndex - _globalColStartIndex;
            var leftElements = new T[_rowCount];

            for (int i = 0; i < _rowCount; i++)
            {
                leftElements[i] = _elements[i][globalColumnIndex];
            }

            return leftElements;
        }

        public void SetColumnValues(int globalColumnIndex, T[] values)
        {
            globalColumnIndex = globalColumnIndex - _globalColStartIndex;
            for (int i = 0; i < values.Length; i++)
            {
                _elements[i][globalColumnIndex] = values[i];
            }
        }

        #endregion

        #region Block Operations

        public T[][] GetBlockValues(Block globalBlock)
        {
            return GetBlockValues(globalBlock.RowRange, globalBlock.ColumnRange);
        }

        public T[][] GetBlockValues(Range globalRowRange, Range globalColumnRange)
        {
            return GetBlockValues(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                                  globalColumnRange.EndIndex);
        }

        public T[][] GetBlockValues(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                                    int globalColumnEndIndex)
        {
            globalRowStartIndex = globalRowStartIndex - _globalRowStartIndex;
            globalRowEndIndex = globalRowEndIndex - _globalRowStartIndex;
            globalColumnStartIndex = globalColumnStartIndex - _globalColStartIndex;
            globalColumnEndIndex = globalColumnEndIndex - _globalColStartIndex;

            T[][] leftElements = CreateElements(globalRowEndIndex - globalRowStartIndex + 1,
                                                globalColumnEndIndex - globalColumnStartIndex + 1);

            for (int i = globalRowStartIndex; i <= globalRowEndIndex; i++)
            {
                for (int j = globalColumnStartIndex; j <= globalColumnEndIndex; j++)
                {
                    leftElements[i - globalRowStartIndex][j - globalColumnStartIndex] = _elements[i][j];
                }
            }


            return leftElements;
        }


        /// <summary>
        /// Gets the block at the given coordinate but transforms the block before returning it.  This is an optimization
        /// method.
        /// </summary>
        /// <param name="block">The block to retrieve</param>
        /// <returns></returns>
        public T[][] GetTransposedBlock(Block globalBlock)
        {
            return GetTransposedBlock(globalBlock.RowRange, globalBlock.ColumnRange);
        }

        /// <summary>
        /// Gets the block at the given coordinate but transforms the block before returning it.  This is an optimization
        /// method.
        /// </summary>
        /// <param name="rowRange">The range specifiying the rows to retrieve</param>
        /// <param name="columnRange">The range specifiying the columns to retrieve</param>
        /// <returns></returns>
        public T[][] GetTransposedBlock(Range globalRowRange, Range globalColumnRange)
        {
            return GetTransposedBlock(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                                      globalColumnRange.EndIndex);
        }

        /// <summary>
        /// Gets the block at the given coordinate but transforms the block before returning it.  This is an optimization
        /// method.
        /// </summary>
        /// <param name="rowStartIndex">Inclusive starting row index.</param>
        /// <param name="rowEndIndex">Inclusive ending row index.</param>
        /// <param name="columnStartIndex">Inclusive starting column index.</param>
        /// <param name="columnEndIndex">Inclusive ending column index.</param>
        /// <returns></returns>
        public T[][] GetTransposedBlock(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                                        int globalColumnEndIndex)
        {
            globalRowStartIndex = globalRowStartIndex - _globalRowStartIndex;
            globalRowEndIndex = globalRowEndIndex - _globalRowStartIndex;
            globalColumnStartIndex = globalColumnStartIndex - _globalColStartIndex;
            globalColumnEndIndex = globalColumnEndIndex - _globalColStartIndex;

            T[][] leftElements = CreateElements(globalColumnEndIndex - globalColumnStartIndex + 1,
                                                globalRowEndIndex - globalRowStartIndex + 1);

            for (int i = globalRowStartIndex; i <= globalRowEndIndex; i++)
            {
                for (int j = globalColumnStartIndex; j <= globalColumnEndIndex; j++)
                {
                    leftElements[j - globalColumnStartIndex][i - globalRowStartIndex] = _elements[i][j];
                }
            }

            return leftElements;
        }

        public void SetBlockValues(Block globalBlock, T[][] blockValues)
        {
            SetBlockValues(globalBlock.RowRange, globalBlock.ColumnRange, blockValues);
        }

        public void SetBlockValues(Range globalRowRange, Range globalColumnRange, T[][] blockValues)
        {
            SetBlockValues(globalRowRange.StartIndex, globalRowRange.EndIndex, globalColumnRange.StartIndex,
                           globalColumnRange.EndIndex, blockValues);
        }

        public void SetBlockValues(int globalRowStartIndex, int globalRowEndIndex, int globalColumnStartIndex,
                                   int globalColumnEndIndex, T[][] blockValues)
        {
            globalRowStartIndex = globalRowStartIndex - _globalRowStartIndex;
            globalRowEndIndex = globalRowEndIndex - _globalRowStartIndex;
            globalColumnStartIndex = globalColumnStartIndex - _globalColStartIndex;
            globalColumnEndIndex = globalColumnEndIndex - _globalColStartIndex;

            for (int i = globalRowStartIndex; i <= globalRowEndIndex; i++)
            {
                for (int j = globalColumnStartIndex; j <= globalColumnEndIndex; j++)
                {
                    try
                    {
                        _elements[i][j] = blockValues[i - globalRowStartIndex][j - globalColumnStartIndex];
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("{0},{1}", i, j);
                        throw (ex);
                    }
                }
            }
        }

        #endregion
    }
}