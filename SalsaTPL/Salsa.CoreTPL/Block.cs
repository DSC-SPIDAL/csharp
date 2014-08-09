using System;

namespace Salsa.Core
{
    /// <summary>
    /// Represents a block within a 2D array
    /// </summary>
    [Serializable]
    public class Block
    {
        public readonly Range ColumnRange;
        public readonly Range RowRange;
        private int _columnBlockNumber;
        private int _rowBlockNumber;

        /// <summary>
        /// Initializes a new instance of the <see cref="BlockPartition"/> class.
        /// </summary>
        /// <param name="rowRange">The row range.</param>
        /// <param name="colRange">The col range.</param>
        public Block(Range rowRange, Range colRange)
        {
            RowRange = rowRange;
            ColumnRange = colRange;
        }

        public bool IsDiagonal { get; set; }

        public int RowBlockNumber
        {
            get { return _rowBlockNumber; }
            set { _rowBlockNumber = value; }
        }

        public int ColumnBlockNumber
        {
            get { return _columnBlockNumber; }
            set { _columnBlockNumber = value; }
        }

        public void SetIndex(int rowBlockNumber, int columnBlockNumber)
        {
            _rowBlockNumber = rowBlockNumber;
            _columnBlockNumber = columnBlockNumber;
        }

        /// <summary>
        /// Transposes the row and column ranges
        /// </summary>
        /// <returns></returns>
        public Block Transpose()
        {
            var b = new Block(ColumnRange, RowRange);
            return b;
        }

        public override string ToString()
        {
            return string.Format("[{0} {1}]", RowRange, ColumnRange);
        }
    }
}