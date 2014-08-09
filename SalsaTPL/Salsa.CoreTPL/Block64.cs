using System;

namespace Salsa.Core
{
    /// <summary>
    /// Represents a block within a 2D array
    /// </summary>
    [Serializable]
    public class Block64
    {
        public readonly Range64 ColumnRange;
        public readonly Range64 RowRange;

        public Block64(Range64 rowRange, Range64 colRange)
        {
            RowRange = rowRange;
            ColumnRange = colRange;
        }

        /// <summary>
        /// Transposes the row and column ranges
        /// </summary>
        /// <returns></returns>
        public Block64 Transpose()
        {
            var b = new Block64(ColumnRange, RowRange);
            return b;
        }

        public override string ToString()
        {
            return string.Format("[{0} {1}]", RowRange, ColumnRange);
        }
    }
}