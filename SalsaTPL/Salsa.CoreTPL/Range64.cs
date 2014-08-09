using System;

namespace Salsa.Core
{
    /// <summary>
    /// Represents a range within an 1D array.
    /// </summary>
    [Serializable]
    public sealed class Range64
    {
        /// <summary>
        /// The inclusive ending index of the BlockPartition.
        /// </summary>
        /// <value>The rangeEnd index.</value>
        public readonly long EndIndex;

        /// <summary>
        /// The total length of the BlockPartition.
        /// </summary>
        /// <value>The length.</value>
        public readonly long Length;

        /// <summary>
        /// The inclusive starting index of the BlockPartition.
        /// </summary>
        /// <value>The rangeStart index.</value>
        public readonly long StartIndex;

        /// <summary>
        /// Initializes a new instance of the <see cref="BlockPartition"/> class.
        /// </summary>
        /// <param name="rangeStart">The starting index of the Range.</param>
        /// <param name="rangeEnd">The ending index of the Range.</param>
        public Range64(long start, long end)
        {
            StartIndex = start;
            EndIndex = end;
            Length = end - start + 1L;
        }

        public bool Contains(long index)
        {
            return (index >= StartIndex && index <= EndIndex);
        }

        /// <summary>
        /// Returns the fully qualified type name of this instance.
        /// </summary>
        /// <returns>
        /// A <see cref="T:System.String"/> containing a fully qualified type name.
        /// </returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", StartIndex.ToString(), EndIndex.ToString());
        }
    }
}