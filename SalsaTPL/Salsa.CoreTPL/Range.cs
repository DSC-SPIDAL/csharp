using System;

namespace Salsa.Core
{
    /// <summary>
    /// Represents a range within an 1D array.
    /// </summary>
    [Serializable]
    public sealed class Range
    {
        /// <summary>
        /// The inclusive ending index of the BlockPartition.
        /// </summary>
        /// <value>The rangeEnd index.</value>
        public readonly int EndIndex;

        /// <summary>
        /// The total length of the BlockPartition.
        /// </summary>
        /// <value>The length.</value>
        public readonly int Length;

        /// <summary>
        /// The inclusive starting index of the BlockPartition.
        /// </summary>
        /// <value>The rangeStart index.</value>
        public readonly int StartIndex;

        /// <summary>
        /// Initializes a new instance of the <see cref="BlockPartition"/> class.
        /// </summary>
        /// <param name="rangeStart">The starting index of the Range.</param>
        /// <param name="rangeEnd">The ending index of the Range.</param>
        public Range(int start, int end)
        {
            StartIndex = start;
            EndIndex = end;
            Length = end - start + 1;
        }

        public string StartSeqName { get; set; }
        public string EndSeqName { get; set; }

        public bool Contains(int index)
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
            return string.Format("({0}:{1})", StartIndex.ToString(), EndIndex.ToString());
        }

        /// <summary>
        /// Returns true if there's an intersection of this range with the given range
        /// </summary>
        /// <param name="range">The range to see if an intersection exists</param>
        /// <returns>True if an intersection exists or false otherwise</returns>
        public bool IntersectsWith(Range range)
        {
            Range lengthiest = range.Length >= Length ? range : this;
            Range other = range == lengthiest ? this : range;

            return lengthiest.Contains(other.StartIndex) || lengthiest.Contains(other.EndIndex);
        }

        /// <summary>
        /// Gets the intersecting range assuming an intersection exists. Use <code>IntersectsWith(Range range)</code>
        /// to check for an existing intersection
        /// </summary>
        /// <param name="range">The range to intersect with</param>
        /// <returns>The intersection with the given range</returns>
        public Range GetIntersectionWith(Range range)
        {
            int start = range.StartIndex >= StartIndex ? range.StartIndex : StartIndex;
            int end = range.EndIndex <= EndIndex ? range.EndIndex : EndIndex;
            if (start <= end)
                return new Range(start, end);
            throw new ApplicationException(
                String.Format("Given range [{0}, {1}] does not intersect with this range [{2}, {3}]", range.StartIndex,
                              range.EndIndex, StartIndex, EndIndex));
        }
    }
}