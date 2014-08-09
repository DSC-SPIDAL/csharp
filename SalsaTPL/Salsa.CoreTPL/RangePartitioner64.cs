using System;

namespace Salsa.Core
{
    public static class RangePartitioner64
    {
        public static Range64[] Partition(long length, long numPartitions)
        {
            return Partition(0, length, numPartitions);
        }

        public static Range64[] Partition(Range64 range, long numPartitions)
        {
            return Partition(range.StartIndex, range.Length, numPartitions);
        }

        public static Range64[] Partition(long startIndex, long length, long numPartitions)
        {
            if (numPartitions < 1)
            {
                throw new ArgumentOutOfRangeException("count",
                                                      "Partitioning requires numPartitions to be greater than zero.");
            }

            if (length < 1)
            {
                throw new ArgumentOutOfRangeException("length", "Partitioning requires length to be greater than zero.");
            }

            if (length < numPartitions)
            {
                throw new InvalidOperationException(
                    "Partitioning cannot be performed when length is less than numPartitions requested.");
            }

            var ranges = new Range64[numPartitions];
            long chunksize = length/numPartitions;
            long remainder = length%numPartitions;

            for (int i = 0; i < numPartitions; i++)
            {
                if (remainder > 0)
                {
                    ranges[i] = new Range64(startIndex, startIndex + chunksize);
                    startIndex += chunksize + 1;
                    remainder--;
                }
                else
                {
                    ranges[i] = new Range64(startIndex, startIndex + chunksize - 1);
                    startIndex += chunksize;
                }
            }

            return ranges;
        }

        public static Range64[] PartitionByLength(long length, long maxPartitionLength)
        {
            return PartitionByLength(0, length, maxPartitionLength);
        }

        public static Range64[] PartitionByLength(Range64 range, long maxPartitionLength)
        {
            return PartitionByLength(range.StartIndex, range.Length, maxPartitionLength);
        }

        public static Range64[] PartitionByLength(long startIndex, long length, long maxPartitionLength)
        {
            if (maxPartitionLength < 1)
            {
                throw new ArgumentException("maxPartitionLength",
                                            "Partitioning requires the maxPartitionLength to be greater than zero.");
            }
            if (length < 1)
            {
                throw new ArgumentException("length", "Partitioning requires the length to be greater than zero.");
            }
            if (length < maxPartitionLength)
            {
                throw new ArgumentException("length",
                                            "Partitioning requires the length to be greater than maxPartitionLength.");
            }

            long rangeCount = length/maxPartitionLength;
            long lastRangeLength = length%maxPartitionLength;

            if (lastRangeLength > 0)
            {
                rangeCount++;
            }
            else
            {
                lastRangeLength = maxPartitionLength;
            }

            var ranges = new Range64[rangeCount];

            for (long i = 0; i < rangeCount; i++)
            {
                long start = i*maxPartitionLength;
                long end = (i == rangeCount - 1 ? start + lastRangeLength - 1 : start + maxPartitionLength - 1);
                ranges[i] = new Range64(start, end);
                start += maxPartitionLength;
            }

            return ranges;
        }
    }
}