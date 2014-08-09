using System;

namespace Salsa.Core
{
    public static class RangePartitioner
    {
        public static Range[] Partition(int length, int numPartitions)
        {
            return Partition(0, length, numPartitions);
        }

        public static Range[] Partition(Range range, int numPartitions)
        {
            return Partition(range.StartIndex, range.Length, numPartitions);
        }

        public static Range[] Partition(int startIndex, int length, int numPartitions)
        {
            if (numPartitions < 1)
            {
                throw new ArgumentOutOfRangeException("count",
                                                      "Partitioning requires numPartitions to be greater than zero.");
            }

            if (length < 1)
            {
//                throw new ArgumentOutOfRangeException("length", "Partitioning requires length to be greater than zero.");
            }

            if (length < numPartitions)
            {
                //              throw new InvalidOperationException("Partitioning cannot be performed when length is less than numPartitions requested.");
            }

            var ranges = new Range[numPartitions];
            int chunksize = length/numPartitions;
            int remainder = length%numPartitions;

            for (int i = 0; i < numPartitions; i++)
            {
                if (remainder > 0)
                {
                    ranges[i] = new Range(startIndex, startIndex + chunksize);
                    startIndex += chunksize + 1;
                    remainder--;
                }
                else
                {
                    ranges[i] = new Range(startIndex, startIndex + chunksize - 1);
                    startIndex += chunksize;
                }
            }

            return ranges;
        }

        public static Range[] PartitionByLength(int length, int maxPartitionLength)
        {
            return PartitionByLength(0, length, maxPartitionLength);
        }

        public static Range[] PartitionByLength(Range range, int maxPartitionLength)
        {
            return PartitionByLength(range.StartIndex, range.Length, maxPartitionLength);
        }

        public static Range[] PartitionByLength(int startIndex, int length, int maxPartitionLength)
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

            int rangeCount = length/maxPartitionLength;
            int lastRangeLength = length%maxPartitionLength;

            if (lastRangeLength > 0)
            {
                rangeCount++;
            }
            else
            {
                lastRangeLength = maxPartitionLength;
            }

            var ranges = new Range[rangeCount];

            for (int i = 0; i < rangeCount; i++)
            {
                int start = i*maxPartitionLength;
                int end = (i == rangeCount - 1 ? start + lastRangeLength - 1 : start + maxPartitionLength - 1);
                ranges[i] = new Range(start, end);
                start += maxPartitionLength;
            }

            return ranges;
        }
    }
}