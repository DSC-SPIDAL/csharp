namespace Salsa.Core
{
    public static class BlockPartitioner
    {
        public static Block[][] Partition(int rowCount, int columnCount, int numRowBlocks, int numColumnBlocks)
        {
            return Partition(0, 0, rowCount, columnCount, numRowBlocks, numColumnBlocks);
        }

        public static Block[][] Partition(int rowStartIndex, int columnStartIndex, int rowCount, int columnCount,
                                          int numRowBlocks, int numColumnBlocks)
        {
            Range[] rowRanges = RangePartitioner.Partition(rowStartIndex, rowCount, numRowBlocks);
            Range[] colRanges = RangePartitioner.Partition(columnStartIndex, columnCount, numColumnBlocks);
            var result = new Block[numRowBlocks][];

            for (int i = 0; i < rowRanges.Length; i++)
            {
                result[i] = new Block[numColumnBlocks];

                for (int j = 0; j < colRanges.Length; j++)
                {
                    result[i][j] = new Block(rowRanges[i], colRanges[j]);
                }
            }

            return result;
        }

        public static Block[][] PartitionByLength(int rowCount, int columnCount, int maxRowPartitionLength,
                                                  int maxColumnPartitionLength)
        {
            return PartitionByLength(0, 0, rowCount, columnCount, maxRowPartitionLength, maxColumnPartitionLength);
        }

        public static Block[][] PartitionByLength(int rowStartIndex, int columnStartIndex, int rowCount, int columnCount,
                                                  int maxRowPartitionLength, int maxColumnPartitionLength)
        {
            Range[] rowRanges = RangePartitioner.PartitionByLength(rowStartIndex, rowCount, maxRowPartitionLength);
            Range[] colRanges = RangePartitioner.PartitionByLength(columnStartIndex, columnCount,
                                                                   maxColumnPartitionLength);
            var result = new Block[rowRanges.Length][];

            for (int i = 0; i < rowRanges.Length; i++)
            {
                result[i] = new Block[colRanges.Length];

                for (int j = 0; j < colRanges.Length; j++)
                {
                    result[i][j] = new Block(rowRanges[i], colRanges[j]);
                }
            }

            return result;
        }
    }
}