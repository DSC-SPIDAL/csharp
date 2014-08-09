using System.IO;
using System.IO.Compression;

namespace Salsa.Core
{
    public static class ArrayExtensions
    {
        public static byte[] Compress(double[] toCompress)
        {
            int count = toCompress.Length;
            using (var ms = new MemoryStream())
            {
                using (var gz = new GZipStream(ms, CompressionMode.Compress, true))
                {
                    var writer = new BinaryWriter(gz);
                    writer.Write(count);

                    for (int i = 0; i < count; i++)
                    {
                        writer.Write(toCompress[i]);
                    }

                    writer.Flush();
                    gz.Flush();
                }
                ms.Flush();
                return ms.ToArray();
            }
        }

        public static double[] DecompressDouble(byte[] toDecompress)
        {
            double[] result = null;

            using (var ms = new MemoryStream(toDecompress))
            {
                using (var gz = new GZipStream(ms, CompressionMode.Decompress, true))
                {
                    using (var br = new BinaryReader(gz))
                    {
                        int count = br.ReadInt32();
                        result = new double[count];

                        for (int i = 0; i < count; i++)
                        {
                            result[i] = br.ReadDouble();
                        }
                    }
                }
            }

            return result;
        }
    }
}