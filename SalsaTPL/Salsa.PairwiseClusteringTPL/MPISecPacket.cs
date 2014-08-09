using System;

namespace Salsa.PairwiseClusteringTPL
{
    [Serializable()]
    public class MPISecPacket
    {
        public int FirstPoint;
        public int NumberofPoints;
        public double[] Marray;
        public double[] Barray;

        public MPISecPacket(int Maxlength1)
        {
            Marray = new double[Maxlength1];
            Barray = new double[Maxlength1];
        }

        public void Clear()
        {
            FirstPoint = 0;
            NumberofPoints = 0;
            Array.Clear(Marray, 0, Marray.Length);
        }

        public static void Membercopy(ref MPISecPacket One, ref MPISecPacket Two)
        {
            Two = (MPISecPacket)One.MemberwiseClone();
            return;
        }
    }
}
