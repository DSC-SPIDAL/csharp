using System;

namespace Salsa.PairwiseClusteringTPL
{
    /// <summary>
    /// Simple MPI Serializable packet for EM Loop
    /// </summary>
    [Serializable()]
    public class MPIPacket<T>
    {
        public int FirstPoint;  
        public int NumberofPoints;
        public T[] Marray;

        public MPIPacket(int Maxlength)
        {
            Marray = new T[Maxlength];
        }

        public void Clear()
        {
            FirstPoint = 0;
            NumberofPoints = 0;
            Array.Clear(Marray, 0, Marray.Length);
        }
    }
}
