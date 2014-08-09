using System;
using SALSALibrary;

namespace Manxcat
{
    /// <summary>
    /// Simple MPI Serializable packet for single double dimensional array
    /// </summary>
    [Serializable]
    public class MPI2DDoubleVectorPacket
    {
        public int FirstPoint;
        public double[][] Marray;
        public int NumberofPoints;

        public MPI2DDoubleVectorPacket(int Maxlength, int VectorDimension)
        {
            Marray = new double[Maxlength][];
            for (int LongIndex = 0; LongIndex < Maxlength; LongIndex++)
            {
                Marray[LongIndex] = new double[VectorDimension];
            }
        }

        public static void CopyMPI2DDoubleVectorPacket(MPI2DDoubleVectorPacket ToObject,
                                                       MPI2DDoubleVectorPacket FromObject)
        {
            ToObject.FirstPoint = FromObject.FirstPoint;
            ToObject.NumberofPoints = FromObject.NumberofPoints;
            SALSABLAS.CopyVector(ToObject.Marray, FromObject.Marray, 0, FromObject.NumberofPoints);
        }

        // End CopyMPI2DDoubleVectorPacket(MPI2DDoubleVectorPacket ToObject, MPI2DDoubleVectorPacket FromObject)
    }

    [Serializable]
    public class MPI1DStringVectorPacket
    {
        public int FirstPoint;
        public string[] Marray;
        public int NumberofPoints;

        public MPI1DStringVectorPacket(int Maxlength, int VectorDimension)
        {
            Marray = new string[Maxlength];
        }

        public static void CopyMPI1DStringVectorPacket(MPI1DStringVectorPacket ToObject,
                                                       MPI1DStringVectorPacket FromObject)
        {
            ToObject.FirstPoint = FromObject.FirstPoint;
            ToObject.NumberofPoints = FromObject.NumberofPoints;
            SALSABLAS.CopyVector(ToObject.Marray, FromObject.Marray, 0, FromObject.NumberofPoints);
        }

        // End CopyMPI2DDoubleVectorPacket(MPI2DDoubleVectorPacket ToObject, MPI2DDoubleVectorPacket FromObject)
    }

    // End MPI1DString
}

// End Namespace Manxcat