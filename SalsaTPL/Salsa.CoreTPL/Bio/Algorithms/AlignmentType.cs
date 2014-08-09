namespace Salsa.Core.Bio.Algorithms
{
    public enum AlignmentType
    {
        Protein,
        Nucleic
    }

    public enum DistanceFunctionType
    {
        PercentIdentity = 0,
        Kimura2 = 1,
        JukesCantor = 2,
        MinMaxNormScore = 3
    }
}