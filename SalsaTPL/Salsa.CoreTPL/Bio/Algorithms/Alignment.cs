#region NAligner Copyright

/*
 * NAligner
 * C# port of JAligner API, http://jaligner.sourceforge.net
 * 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#endregion

using System;
using System.Text;

namespace Salsa.Core.Bio.Algorithms
{
    /// <summary> Holds the output of a pairwise sequences alignment. </summary>
    [Serializable]
    public class PairwiseAlignment
    {
        #region Members

        /// <summary> Gap character</summary>
        public const char GAP = '-';

        /// <summary> Default name for sequence #1</summary>
        private const String DEFAULT_SEQUENCE1_NAME = "swg_1";

        /// <summary> Default name for sequence #2</summary>
        private const String DEFAULT_SEQUENCE2_NAME = "swg_2";

        /// <summary> Gap extend cost</summary>
        private float extend;

        /// <summary> Count of gap locations</summary>
        private int gaps;

        /// <summary> Count of identical locations</summary>
        private int identity;

        /// <summary> Markup line</summary>
        private char[] markupLine;

        /// <summary> Scoring matrix</summary>
        private string matrix;

        /// <summary> Name of sequence #1</summary>
        private String name1;

        /// <summary> Name of sequence #2</summary>
        private String name2;

        /// <summary> Gap open cost</summary>
        private float open;

        /// <summary> Alignment score</summary>
        private float score;

        /// <summary> Aligned sequence #1</summary>
        private char[] sequence1;

        /// <summary> Aligned sequence #2</summary>
        private char[] sequence2;

        /// <summary> Count of similar locations</summary>
        private int similarity;

        /// <summary> Alignment start location in sequence #1</summary>
        private int start1;

        /// <summary> Alignment start location in sequence #2</summary>
        private int start2;

        #endregion

        /// <returns> Returns the matrix. </returns>
        public string Matrix
        {
            get { return matrix; }
            set { matrix = value; }
        }

        /// <returns> Returns the name1. </returns>
        public string Name1
        {
            get { return name1 == null || name1.Length == 0 ? DEFAULT_SEQUENCE1_NAME : name1; }
            set { name1 = value; }
        }

        /// <returns> Returns the name2. </returns>
        public string Name2
        {
            get { return name2 == null || name2.Length == 0 ? DEFAULT_SEQUENCE2_NAME : name2; }

            set { name2 = value; }
        }

        /// <returns> Returns the open. </returns>
        public float GapOpenPenalty
        {
            get { return open; }
            set { open = value; }
        }

        /// <returns> Returns the extend. </returns>
        public float GapExtendPenalty
        {
            get { return extend; }
            set { extend = value; }
        }

        /// <returns> Returns the score. </returns>
        public float Score
        {
            get { return score; }
            set { score = value; }
        }

        /// <returns> Returns the sequence1. </returns>
        public char[] Sequence1
        {
            get { return sequence1; }
            set { sequence1 = value; }
        }

        /// <returns> Returns the sequence2. </returns>
        public char[] Sequence2
        {
            get { return sequence2; }
            set { sequence2 = value; }
        }

        /// <returns> Returns the start1. </returns>
        public int Start1
        {
            get { return start1; }
            set { start1 = value; }
        }

        /// <returns> Returns the start2. </returns>
        public int Start2
        {
            get { return start2; }
            set { start2 = value; }
        }

        /// <returns> Returns the gaps. </returns>
        public int Gaps
        {
            get { return gaps; }
            set { gaps = value; }
        }

        /// <returns> Returns the markupLine. </returns>
        public char[] MarkupLine
        {
            get { return markupLine; }
            set { markupLine = value; }
        }

        /// <returns> Returns the identity. </returns>
        public int Identity
        {
            get { return identity; }
            set { identity = value; }
        }

        /// <returns> Returns the similarity. </returns>
        public int Similarity
        {
            get { return similarity; }
            set { similarity = value; }
        }

        public float PercentIdentity
        {
            get { return (Identity/(Sequence1.Length*1.0f)); }
        }

        public float PercentSimilarity
        {
            get { return (Similarity/(Sequence1.Length*1.0f)); }
        }

        public float PercentGaps
        {
            get { return (Gaps/(Sequence1.Length*1.0f)); }
        }

        /// <summary> Returns a summary for alignment</summary>
        public string Summary
        {
            get
            {
                var buffer = new StringBuilder();
                buffer.Append("Sequence #1: " + Name1);
                buffer.AppendLine();
                buffer.Append("Sequence #2: " + Name2);
                buffer.AppendLine();
                buffer.Append("Matrix: " + Matrix);
                buffer.AppendLine();
                buffer.Append("Gap open: " + open);
                buffer.AppendLine();
                buffer.Append("Gap extend: " + extend);
                buffer.AppendLine();
                buffer.Append("Start1: " + start1);
                buffer.AppendLine();
                buffer.Append("Start2: " + start2);
                buffer.AppendLine();
                buffer.Append("Length: " + Sequence1.Length);
                buffer.AppendLine();
                buffer.AppendFormat("Identity: {0}/{1} ({2}) ({3})", identity, Sequence1.Length,
                                    PercentIdentity.ToString("F5"), (1.0f - PercentIdentity).ToString("F5"));
                buffer.AppendLine();
                buffer.AppendFormat("Similarity: {0}/{1} ({2}) ({3})", similarity, Sequence1.Length,
                                    PercentSimilarity.ToString("F5"), (1.0f - PercentSimilarity).ToString("F5"));
                buffer.AppendLine();
                buffer.AppendFormat("JunkesCantor: {0}", ComputeJukesCantorDistance().ToString("F5"));
                buffer.AppendLine();
                buffer.AppendFormat("Kimera2: {0}", ComputeKimuraDistance().ToString("F5"));
                buffer.AppendLine();
                buffer.Append("Gaps: " + gaps + "/" + Sequence1.Length + " (" + PercentGaps.ToString("F5") + ")");
                buffer.AppendLine();
                buffer.Append("Score: " + score.ToString("F2"));
                buffer.AppendLine();
                buffer.Append("Score/Length: " + (score/(sequence1.Length*1.0)).ToString("F3"));
                buffer.AppendLine();
                buffer.AppendFormat(">Name={0} Start={1}", Name1, Start1);
                buffer.AppendLine();
                buffer.Append(Sequence1);
                buffer.AppendLine();
                buffer.AppendFormat(">Name={0} Start={1}", Name2, Start2);
                buffer.AppendLine();
                buffer.Append(Sequence2);
                buffer.AppendLine();
                buffer.Append(MarkupLine);
                buffer.AppendLine();
                return buffer.ToString();
            }
        }

        public double ComputeKimuraDistance()
        {
            double length = 0;
            double gapCount = 0;
            double transitionCount = 0; // P = A -> G | G -> A | C -> T | T -> C
            double transversionCount = 0; // Q = A -> C | A -> T | C -> A | C -> G | T -> A  | T -> G | G -> T | G -> C

            for (int i = 0; i < sequence1.Length; i++)
            {
                length++;
                char nt1 = Sequence1[i]; //nucleotide 1
                char nt2 = Sequence2[i]; //nucelotide 2

                if (nt1 != nt2)
                {
                    // Don't consider gaps at all in this computation;
                    if (nt1 == GAP || nt2 == GAP)
                    {
                        gapCount++;
                    }
                    else if ((nt1 == 'A' && nt2 == 'G') || (nt1 == 'G' && nt2 == 'A') || (nt1 == 'C' && nt2 == 'T') ||
                             (nt1 == 'T' && nt2 == 'C'))
                    {
                        transitionCount++;
                    }
                    else
                    {
                        transversionCount++;
                    }
                }
            }

            double P = transitionCount/(length - gapCount);
            double Q = transversionCount/(length - gapCount);

            double artificialDistance = 10;
            if (1.0 - (2.0*P + Q) <= double.Epsilon)
            {
                PrintArtificialDistanceAlignments(Sequence1, Sequence2, artificialDistance, "Kimura2");
                return artificialDistance;
            }
            if (1.0 - (2.0*Q) <= double.Epsilon)
            {
                PrintArtificialDistanceAlignments(Sequence1, Sequence2, artificialDistance, "Kimura2");
                return artificialDistance;
            }

            return (-0.5*Math.Log(1.0 - 2.0*P - Q) - 0.25*Math.Log(1.0 - 2.0*Q));
        }

        private static void PrintArtificialDistanceAlignments(char[] alignedSeqA, char[] alignedSeqB, double aDistance,
                                                              string distanceType)
        {
            Console.WriteLine("*******************************************");
            Console.WriteLine(alignedSeqB);
            Console.WriteLine(alignedSeqB);
            Console.WriteLine("Artificial " + distanceType + "Distance: " + aDistance);
            Console.WriteLine("*******************************************");
        }

        public double ComputeJukesCantorDistance()
        {
            double length = 0;
            double gapCount = 0;
            double differenceCount = 0;

            for (int i = 0; i < sequence1.Length; i++)
            {
                length++;
                char nt1 = Sequence1[i]; //nucleotide 1
                char nt2 = Sequence2[i]; //nucelotide 2

                if (nt1 != nt2)
                {
                    // Don't consider gaps at all in this computation;
                    if (nt1 == GAP || nt2 == GAP)
                    {
                        gapCount++;
                    }
                    else
                    {
                        differenceCount++;
                    }
                }
            }

            double d = 1.0 - ((4.0/3.0)*(differenceCount/(length - gapCount)));

            if (d <= double.Epsilon)
            {
                throw new Exception("Jukes and Cantor Distance Undefined - Log(Zero)");
            }

            return (-0.75*Math.Log(d));
        }
    }
}