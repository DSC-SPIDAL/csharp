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

namespace Salsa.Core.Bio.Algorithms.SmithWaterman
{
    /// <summary> An implementation of the Smith-Waterman algorithm with
    /// Gotoh's improvement for biological local pairwise sequence alignment. </summary>
    public class SmithWatermanGotoh
    {
        private AlignmentType _alignmentType;
        private float _gapExtensionPenalty;
        private float _gapOpenPenalty;
        private ScoringMatrix _scoringMatrix;

        public SmithWatermanGotoh(AlignmentType alignmentType)
        {
            _alignmentType = alignmentType;
        }

        /// <summary> Aligns two sequences by Smith-Waterman algorithm</summary>
        /// <param name="s1">sequene #1 </param>
        /// <param name="s2">sequene #2 </param>
        /// <param name="matrix">scoring matrix </param>
        /// <param name="o">open gap penalty </param>
        /// <param name="e">extend gap penalty </param>
        /// <returns> alignment object contains the two aligned sequences, 
        /// the alignment score and alignment statistics</returns>
        /// <seealso cref="Sequence"/>
        /// <seealso cref="Matrix"/>
        public PairwiseAlignment Align(Sequence s1, Sequence s2, ScoringMatrix matrix, float gapOpenPenalty,
                                       float gapExtensionPenalty)
        {
            if (s1 == null)
            {
                Console.WriteLine("Error: S1 is null");
                throw new SmithWatermanGotohException("S1 is null");
            }

            if (s2 == null)
            {
                Console.WriteLine("Error: S2 is null");
                throw new SmithWatermanGotohException("S2 is null");
            }

            int m = s1.Length + 1;
            int n = s2.Length + 1;

            var pointers = new Directions[m*n];

            // Initializes the boundaries of the traceback matrix to STOP.
            for (int i = 0, k = 0; i < m; i++, k += n)
            {
                pointers[k] = Directions.STOP;
            }

            for (int j = 1; j < n; j++)
            {
                pointers[j] = Directions.STOP;
            }

            var sizesOfVerticalGaps = new short[m*n];
            var sizesOfHorizontalGaps = new short[m*n];

            for (int i = 0, k = 0; i < m; i++, k += n)
            {
                for (int j = 0; j < n; j++)
                {
                    sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
                }
            }

            _scoringMatrix = matrix;
            _gapOpenPenalty = gapOpenPenalty;
            _gapExtensionPenalty = gapExtensionPenalty;


            Cell cell = Construct(s1, s2, pointers, sizesOfVerticalGaps, sizesOfHorizontalGaps);
            PairwiseAlignment alignment = Traceback(s1.Residues, s2.Residues, pointers, cell, sizesOfVerticalGaps,
                                                    sizesOfHorizontalGaps);

            alignment.Name1 = s1.Label;
            alignment.Name2 = s2.Label;
            alignment.GapOpenPenalty = _gapOpenPenalty;
            alignment.GapExtendPenalty = _gapExtensionPenalty;
            alignment.Matrix = matrix.Name;

            return alignment;
        }

        /// <summary> Constructs directions matrix for the traceback </summary>
        /// <param name="s1">sequence #1 </param>
        /// <param name="s2">sequence #2 </param>
        /// <param name="matrix">scoring matrix </param>
        /// <param name="o">open gap penalty </param>
        /// <param name="e">extend gap penalty </param>
        /// <returns> The cell where the traceback starts. </returns>
        private Cell Construct(Sequence s1, Sequence s2, Directions[] pointers, short[] sizesOfVerticalGaps,
                               short[] sizesOfHorizontalGaps)
        {
            int m = s1.Length + 1;
            int n = s2.Length + 1;

            float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
            var g = new float[n]; // score if xi aligns to a gap after yi
            float h; // score if yi aligns to a gap after xi
            var v = new float[n]; // best score of alignment x1...xi to y1...yi
            float vDiagonal;

            g[0] = float.NegativeInfinity;
            h = float.NegativeInfinity;
            v[0] = 0;

            for (int j = 1; j < n; j++)
            {
                g[j] = float.NegativeInfinity;
                v[j] = 0;
            }

            float[,] scores = _scoringMatrix.Scores;
            float similarityScore, g1, g2, h1, h2;


            var cell = new Cell();

            for (int i = 1, k = n; i < m; i++, k += n)
            {
                h = float.NegativeInfinity;
                vDiagonal = v[0];

                for (int j = 1, l = k + 1; j < n; j++, l++)
                {
                    similarityScore = scores[s1[i - 1], s2[j - 1]];

                    // Fill the matrices
                    f = vDiagonal + similarityScore;

                    g1 = g[j] - _gapExtensionPenalty;
                    g2 = v[j] - _gapOpenPenalty;

                    if (g1 > g2)
                    {
                        // Gap extension penalty
                        g[j] = g1;
                        sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
                    }
                    else
                    {
                        // Gap open penalty
                        g[j] = g2;
                    }

                    h1 = h - _gapExtensionPenalty;
                    h2 = v[j - 1] - _gapOpenPenalty;

                    if (h1 > h2)
                    {
                        // Gap extension penalty
                        h = h1;
                        sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
                    }
                    else
                    {
                        // Gap open penalty
                        h = h2;
                    }

                    vDiagonal = v[j];
                    v[j] = Max(f, g[j], h, 0);

                    // Determine the traceback direction
                    if (v[j] == 0)
                    {
                        pointers[l] = Directions.STOP;
                    }
                    else if (v[j] == f)
                    {
                        pointers[l] = Directions.DIAGONAL;
                    }
                    else if (v[j] == g[j])
                    {
                        pointers[l] = Directions.UP;
                    }
                    else
                    {
                        pointers[l] = Directions.LEFT;
                    }

                    // Set the traceback start at the current cell i, j and score
                    if (v[j] > cell.Score)
                    {
                        cell.Set(i, j, v[j]);
                    }
                }
            }

            return cell;
        }

        /// <summary> Returns the alignment of two sequences based on the passed array of pointers</summary>
        /// <param name="s1">sequence #1 </param>
        /// <param name="s2">sequence #2 </param>
        /// <param name="m">scoring matrix </param>
        /// <param name="cell">The cell where the traceback starts. </param>
        /// <returns> <see cref="Alignment"/> with the two aligned sequences and alignment score. </returns>
        /// <seealso cref="Cell"/>
        /// <seealso cref="Alignment"/>
        private PairwiseAlignment Traceback(string s1, string s2, Directions[] pointers, Cell cell,
                                            short[] sizesOfVerticalGaps, short[] sizesOfHorizontalGaps)
        {
            int n = s2.Length + 1;

            var alignment = new PairwiseAlignment();
            alignment.Score = cell.Score;

            float[,] scores = _scoringMatrix.Scores;

            int maxlen = s1.Length + s2.Length; // maximum length after the aligned sequences

            var reversed1 = new char[maxlen]; // reversed sequence #1
            var reversed2 = new char[maxlen]; // reversed sequence #2
            var reversed3 = new char[maxlen]; // reversed markup


            int len1 = 0; // length of sequence #1 after alignment
            int len2 = 0; // length of sequence #2 after alignment
            int len3 = 0; // length of the markup line

            int identity = 0; // count of identitcal pairs
            int similarity = 0; // count of similar pairs
            int gaps = 0; // count of gaps

            char c1, c2;

            int i = cell.Row; // traceback start row
            int j = cell.Column; // traceback start col
            int k = i*n;

            bool stillGoing = true; // traceback flag: true -> continue & false -> stop

            while (stillGoing)
            {
                switch (pointers[k + j])
                {
                    case Directions.UP:
                        for (int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++)
                        {
                            reversed1[len1++] = s1[--i];
                            reversed2[len2++] = PairwiseAlignment.GAP;
                            reversed3[len3++] = Markups.GAP;
                            k -= n;
                            gaps++;
                        }
                        break;

                    case Directions.DIAGONAL:
                        c1 = s1[--i];
                        c2 = s2[--j];
                        k -= n;
                        reversed1[len1++] = c1;
                        reversed2[len2++] = c2;

                        if (c1 == c2)
                        {
                            reversed3[len3++] = Markups.IDENTITY;
                            identity++;
                            similarity++;
                        }
                        else if (scores[c1, c2] > 0)
                        {
                            reversed3[len3++] = Markups.SIMILARITY;
                            similarity++;
                        }
                        else
                        {
                            reversed3[len3++] = Markups.MISMATCH;
                        }
                        break;

                    case Directions.LEFT:
                        for (int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++)
                        {
                            reversed1[len1++] = PairwiseAlignment.GAP;
                            reversed2[len2++] = s2[--j];
                            reversed3[len3++] = Markups.GAP;
                            gaps++;
                        }
                        break;

                    case Directions.STOP:
                        stillGoing = false;
                        break;
                }
            }


            alignment.Matrix = _scoringMatrix.Name;
            alignment.Gaps = gaps;
            alignment.GapOpenPenalty = _gapOpenPenalty;
            alignment.GapExtendPenalty = _gapExtensionPenalty;
            alignment.Score = cell.Score;
            alignment.Sequence1 = Reverse(reversed1, len1);
            alignment.MarkupLine = Reverse(reversed3, len3);
            alignment.Sequence2 = Reverse(reversed2, len2);
            alignment.Start1 = i;
            alignment.Start2 = j;
            alignment.Identity = identity;
            alignment.Similarity = similarity;

            return alignment;
        }

        /// <summary> 
        /// Returns the maximum of 4 float numbers.
        /// </summary>
        private static float Max(float a, float b, float c, float d)
        {
            return Math.Max(Math.Max(a, b), Math.Max(c, d));
        }

        /// <summary> 
        /// Reverses an array of chars
        /// </summary>
        private static char[] Reverse(char[] a, int len)
        {
            // TODO: replace this method by System.Array.Reverse
            var b = new char[len];

            for (int i = len - 1, j = 0; i >= 0; i--, j++)
            {
                b[j] = a[i];
            }
            return b;
        }
    }
}