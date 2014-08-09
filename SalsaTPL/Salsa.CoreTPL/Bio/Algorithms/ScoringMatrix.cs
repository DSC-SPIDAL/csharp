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
using System.IO;
using System.Reflection;

namespace Salsa.Core.Bio.Algorithms
{
    [Serializable]
    public class ScoringMatrix
    {
        private const char COMMENT_STARTER = '#';
        private const int SIZE = 127;
        private static readonly string[] _matrixNames;
        private readonly string _matrixName;
        private readonly float[,] _scores;

        static ScoringMatrix()
        {
            _matrixNames = new[]
                               {
                                   "BLOSUM100",
                                   "BLOSUM30",
                                   "BLOSUM35",
                                   "BLOSUM40",
                                   "BLOSUM45",
                                   "BLOSUM50",
                                   "BLOSUM55",
                                   "BLOSUM60",
                                   "BLOSUM62",
                                   "BLOSUM65",
                                   "BLOSUM70",
                                   "BLOSUM75",
                                   "BLOSUM80",
                                   "BLOSUM85",
                                   "BLOSUM90",
                                   "BLOSUMN",
                                   "DAYHOFF",
                                   "EDNAFULL",
                                   "GONNET",
                                   "IDENTITY",
                                   "MATCH",
                                   "NUC44",
                                   "PAM10",
                                   "PAM100",
                                   "PAM110",
                                   "PAM120",
                                   "PAM130",
                                   "PAM140",
                                   "PAM150",
                                   "PAM160",
                                   "PAM170",
                                   "PAM180",
                                   "PAM190",
                                   "PAM20",
                                   "PAM200",
                                   "PAM210",
                                   "PAM220",
                                   "PAM230",
                                   "PAM240",
                                   "PAM250",
                                   "PAM260",
                                   "PAM270",
                                   "PAM280",
                                   "PAM290",
                                   "PAM30",
                                   "PAM300",
                                   "PAM310",
                                   "PAM320",
                                   "PAM330",
                                   "PAM340",
                                   "PAM350",
                                   "PAM360",
                                   "PAM370",
                                   "PAM380",
                                   "PAM390",
                                   "PAM40",
                                   "PAM400",
                                   "PAM410",
                                   "PAM420",
                                   "PAM430",
                                   "PAM440",
                                   "PAM450",
                                   "PAM460",
                                   "PAM470",
                                   "PAM480",
                                   "PAM490",
                                   "PAM50",
                                   "PAM500",
                                   "PAM60",
                                   "PAM70",
                                   "PAM80",
                                   "PAM90",
                                   "TEST1"
                               };
        }

        public ScoringMatrix(string matrixName, float[,] scores)
        {
            _matrixName = matrixName;
            _scores = scores;
        }

        public string Name
        {
            get { return _matrixName; }
        }

        public float[,] Scores
        {
            get { return _scores; }
        }

        public static string[] MatrixNames
        {
            get { return _matrixNames; }
        }

        public float GetScore(char a, char b)
        {
            return _scores[a, b];
        }

        public static ScoringMatrix Load(string matrixName)
        {
            var acids = new char[SIZE];

            // Initialize the acids array to null values (ascii = 0)
            for (int i = 0; i < SIZE; i++)
            {
                acids[i] = (char) (0);
            }

            var scores = new float[SIZE,SIZE];

            using (
                Stream stream =
                    Assembly.GetExecutingAssembly().GetManifestResourceStream("Salsa.Core.Bio.Algorithms.Matrices." +
                                                                              matrixName))
            {
                using (var reader = new StreamReader(stream))
                {
                    // Skip the comment lines
                    string line;

                    while ((line = reader.ReadLine()) != null && line.Trim()[0] == COMMENT_STARTER)
                        ;

                    // Read the headers line (the letters of the acids)
                    var tokenizer = new ScoringMatrixTokenizer(line.Trim());

                    for (int j = 0; tokenizer.HasMoreTokens(); j++)
                    {
                        acids[j] = tokenizer.NextToken()[0];
                    }

                    // Read the scores
                    while ((line = reader.ReadLine()) != null)
                    {
                        tokenizer = new ScoringMatrixTokenizer(line.Trim());
                        char acid = tokenizer.NextToken()[0];

                        for (int i = 0; i < SIZE; i++)
                        {
                            if (acids[i] != 0)
                            {
                                scores[acid, acids[i]] = Single.Parse(tokenizer.NextToken());
                            }
                        }
                    }
                }
            }

            return new ScoringMatrix(matrixName, scores);
        }
    }
}