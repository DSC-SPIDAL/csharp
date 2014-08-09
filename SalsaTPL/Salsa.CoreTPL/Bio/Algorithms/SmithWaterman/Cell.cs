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

namespace Salsa.Core.Bio.Algorithms.SmithWaterman
{
    /// <summary> 
    /// A cell in a similarity matrix, to hold row, column and score.
    /// </summary>
    internal sealed class Cell
    {
        #region Members

        private int column;
        private int row;
        private float score;

        #endregion

        public Cell()
        {
            row = 0;
            column = 0;
            score = float.NegativeInfinity;
        }

        public int Column
        {
            get { return column; }
            set { column = value; }
        }

        public int Row
        {
            get { return row; }
            set { row = value; }
        }

        public float Score
        {
            get { return score; }
            set { score = value; }
        }

        public void Set(int row, int column, float score)
        {
            this.row = row;
            this.column = column;
            this.score = score;
        }
    }
}