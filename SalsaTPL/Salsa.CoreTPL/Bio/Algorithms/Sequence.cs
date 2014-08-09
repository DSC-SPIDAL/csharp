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
using Salsa.Core.Bio.IO;

namespace Salsa.Core.Bio.Algorithms
{
    [Serializable]
    public class Sequence : IFastaRecord
    {
        #region Members

        private string _label = string.Empty;
        private string _residues = string.Empty;

        #endregion

        /// <summary> Constructor</summary>
        public Sequence()
        {
        }

        /// <summary> Constructor</summary>
        public Sequence(string label, string sequence)
        {
            _label = label;
            _residues = sequence;
        }

        public char this[int residueIndex]
        {
            get { return _residues[residueIndex]; }
        }

        public int Length
        {
            get { return _residues.Length; }
        }

        #region IFastaRecord Members

        /// <summary> Returns the sequence label</summary>
        public string Label
        {
            get { return _label; }
            set { _label = value; }
        }

        /// <summary> Gets or sets the amino acid sequence</summary>
        public string Residues
        {
            get { return _residues; }
            set { _residues = value; }
        }

        #endregion

        /// <summary> Returns a subsequence</summary>
        /// <param name="index">start index </param>
        /// <param name="length">length of subsequence </param>
        /// <returns> subsequence </returns>
        public string Subsequence(int index, int length)
        {
            return _residues.Substring(index, (index + length) - (index));
        }
    }
}