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
    /// <summary> JAligner exception </summary>
    [Serializable]
    public class SmithWatermanGotohException : Exception
    {
        public SmithWatermanGotohException()
        {
        }

        public SmithWatermanGotohException(string message) : base(message)
        {
        }

        public SmithWatermanGotohException(string message, Exception cause) : base(message, cause)
        {
        }

        public SmithWatermanGotohException(Exception cause) : base("", cause)
        {
        }
    }
}