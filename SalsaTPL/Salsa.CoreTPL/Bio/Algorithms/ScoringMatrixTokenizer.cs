using System;
using System.Collections;

namespace Salsa.Core.Bio.Algorithms
{
    /// <summary>The class performs token processing in strings</summary>
    /// <remarks>Class create by the Java Language Conversion Assistant.</remarks>
    internal class ScoringMatrixTokenizer : IEnumerator
    {
        /// Char representation of the String to tokenize.
        private readonly char[] chars;

        /// Include demiliters in the results.
        private readonly bool includeDelims;

        /// Position over the string
        private long currentPos;

        //The tokenizer uses the default delimiter set: the space character, the tab character, the newline character, and the carriage-return character and the form-feed character
        private string delimiters = " \t\n\r\f";

        /// <summary>
        /// Initializes a new class instance with a specified string to process
        /// </summary>
        /// <param name="source">String to tokenize</param>
        public ScoringMatrixTokenizer(string source)
        {
            chars = source.ToCharArray();
        }

        /// <summary>
        /// Initializes a new class instance with a specified string to process
        /// and the specified token delimiters to use
        /// </summary>
        /// <param name="source">String to tokenize</param>
        /// <param name="delimiters">String containing the delimiters</param>
        public ScoringMatrixTokenizer(string source, string delimiters) : this(source)
        {
            this.delimiters = delimiters;
        }


        /// <summary>
        /// Initializes a new class instance with a specified string to process, the specified token 
        /// delimiters to use, and whether the delimiters must be included in the results.
        /// </summary>
        /// <param name="source">String to tokenize</param>
        /// <param name="delimiters">String containing the delimiters</param>
        /// <param name="includeDelims">Determines if delimiters are included in the results.</param>
        public ScoringMatrixTokenizer(string source, string delimiters, bool includeDelims) : this(source, delimiters)
        {
            this.includeDelims = includeDelims;
        }


        /// <summary>
        /// Remaining tokens count
        /// </summary>
        public int Count
        {
            get
            {
                //keeping the current pos
                long pos = currentPos;
                int i = 0;

                try
                {
                    while (true)
                    {
                        NextToken();
                        i++;
                    }
                }
                catch (ArgumentOutOfRangeException)
                {
                    currentPos = pos;
                    return i;
                }
            }
        }

        #region IEnumerator Members

        /// <summary>
        ///  Performs the same action as NextToken.
        /// </summary>
        public Object Current
        {
            get { return NextToken(); }
        }

        /// <summary>
        ///  Performs the same action as HasMoreTokens.
        /// </summary>
        /// <returns>True or false, depending if there are more tokens</returns>
        public bool MoveNext()
        {
            return HasMoreTokens();
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void Reset()
        {
            ;
        }

        #endregion

        /// <summary>
        /// Returns the next token from the token list
        /// </summary>
        /// <returns>The string value of the token</returns>
        public string NextToken()
        {
            return NextToken(delimiters);
        }

        /// <summary>
        /// Returns the next token from the source string, using the provided
        /// token delimiters
        /// </summary>
        /// <param name="delimiters">String containing the delimiters to use</param>
        /// <returns>The string value of the token</returns>
        public string NextToken(string delimiters)
        {
            //According to documentation, the usage of the received delimiters should be temporary (only for this call).
            //However, it seems it is not true, so the following line is necessary.
            this.delimiters = delimiters;

            //at the end 
            if (currentPos == chars.Length)
                throw new ArgumentOutOfRangeException();
                //if over a delimiter and delimiters must be returned
            else if ((Array.IndexOf(delimiters.ToCharArray(), chars[currentPos]) != -1)
                     && includeDelims)
                return "" + chars[currentPos++];
                //need to get the token wo delimiters.
            else
                return nextToken(delimiters.ToCharArray());
        }

        //Returns the nextToken wo delimiters
        private string nextToken(char[] delimiters)
        {
            string token = "";
            long pos = currentPos;

            //skip possible delimiters
            while (Array.IndexOf(delimiters, chars[currentPos]) != -1)
                //The last one is a delimiter (i.e there is no more tokens)
                if (++currentPos == chars.Length)
                {
                    currentPos = pos;
                    throw new ArgumentOutOfRangeException();
                }

            //getting the token
            while (Array.IndexOf(delimiters, chars[currentPos]) == -1)
            {
                token += chars[currentPos];
                //the last one is not a delimiter
                if (++currentPos == chars.Length)
                    break;
            }
            return token;
        }


        /// <summary>
        /// Determines if there are more tokens to return from the source string
        /// </summary>
        /// <returns>True or false, depending if there are more tokens</returns>
        public bool HasMoreTokens()
        {
            //keeping the current pos
            long pos = currentPos;

            try
            {
                NextToken();
            }
            catch (ArgumentOutOfRangeException)
            {
                return false;
            }
            finally
            {
                currentPos = pos;
            }
            return true;
        }
    }
}