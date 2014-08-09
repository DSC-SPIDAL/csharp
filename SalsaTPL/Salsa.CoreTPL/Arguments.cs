/*
* Arguments class: application arguments interpreter
*
* Authors:		R. LOPES
* Contributors:	R. LOPES
* Created:		25 October 2002
* Modified:		28 October 2002
*
* Version:		1.0
 * 
 * http://www.codeproject.com/KB/recipes/command_line.aspx
*/

using System;
using System.Collections.Specialized;
using System.Text.RegularExpressions;

namespace Salsa.Core
{
    /// <summary>
    /// Arguments class
    /// </summary>
    /// Valid parameters forms:
    /// {-,/,--}param{ ,=,:}((",')value(",'))
    /// Examples: -param1 value1 --param2 /param3:"Test-:-work" /param4=happy -param5 '--=nice=--'
    public class Arguments
    {
        // Variables
        private readonly StringDictionary _parameters;
        private string _usage = "Usage:";

        // Constructor
        public Arguments(string[] Args)
        {
            _parameters = new StringDictionary();
            var Spliter = new Regex(@"^-{1,2}|^/|=|:", RegexOptions.IgnoreCase | RegexOptions.Compiled);
            var Remover = new Regex(@"^['""]?(.*?)['""]?$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
            string parameter = null;
            string[] parts;

            // Valid parameters forms:
            // {-,/,--}param{ ,=,:}((",')value(",'))
            // Examples: -param1 value1 --param2 /param3:"Test-:-work" /param4=happy -param5 '--=nice=--'
            foreach (string argument in Args)
            {
                // Look for new parameters (-,/ or --) and a possible enclosed value (=,:)
                parts = Spliter.Split(argument, 3);

                switch (parts.Length)
                {
                        // Found a value (for the last parameter found (space separator))
                    case 1:
                        if (parameter != null)
                        {
                            if (!_parameters.ContainsKey(parameter))
                            {
                                parts[0] = Remover.Replace(parts[0], "$1");
                                _parameters.Add(parameter, parts[0]);
                            }
                            parameter = null;
                        }
                        // else Error: no parameter waiting for a value (skipped)
                        break;
                        // Found just a parameter
                    case 2:
                        // The last parameter is still waiting. With no value, set it to true.
                        if (parameter != null)
                        {
                            if (!_parameters.ContainsKey(parameter))
                                _parameters.Add(parameter, "true");
                        }
                        parameter = parts[1];
                        break;
                        // Parameter with enclosed value
                    case 3:
                        // The last parameter is still waiting. With no value, set it to true.
                        if (parameter != null)
                        {
                            if (!_parameters.ContainsKey(parameter))
                                _parameters.Add(parameter, "true");
                        }
                        parameter = parts[1];
                        // Remove possible enclosing characters (",')
                        if (!_parameters.ContainsKey(parameter))
                        {
                            parts[2] = Remover.Replace(parts[2], "$1");
                            _parameters.Add(parameter, parts[2]);
                        }
                        parameter = null;
                        break;
                }
            }
            // In case a parameter is still waiting
            if (parameter != null)
            {
                if (!_parameters.ContainsKey(parameter))
                    _parameters.Add(parameter, "true");
            }
        }

        public string this[string parameterName]
        {
            get { return (_parameters[parameterName]); }
        }

        public string Usage
        {
            get { return _usage; }
            set { _usage = value; }
        }

        public bool Contains(string parameterName)
        {
            return _parameters.ContainsKey(parameterName);
        }

        public T GetValue<T>(string parameterName)
        {
            return (T) Convert.ChangeType(_parameters[parameterName], typeof (T));
        }

        public T GetValue<T>(string parameterName, T defaultValue)
        {
            if (Contains(parameterName) == false)
            {
                return defaultValue;
            }
            else
            {
                return (T) Convert.ChangeType(_parameters[parameterName], typeof (T));
            }
        }

        public bool CheckRequired(string[] requiredParameters)
        {
            foreach (string parameter in requiredParameters)
            {
                if (_parameters.ContainsKey(parameter) == false)
                {
                    return false;
                }
            }

            return true;
        }
    }
}