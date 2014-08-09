namespace Salsa.Core.Configuration
{
    public class MacroReplacement
    {
        private readonly string _macroTemplate;
        private readonly string _replacementValue;

        public MacroReplacement(string macroTemplate, string replacementValue)
        {
            _macroTemplate = macroTemplate;
            _replacementValue = replacementValue;
        }

        public string MacroTemplate
        {
            get { return _macroTemplate; }
        }

        public string ReplacementValue
        {
            get { return _replacementValue; }
        }

        public string Expand(string target)
        {
            return target.Replace(_macroTemplate, _replacementValue);
        }
    }
}