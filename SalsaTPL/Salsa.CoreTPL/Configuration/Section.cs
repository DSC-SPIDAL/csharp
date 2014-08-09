using System.ComponentModel;

namespace Salsa.Core.Configuration
{
    [Category("Configuration")]
    [TypeConverter(typeof (ExpandableObjectConverter))]
    public abstract class Section : INotifyPropertyChanged
    {
        #region INotifyPropertyChanged Members

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion

        protected void OnPropertyChanged(string propertyName)
        {
            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs(propertyName));
            }
        }

        internal virtual void ExpandMacro(MacroReplacement macroReplacement)
        {
        }

        internal virtual void ExpandEnvVars()
        {
        }
    }
}