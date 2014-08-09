using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Drawing.Design;
using System.Windows.Forms.Design;

namespace Salsa.Core.Configuration.Sections
{
    public class GlobalSection : Section
    {
        #region Constructor

        #endregion

        #region Members

        private readonly Collection<string> _emailAddresses = new Collection<string>();
        private string _configRootPath = string.Empty;
        private string _inputRootPath = string.Empty;
        private string _outputRootPath = string.Empty;

        #endregion

        #region Properties

        public Collection<string> EmailAddresses
        {
            get { return _emailAddresses; }
        }

        [Category("I/O")]
        [Description("The full path to the config root folder.  Macro: $(ConfigRootPath)")]
        [Editor(typeof (FolderNameEditor), typeof (UITypeEditor))]
        public string ConfigRootPath
        {
            get { return _configRootPath; }
            set
            {
                _configRootPath = value.Trim();
                OnPropertyChanged("ConfigRootPath");
            }
        }


        [Category("I/O")]
        [Description("The full path to the input root folder.   Macro: $(InputRootPath)")]
        [Editor(typeof (FolderNameEditor), typeof (UITypeEditor))]
        public string InputRootPath
        {
            get { return _inputRootPath; }
            set
            {
                _inputRootPath = value.Trim();
                OnPropertyChanged("InputRootPath");
            }
        }


        [Category("I/O")]
        [Description("The full path to the output root folder.  Macro: $(OutputRootPath)")]
        [Editor(typeof (FolderNameEditor), typeof (UITypeEditor))]
        public string OutputRootPath
        {
            get { return _outputRootPath; }
            set
            {
                _outputRootPath = value.Trim();
                OnPropertyChanged("OutputRootPath");
            }
        }

        #endregion
    }
}