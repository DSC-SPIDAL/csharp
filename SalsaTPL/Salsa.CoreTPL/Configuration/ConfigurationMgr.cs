using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;
using Salsa.Core.Configuration.Sections;

namespace Salsa.Core.Configuration
{
    [XmlRoot(ElementName = "Configuration")]
    public class ConfigurationMgr
    {
        #region Constructors

        public ConfigurationMgr()
        {
            _global = new GlobalSection();
            _manxcat = new ManxcatSection();
            _pairwise = new PairwiseSection();
            _smithWaterman = new SmithWatermanSection();
            _needlemanWunsch = new NeedlemanWunschSection();
            _smithWatermanMS = new SmithWatermanMS();
            _daVectorSpongeSection = new DAVectorSpongeSection();
        }

        #endregion

        #region Members

        private DAVectorSpongeSection _daVectorSpongeSection;
        private GlobalSection _global;
        private ManxcatSection _manxcat;
        private NeedlemanWunschSection _needlemanWunsch;
        private PairwiseSection _pairwise;
        private SmithWatermanSection _smithWaterman;
        private SmithWatermanMS _smithWatermanMS;

        #endregion

        #region Properties

        public GlobalSection GlobalSection
        {
            get { return _global; }
            set { _global = value; }
        }

        public ManxcatSection ManxcatSection
        {
            get { return _manxcat; }
            set { _manxcat = value; }
        }

        public PairwiseSection PairwiseSection
        {
            get { return _pairwise; }
            set { _pairwise = value; }
        }

        public SmithWatermanSection SmithWatermanSection
        {
            get { return _smithWaterman; }
            set { _smithWaterman = value; }
        }

        public SmithWatermanMS SmithWatermanMS
        {
            get { return _smithWatermanMS; }
            set { _smithWatermanMS = value; }
        }

        public NeedlemanWunschSection NeedlemanWunschSection
        {
            get { return _needlemanWunsch; }
            set { _needlemanWunsch = value; }
        }

        public DAVectorSpongeSection DAVectorSpongeSection
        {
            get { return _daVectorSpongeSection; }
            set { _daVectorSpongeSection = value; }
        }

        #endregion

        public string[] GetSectionNames()
        {
            return new[]
                       {
                           "Global", "Manxcat", "Pairwise", "SmithWaterman", "SmithWatermanMS", "NeedlemanWunsch",
                           "DAVectorSponge"
                       };
        }

        public Section GetSection(string sectionName)
        {
            Section section = null;

            switch (sectionName.ToUpper())
            {
                case "GLOBAL":
                    section = _global;
                    break;
                case "PAIRWISE":
                case "PWC":
                    section = _pairwise;
                    break;
                case "MANXCAT":
                case "MDS":
                    section = _manxcat;
                    break;
                case "SMITHWATERMAN":
                case "SWG":
                    section = _smithWaterman;
                    break;
                case "SMITHWATERMANMS":
                case "SWMS":
                    section = _smithWatermanMS;
                    break;
                case "NEEDLEMANWUNSCH":
                case "NW":
                    section = _needlemanWunsch;
                    break;
                case "SPONGE":
                case "DAVECTORSPONGE":
                    section = _daVectorSpongeSection;
                    break;
            }


            return section;
        }

        public void SaveAs(string fileName)
        {
            using (var writer = new StreamWriter(fileName))
            {
                var serializer = new XmlSerializer(typeof (ConfigurationMgr));
                serializer.Serialize(writer, this);
                writer.Close();
            }
        }

        public static ConfigurationMgr LoadConfiguration(string fileName, bool expandMacros)
        {
            ConfigurationMgr manager = null;

            using (var reader = new StreamReader(fileName))
            {
                var serializer = new XmlSerializer(typeof (ConfigurationMgr));
                manager = serializer.Deserialize(reader) as ConfigurationMgr;
                reader.Close();
            }


            if (expandMacros)
            {
                ExpandMacros(manager);
            }

            return manager;
        }

        private static void ExpandMacros(ConfigurationMgr manager)
        {
            var macros = new List<MacroReplacement>();
            macros.Add(new MacroReplacement("$(ConfigRootPath)", manager.GlobalSection.ConfigRootPath));
            macros.Add(new MacroReplacement("$(InputRootPath)", manager.GlobalSection.InputRootPath));
            macros.Add(new MacroReplacement("$(OutputRootPath)", manager.GlobalSection.OutputRootPath));

            foreach (string section in manager.GetSectionNames())
            {
                foreach (MacroReplacement macro in macros)
                {
                    Section sec = manager.GetSection(section);
                    sec.ExpandMacro(macro);
                    sec.ExpandEnvVars();
                }
            }
        }
    }
}