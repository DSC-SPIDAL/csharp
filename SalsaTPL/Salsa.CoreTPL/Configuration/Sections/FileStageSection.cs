namespace Salsa.Core.Configuration.Sections
{
    public class FileStageSection : Section
    {
    }

    public class FileStageItem
    {
        public string Source { get; set; }

        public string Target { get; set; }
    }
}