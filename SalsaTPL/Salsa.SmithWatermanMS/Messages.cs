using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Salsa.SmithWatermanMS
{
    class Messages
    {
        public static string ParameterConfigFile = @"configFile";
        public static string ParameterNodeCount = @"nodeCount";
        public static string ParameterJobId = @"%CCP_JOBID%";
        public static string ParameterTaskId = @"%CCP_TASKID%";

        public static string ErrorCreatingFile = @"Error: Unable to create file {0} due to {1}";
        public static string ErrorUnableToWriteIamDone = @"Error: Rank {0} unable to write iamdone status due to {1}";
        public static string ErrorUnableToDisposeEnv = @"Error: Rank {0} unable dispose env due to {1}";
        public static string ErrorGeneral = @"Error: Rank {0} got error {1}";

        public static string WarningIgnoreAlignmentWrite =
            @"Warning: Ignoring alignment write because couldn't find file: {0}";

        public static string InfoUsage = @"Usage: Salsa.SmithWatermanMS.exe /configFile=<string> /nodeCount=<int>";
        public static string Restart = @"Restart: Rank {0} restarting as fake rank {1} in fake world size {2} with last saved cb# {3}";
        public static string InfoReadFasta = @"Info: Rank {0} read {1} sequence from Fasta File {2}";
        public static string InfoBeginBlocksCompute = @"Info: Rank {0} beginning blocks compute: {1}";
        public static string InfoEndBlocksCompute = @"Info: Rank {0} finished blocks compute: {1}";
        public static string DoneMsg = @"Done: Rank {0} done.";
        public static string StatusInitiallyComputed = @"    Status: Rank {0} initially computed blocks {1} of {2} - remaining {3}";
        public static string StatusComputingBlocks = @"    Status: Rank {0} {1}% done, computing block {2} of {3} - (rb#: {4}, cb#: {5})";
        public static string StatusComputingBlocksCompleted = @"    Status: Rank {0} completed {1}/{2} blocks in {3} at {4}";
    }
}
