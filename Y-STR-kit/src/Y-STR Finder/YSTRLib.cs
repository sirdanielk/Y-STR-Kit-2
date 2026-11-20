
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Y_STR_Finder
{
    public class YSTRLib
    {
        static string _input_vcf_gz = null;
        static string bcftools_exe = null;
        static YSTRLib _lib = null;

        private const char GREEN_MISSING = 'X';
        private const char RED_DELETION = 'R';

        private const string BLUE_STYLE =
            "<span style='border-color:lightblue;border-style: solid; border-width: 1px;background: -webkit-linear-gradient(left top, white , lightblue); background: -o-linear-gradient(bottom right, white, lightblue); background: -moz-linear-gradient(bottom right, white, lightblue); background: linear-gradient(to bottom right, white , lightblue);color:blue'>";

        private const string RED_STYLE =
            "<span style='border-color:lightcoral;border-style: solid; border-width: 1px;background: -webkit-linear-gradient(left top, white , lightcoral); background: -o-linear-gradient(bottom right, white, lightcoral); background: -moz-linear-gradient(bottom right, white, lightcoral); background: linear-gradient(to bottom right, white , lightcoral);color:red'>";

        private const string GREEN_STYLE =
            "<span style='border-color:lightgreen;border-style: solid; border-width: 1px;background: -webkit-linear-gradient(left top, white , lightgreen); background: -o-linear-gradient(bottom right, white, lightgreen); background: -moz-linear-gradient(bottom right, white, lightgreen); background: linear-gradient(to bottom right, white , lightgreen);color:green'>";

        private const string GREY_STYLE =
            "<span style='color:#a0a0a0;'>";

        public static YSTRLib getYSTRLib()
        {
            return _lib;
        }

        public YSTRLib(string input_vcf_gz)
        {
            _input_vcf_gz = input_vcf_gz;
            string exe_path = @"bin\ubin\";
            bcftools_exe = exe_path + "bcftools.exe";
            _lib = this;
        }

        private static string MakeUnknownTolerantRegex(string pattern)
        {
            StringBuilder result = new StringBuilder();

            for (int i = 0; i < pattern.Length; i++)
            {
                char c = pattern[i];

                if (c == '\\' && i + 1 < pattern.Length)
                {
                    result.Append(c);
                    result.Append(pattern[++i]);
                    continue;
                }

                if (c == '[')
                {
                    int end = pattern.IndexOf(']', i + 1);

                    if (end >= 0)
                    {
                        string characterClass =
                            pattern.Substring(i, end - i + 1);

                        result.Append(
                            "(?:" + characterClass +
                            "|(?<!X)X(?!X))");

                        i = end;
                        continue;
                    }
                }

                if (c == 'A' || c == 'C' ||
                    c == 'G' || c == 'T')
                {
                    result.Append(
                        "(?:" + c +
                        "|(?<!X)X(?!X))");
                }
                else
                {
                    result.Append(c);
                }
            }

            return result.ToString();
        }

        private static void MarkDeletion(
            long pos,
            string reference,
            string alternate,
            HashSet<long> deletionPositions)
        {
            if (string.IsNullOrEmpty(reference) ||
                string.IsNullOrEmpty(alternate))
                return;

            int prefix = 0;

            while (prefix < reference.Length &&
                   prefix < alternate.Length &&
                   reference[prefix] == alternate[prefix])
                prefix++;

            int refEnd = reference.Length;
            int altEnd = alternate.Length;

            while (refEnd > prefix &&
                   altEnd > prefix &&
                   reference[refEnd - 1] ==
                   alternate[altEnd - 1])
            {
                refEnd--;
                altEnd--;
            }

            int deletedLength =
                refEnd - prefix - (altEnd - prefix);

            if (deletedLength <= 0 || altEnd != prefix)
                return;

            for (int i = prefix; i < refEnd; i++)
                deletionPositions.Add(pos + i);
        }

        private static void MapAllele(
            long pos,
            string reference,
            string allele,
            Dictionary<long, char> sequence,
            HashSet<long> deletionPositions,
            Dictionary<long, string> insertionsAfter)
        {
            if (string.IsNullOrEmpty(reference) ||
                string.IsNullOrEmpty(allele))
                return;

            int prefix = 0;

            while (prefix < reference.Length &&
                   prefix < allele.Length &&
                   reference[prefix] == allele[prefix])
                prefix++;

            int refEnd = reference.Length;
            int altEnd = allele.Length;

            while (refEnd > prefix &&
                   altEnd > prefix &&
                   reference[refEnd - 1] ==
                   allele[altEnd - 1])
            {
                refEnd--;
                altEnd--;
            }

            for (int i = 0; i < prefix; i++)
                sequence[pos + i] = allele[i];

            int refMiddle = refEnd - prefix;
            int altMiddle = altEnd - prefix;
            int common = Math.Min(refMiddle, altMiddle);

            if (refMiddle == 0 && altMiddle > 0)
            {
                long anchorPosition = pos + prefix - 1;

                string inserted =
                    allele.Substring(prefix, altMiddle);

                if (insertionsAfter.ContainsKey(anchorPosition))
                    insertionsAfter[anchorPosition] += inserted;
                else
                    insertionsAfter.Add(anchorPosition, inserted);
            }

            for (int i = 0; i < common; i++)
                sequence[pos + prefix + i] =
                    allele[prefix + i];

            if (refMiddle > altMiddle)
            {
                for (int i = prefix + altMiddle;
                     i < refEnd; i++)
                    deletionPositions.Add(pos + i);
            }

            int suffixLength = reference.Length - refEnd;

            for (int i = 0; i < suffixLength; i++)
                sequence[pos + refEnd + i] =
                    allele[altEnd + i];
        }
		
        private static int CountAnchoredMotifs(
            string sequence,
            Regex motifRegex,
            List<int> selectedLogicalIndexes,
            Dictionary<int, int> blueLogicalClusterIds)
        {
            int count = 0;
            int position = 0;

            while (position < sequence.Length)
            {
                Match match =
                    motifRegex.Match(sequence, position);

                if (!match.Success ||
                    match.Index != position ||
                    match.Length == 0)
                {
                    position++;
                    continue;
                }

                bool hasRealNucleotide = false;

                for (int i = 0; i < match.Length; i++)
                {
                    char c = sequence[match.Index + i];

                    if (c == 'A' || c == 'C' ||
                        c == 'G' || c == 'T')
                    {
                        hasRealNucleotide = true;
                        break;
                    }
                }

                if (!hasRealNucleotide)
                {
                    position++;
                    continue;
                }

                int clusterId = count;
                count++;

                for (int i = 0; i < match.Length; i++)
                {
                    int index = match.Index + i;

                    if (index < selectedLogicalIndexes.Count &&
                        sequence[index] != GREEN_MISSING)
                    {
                        int logicalIndex =
                            selectedLogicalIndexes[index];

                        blueLogicalClusterIds[logicalIndex] =
                            clusterId;
                    }
                }

                position += match.Length;
            }

            return count;
        }

        public YStrResult getYSTR(
            long pos_start,
            long pos_end,
            string motif_regex,
            string filter_regex,
            string replace_regex)
        {
            Process p = new Process();

            p.StartInfo.FileName = bcftools_exe;

            p.StartInfo.Arguments =
                "query -f '%POS\\t%REF\\t%ALT\\t[%TGT]\\n' -r chrY:"
                + pos_start + "-" + pos_end
                + " \"" + _input_vcf_gz + "\"";

            p.StartInfo.WorkingDirectory = ".";
            p.StartInfo.CreateNoWindow = true;
            p.StartInfo.RedirectStandardOutput = true;
            p.StartInfo.UseShellExecute = false;

            p.StartInfo.EnvironmentVariables["CYGWIN"] =
                "nodosfilewarning";

            p.Start();

            string output = p.StandardOutput.ReadToEnd();
            p.WaitForExit();

            if (p.ExitCode != 0)
                throw new Exception("bcftools error: " + p.ExitCode);

            string tmpFile = Path.Combine(
                Path.GetTempPath(), "ystr.txt");

            File.WriteAllText(tmpFile, output);

            return getYSTR(
                tmpFile,
                motif_regex,
                filter_regex,
                replace_regex,
                pos_start,
                pos_end);
        }

        private YStrResult getYSTR(
            string bcfout_file,
            string motif_regex,
            string filter_regex,
            string replace_regex,
            long pos_start,
            long pos_end)
        {
            YStrResult result = new YStrResult();

            try
            {
                Dictionary<long, char> sequence =
                    new Dictionary<long, char>();

                Dictionary<long, string> insertionsAfter =
                    new Dictionary<long, string>();

                HashSet<long> observedPositions =
                    new HashSet<long>();

                HashSet<long> deletionPositions =
                    new HashSet<long>();

                long minPos = long.MaxValue;
                long maxPos = long.MinValue;

                foreach (string line in File.ReadLines(bcfout_file))
                {
                    if (string.IsNullOrWhiteSpace(line))
                        continue;

                    string[] data = line.Split('\t');

                    if (data.Length < 4)
                        continue;

                    long pos = long.Parse(data[0]);

                    string reference = data[1];
                    string alternate = data[2];
                    string genotype = data[3];

                    if (string.IsNullOrEmpty(reference) ||
                        reference == ".")
                        continue;

                    string[] alleles = genotype.Split(
                        new char[] { '/', '|' });

                    string allele = alleles[0].Trim();

                    if (string.IsNullOrEmpty(allele) ||
                        allele == ".")
                        allele = reference;

                    minPos = Math.Min(minPos, pos);

                    maxPos = Math.Max(
                        maxPos,
                        pos + reference.Length - 1);

                    for (int i = 0; i < reference.Length; i++)
                        observedPositions.Add(pos + i);

                    MarkDeletion(
                        pos,
                        reference,
                        allele,
                        deletionPositions);

                    MapAllele(
                        pos,
                        reference,
                        allele,
                        sequence,
                        deletionPositions,
                        insertionsAfter);
                }

                if (minPos == long.MaxValue)
                {
                    result.value = 0;
                    result.reliability = 0;
                    result.raw_seq = "";
                    result.html_seq = "";
                    return result;
                }

                StringBuilder logicalBuilder =
                    new StringBuilder();

                StringBuilder rawBuilder =
                    new StringBuilder();

                List<int> logicalCharacterIndexes =
                    new List<int>();

                for (long i = pos_start; i <= pos_end; i++)
                {
                    if (deletionPositions.Contains(i))
                    {
                        rawBuilder.Append('N');
                        continue;
                    }

                    if (sequence.ContainsKey(i))
                    {
                        char nucleotide = sequence[i];

                        rawBuilder.Append(nucleotide);

                        logicalCharacterIndexes.Add(
                            logicalBuilder.Length);

                        logicalBuilder.Append(nucleotide);

                        if (insertionsAfter.ContainsKey(i))
                        {
                            string inserted = insertionsAfter[i];

                            rawBuilder.Append(inserted);

                            for (int j = 0; j < inserted.Length; j++)
                            {
                                logicalCharacterIndexes.Add(
                                    logicalBuilder.Length);

                                logicalBuilder.Append(inserted[j]);
                            }
                        }
                    }
                    else if (i >= minPos && i <= maxPos)
                    {
                        rawBuilder.Append('N');

                        logicalCharacterIndexes.Add(
                            logicalBuilder.Length);

                        logicalBuilder.Append(GREEN_MISSING);
                    }
                    else
                    {
                        rawBuilder.Append('N');
                    }
                }

                string logicalSeq = logicalBuilder.ToString();
                result.raw_seq = rawBuilder.ToString();

                long totalPositions =
                    pos_end - pos_start + 1;

                if (totalPositions > 0)
                {
                    result.reliability = (int)(
                        observedPositions.Count * 100 /
                        totalPositions);

                    if (result.reliability > 100)
                        result.reliability = 100;
                }
                else
                {
                    result.reliability = 0;
                }

                string tolerantMotif =
                    MakeUnknownTolerantRegex(motif_regex);

                string tolerantFilter =
                    MakeUnknownTolerantRegex(filter_regex);

                Regex motifRegex = new Regex(
                    tolerantMotif,
                    RegexOptions.IgnoreCase);

                Regex filterRegex = new Regex(
                    tolerantFilter,
                    RegexOptions.IgnoreCase);

                StringBuilder selectedBuilder =
                    new StringBuilder();

                List<int> selectedLogicalIndexes =
                    new List<int>();

                if (replace_regex == "CONCAT")
                {
                    MatchCollection filterMatches =
                        filterRegex.Matches(logicalSeq);

                    foreach (Match match in filterMatches)
                    {
                        selectedBuilder.Append(match.Value);

                        for (int i = match.Index;
                             i < match.Index + match.Length;
                             i++)
                        {
                            selectedLogicalIndexes.Add(
                                logicalCharacterIndexes[i]);
                        }
                    }
                }
                else
                {
                    selectedBuilder.Append(logicalSeq);

                    selectedLogicalIndexes.AddRange(
                        logicalCharacterIndexes);
                }

                string selectedSeq = selectedBuilder.ToString();

                Dictionary<int, int> blueLogicalClusterIds =
                    new Dictionary<int, int>();

                result.value = CountAnchoredMotifs(
                    selectedSeq,
                    motifRegex,
                    selectedLogicalIndexes,
                    blueLogicalClusterIds);

                StringBuilder html = new StringBuilder();

                int logicalIndex = 0;
                int? activeBlueClusterId = null;

                for (long i = pos_start; i <= pos_end; i++)
                {

                    if (deletionPositions.Contains(i))
                    {
                        if (activeBlueClusterId.HasValue)
                        {
                            html.Append("</span>");
                            activeBlueClusterId = null;
                        }

                        html.Append(RED_STYLE);
                        html.Append("&middot;</span>");
                        continue;
                    }

                    if (sequence.ContainsKey(i))
                    {
                        char nucleotide = sequence[i];

                        int clusterId;

                        if (blueLogicalClusterIds.TryGetValue(
                            logicalIndex, out clusterId))
                        {
                            if (activeBlueClusterId != clusterId)
                            {
                                if (activeBlueClusterId.HasValue)
                                    html.Append("</span>");

                                html.Append(BLUE_STYLE);
                                activeBlueClusterId = clusterId;
                            }

                            html.Append(nucleotide);
                        }
                        else
                        {
                            if (activeBlueClusterId.HasValue)
                            {
                                html.Append("</span>");
                                activeBlueClusterId = null;
                            }

                            html.Append(nucleotide);
                        }

                        logicalIndex++;

                        if (insertionsAfter.ContainsKey(i))
                        {
                            string inserted = insertionsAfter[i];

                            for (int j = 0; j < inserted.Length; j++)
                            {
                                char insertedNucleotide = inserted[j];

                                if (blueLogicalClusterIds.TryGetValue(
                                    logicalIndex, out clusterId))
                                {
                                    if (activeBlueClusterId != clusterId)
                                    {
                                        if (activeBlueClusterId.HasValue)
                                            html.Append("</span>");

                                        html.Append(BLUE_STYLE);
                                        activeBlueClusterId = clusterId;
                                    }

                                    html.Append(insertedNucleotide);
                                }
                                else
                                {
                                    if (activeBlueClusterId.HasValue)
                                    {
                                        html.Append("</span>");
                                        activeBlueClusterId = null;
                                    }

                                    html.Append(insertedNucleotide);
                                }

                                logicalIndex++;
                            }
                        }

                        continue;
                    }

                    if (i >= minPos && i <= maxPos)
                    {
                        if (activeBlueClusterId.HasValue)
                        {
                            html.Append("</span>");
                            activeBlueClusterId = null;
                        }

                        html.Append(GREEN_STYLE);
                        html.Append("&middot;</span>");

                        logicalIndex++;
                    }
                    else
                    {

                        if (activeBlueClusterId.HasValue)
                        {
                            html.Append("</span>");
                            activeBlueClusterId = null;
                        }

                        html.Append(GREY_STYLE);
                        html.Append("&middot;</span>");
                    }
                }

                if (activeBlueClusterId.HasValue)
                    html.Append("</span>");

                result.html_seq = html.ToString();
            }
            catch (Exception ex)
            {
                result.value = 0;
                Debug.WriteLine("YSTR ERROR: " + ex.ToString());
            }

            return result;
        }
    }
}
