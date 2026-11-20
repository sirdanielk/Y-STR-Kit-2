using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Y_STR_Finder
{
    public class YSTRLib
    {
        private readonly string _inputVcf;
        private readonly string _bcftoolsExe;
        private static YSTRLib _lib = null;

        private const string MotifStyle =
            "border-color:lightblue;" +
            "border-style:solid;" +
            "border-width:1px;" +
            "background:linear-gradient(to bottom right, white, lightblue);" +
            "color:blue";

        private const string DeletionStyle =
            "border-color:lightcoral;" +
            "border-style:solid;" +
            "border-width:1px;" +
            "background:linear-gradient(to bottom right, white, lightcoral);" +
            "color:red";

        private const string GreyStyle =
            "color:#a0a0a0;";

        public static YSTRLib getYSTRLib()
        {
            return _lib;
        }

        public YSTRLib(string input_vcf_gz)
        {
            _inputVcf = Path.GetFullPath(input_vcf_gz);

            _bcftoolsExe = Path.GetFullPath(
                @"bin\ubin\bcftools.exe"
            );

            if (!File.Exists(_inputVcf))
                throw new FileNotFoundException(
                    "Nie znaleziono pliku VCF", _inputVcf
                );

            if (!File.Exists(_bcftoolsExe))
                throw new FileNotFoundException(
                    "Nie znaleziono bcftools.exe", _bcftoolsExe
                );

            _lib = this;
        }

        public YStrResult getYSTR(
            long pos_start,
            long pos_end,
            string motif_regex,
            string filter_regex,
            string replace_regex)
        {
            YStrResult result = new YStrResult();

            try
            {
                if (pos_start <= 0 || pos_end < pos_start)
                    throw new ArgumentException("Nieprawidłowy zakres pozycji.");

                string output = QueryVcf(pos_start, pos_end);

                return ReconstructSequence(
                    output,
                    pos_start,
                    pos_end,
                    motif_regex,
                    filter_regex,
                    replace_regex
                );
            }
            catch (Exception ex)
            {
                Console.WriteLine("YSTR error: " + ex.Message);

                result.value = 0;
                result.reliability = 0;
                result.raw_seq = "";
                result.html_seq = "";

                return result;
            }
        }

        private string QueryVcf(long pos_start, long pos_end)
        {
            ProcessStartInfo psi = new ProcessStartInfo();

            psi.FileName = _bcftoolsExe;

            psi.Arguments =
                "query " +
                "-r " + Quote("chrY:" + pos_start + "-" + pos_end) + " " +
                "-f " + Quote("%POS\\t%REF\\t%ALT\\t[%GT]\\n") + " " +
                Quote(_inputVcf);

            psi.WorkingDirectory = Environment.CurrentDirectory;

            psi.CreateNoWindow = true;
            psi.UseShellExecute = false;
            psi.RedirectStandardOutput = true;
            psi.RedirectStandardError = true;

            psi.EnvironmentVariables["CYGWIN"] = "nodosfilewarning";

            using (Process p = new Process())
            {
                p.StartInfo = psi;

                StringBuilder errorOutput = new StringBuilder();

                p.ErrorDataReceived += delegate (
                    object sender,
                    DataReceivedEventArgs e)
                {
                    if (e.Data != null)
                        errorOutput.AppendLine(e.Data);
                };

                p.Start();
                p.BeginErrorReadLine();

                string output = p.StandardOutput.ReadToEnd();

                p.WaitForExit();

                if (p.ExitCode != 0)
                {
                    throw new Exception(
                        "Błąd bcftools:\r\n" + errorOutput.ToString()
                    );
                }

                return output;
            }
        }

        private string Quote(string value)
        {
            return "\"" + value.Replace("\"", "\\\"") + "\"";
        }

        private class VcfRecord
        {
            public long Position;
            public string Ref;
            public string Alt;
            public string GT;

            public int AlleleIndex;
            public string SelectedAllele;

            public bool IsMissing;
            public bool IsStar;
            public bool IsSymbolic;
        }

        private int ParseGenotype(string gt)
        {
            if (String.IsNullOrWhiteSpace(gt))
                return -1;

            gt = gt.Trim();

            if (gt == "." || gt == "./." || gt == ".|.")
                return -1;

            string[] alleles = gt.Split(
                new char[] { '/', '|' },
                StringSplitOptions.None
            );

            if (alleles.Length == 0)
                return -1;

            if (alleles.Any(x => x == "."))
                return -1;

            int firstAllele;

            if (!Int32.TryParse(alleles[0], out firstAllele))
                return -1;

            return firstAllele;
        }

        private VcfRecord ParseRecord(string line)
        {
            string[] fields = line.Split('\t');

            if (fields.Length < 4)
                return null;

            VcfRecord record = new VcfRecord();

            record.Position = Int64.Parse(
                fields[0],
                CultureInfo.InvariantCulture
            );

            record.Ref = fields[1].Trim().ToUpperInvariant();
            record.Alt = fields[2].Trim().ToUpperInvariant();
            record.GT = fields[3].Trim();

            record.AlleleIndex = ParseGenotype(record.GT);

            if (String.IsNullOrEmpty(record.Ref) ||
                record.Ref == ".")
            {
                record.IsMissing = true;
                return record;
            }

            if (record.AlleleIndex < 0)
            {
                record.IsMissing = true;
                return record;
            }

            if (record.AlleleIndex == 0)
            {
                record.SelectedAllele = record.Ref;
                return record;
            }

            string[] alternatives = record.Alt.Split(',');

            int altIndex = record.AlleleIndex - 1;

            if (altIndex < 0 || altIndex >= alternatives.Length)
            {
                record.IsMissing = true;
                return record;
            }

            record.SelectedAllele = alternatives[altIndex];

            if (record.SelectedAllele == "*")
            {
                record.IsStar = true;
                return record;
            }

            if (record.SelectedAllele.StartsWith("<") ||
                record.SelectedAllele.Contains("[") ||
                record.SelectedAllele.Contains("]"))
            {
                record.IsSymbolic = true;
                record.IsMissing = true;
            }

            return record;
        }

        private void AddDeletionMarkers(
            string reference,
            string alternate,
            int sequenceStartIndex,
            Dictionary<int, int> deletionMarkers)
        {
            if (String.IsNullOrEmpty(reference) ||
                String.IsNullOrEmpty(alternate) ||
                alternate == "*")
                return;

            int prefix = 0;

            while (prefix < reference.Length &&
                   prefix < alternate.Length &&
                   reference[prefix] == alternate[prefix])
            {
                prefix++;
            }

            int refEnd = reference.Length;
            int altEnd = alternate.Length;

            while (refEnd > prefix &&
                   altEnd > prefix &&
                   reference[refEnd - 1] == alternate[altEnd - 1])
            {
                refEnd--;
                altEnd--;
            }

            int refMiddle = refEnd - prefix;
            int altMiddle = altEnd - prefix;

            int deletedLength = refMiddle - altMiddle;

            if (deletedLength <= 0)
                return;

            int commonLength = Math.Min(refMiddle, altMiddle);

            int boundary = sequenceStartIndex + prefix + commonLength;

            if (deletionMarkers.ContainsKey(boundary))
                deletionMarkers[boundary] += deletedLength;
            else
                deletionMarkers.Add(boundary, deletedLength);
        }

        private void AppendHtmlRange(
            StringBuilder html,
            string sequence,
            int start,
            int length,
            string motifRegex,
            bool colorMotifs,
            Dictionary<int, int> deletionMarkers)
        {
            int end = start + length;

            Dictionary<int, int> blueClusterIds =
                new Dictionary<int, int>();

            if (colorMotifs && length > 0)
            {
                MatchCollection matches = Regex.Matches(
                    sequence.Substring(start, length),
                    motifRegex
                );

                int clusterId = 0;

                foreach (Match match in matches)
                {
                    if (!match.Success || match.Length == 0)
                        continue;

                    for (int i = 0; i < match.Length; i++)
                    {
                        int sequenceIndex = start + match.Index + i;

                        if (sequence[sequenceIndex] != 'N')
                        {
                            blueClusterIds[sequenceIndex] = clusterId;
                        }
                    }

                    clusterId++;
                }
            }

            int? activeBlueClusterId = null;

            for (int i = start; i < end; i++)
            {

                int deletionCount;

                if (deletionMarkers.TryGetValue(i, out deletionCount))
                {
                    if (activeBlueClusterId.HasValue)
                    {
                        html.Append("</span>");
                        activeBlueClusterId = null;
                    }

                    for (int d = 0; d < deletionCount; d++)
                    {
                        html.Append(
                            "<span style='" +
                            DeletionStyle +
                            "'>&middot;</span>"
                        );
                    }
                }

                char nucleotide = sequence[i];

                if (nucleotide == 'N')
                {
                    if (activeBlueClusterId.HasValue)
                    {
                        html.Append("</span>");
                        activeBlueClusterId = null;
                    }

                    html.Append(
                        "<span style='" +
                        GreyStyle +
                        "'>&middot;</span>"
                    );

                    continue;
                }

                int clusterId;

                if (blueClusterIds.TryGetValue(i, out clusterId))
                {
                    if (activeBlueClusterId != clusterId)
                    {
                        if (activeBlueClusterId.HasValue)
                            html.Append("</span>");

                        html.Append(
                            "<span style='" +
                            MotifStyle +
                            "'>"
                        );

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
            }

            if (activeBlueClusterId.HasValue)
                html.Append("</span>");
        }

        private void AppendFinalDeletionMarkers(
    StringBuilder html,
    Dictionary<int, int> deletionMarkers,
    int sequenceLength)
{
    int deletionCount;

    if (deletionMarkers.TryGetValue(
        sequenceLength, out deletionCount))
    {
        for (int i = 0; i < deletionCount; i++)
        {
            html.Append(
                "<span style='" +
                DeletionStyle +
                "'>&middot;</span>"
            );
        }
    }
}

        private YStrResult ReconstructSequence(
            string vcfOutput,
            long pos_start,
            long pos_end,
            string motif_regex,
            string filter_regex,
            string replace_regex)
        {
            YStrResult result = new YStrResult();

            Dictionary<long, VcfRecord> records =
                new Dictionary<long, VcfRecord>();

            HashSet<long> duplicatePositions =
                new HashSet<long>();

            using (StringReader reader = new StringReader(vcfOutput))
            {
                string line;

                while ((line = reader.ReadLine()) != null)
                {
                    if (String.IsNullOrWhiteSpace(line))
                        continue;

                    VcfRecord record = ParseRecord(line);

                    if (record == null)
                        continue;

                    if (record.Position < pos_start ||
                        record.Position > pos_end)
                        continue;

                    if (records.ContainsKey(record.Position))
                    {
                        duplicatePositions.Add(record.Position);
                        continue;
                    }

                    records.Add(record.Position, record);
                }
            }

            StringBuilder sequence = new StringBuilder();

            Dictionary<int, int> deletionMarkers =
                new Dictionary<int, int>();

            long regionLength = pos_end - pos_start + 1;
            long coveredPositions = 0;
            long currentPosition = pos_start;

            while (currentPosition <= pos_end)
            {
                VcfRecord record;

                if (!records.TryGetValue(currentPosition, out record))
                {
                    sequence.Append('N');
                    currentPosition++;
                    continue;
                }

                if (duplicatePositions.Contains(currentPosition))
                {
                    sequence.Append('N');
                    currentPosition++;
                    continue;
                }


                if (record.IsStar)
                {
                    int boundary = sequence.Length;

                    if (deletionMarkers.ContainsKey(boundary))
                        deletionMarkers[boundary]++;
                    else
                        deletionMarkers.Add(boundary, 1);

                    currentPosition++;
                    continue;
                }

                if (record.IsMissing)
                {
                    int refLength = record.Ref.Length;

                    for (int i = 0; i < refLength; i++)
                    {
                        if (currentPosition + i > pos_end)
                            break;

                        sequence.Append('N');
                    }

                    currentPosition += refLength;
                    continue;
                }

                string selectedAllele = record.SelectedAllele;

                if (String.IsNullOrEmpty(selectedAllele))
                {
                    sequence.Append('N');
                    currentPosition++;
                    continue;
                }

                int referenceLength = record.Ref.Length;

                long positionsInRange = Math.Min(
                    referenceLength,
                    pos_end - currentPosition + 1
                );

                coveredPositions += positionsInRange;

                int sequenceStartIndex = sequence.Length;

                AddDeletionMarkers(
                    record.Ref,
                    selectedAllele,
                    sequenceStartIndex,
                    deletionMarkers
                );

                sequence.Append(selectedAllele);

                currentPosition += referenceLength;
            }

            string seq = sequence.ToString();

            if (regionLength > 0)
            {
                result.reliability = (int)(
                    coveredPositions * 100 / regionLength
                );
            }
            else
            {
                result.reliability = 0;
            }

            result.raw_seq = seq;

            if (String.IsNullOrEmpty(motif_regex))
                motif_regex = "(?!)";

            if (String.IsNullOrEmpty(filter_regex))
                filter_regex = ".*";

            StringBuilder htmlBuilder = new StringBuilder();

            if (replace_regex == "CONCAT")
            {
                MatchCollection matches = Regex.Matches(
                    seq,
                    filter_regex
                );

                StringBuilder filteredSequence = new StringBuilder();

                int lastIndex = 0;

                foreach (Match match in matches)
                {
                    if (!match.Success)
                        continue;

                    AppendHtmlRange(
                        htmlBuilder,
                        seq,
                        lastIndex,
                        match.Index - lastIndex,
                        motif_regex,
                        false,
                        deletionMarkers
                    );

                    filteredSequence.Append(match.Value);

                    AppendHtmlRange(
                        htmlBuilder,
                        seq,
                        match.Index,
                        match.Length,
                        motif_regex,
                        true,
                        deletionMarkers
                    );

                    lastIndex = match.Index + match.Length;
                }

                AppendHtmlRange(
                    htmlBuilder,
                    seq,
                    lastIndex,
                    seq.Length - lastIndex,
                    motif_regex,
                    false,
                    deletionMarkers
                );

                seq = filteredSequence.ToString();
            }
            else
            {

                AppendHtmlRange(
                    htmlBuilder,
                    seq,
                    0,
                    seq.Length,
                    motif_regex,
                    true,
                    deletionMarkers
                );

                seq = Regex.Replace(
                    seq,
                    filter_regex,
                    replace_regex
                );
            }
            AppendFinalDeletionMarkers(
    htmlBuilder,
    deletionMarkers,
    sequence.Length
);
            result.html_seq = htmlBuilder.ToString();

            MatchCollection motifMatches = Regex.Matches(
                seq,
                motif_regex
            );

            result.value = motifMatches.Count;

            return result;
        }
    }
}
