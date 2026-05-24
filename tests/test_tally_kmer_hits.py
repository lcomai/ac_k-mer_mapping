import csv
import tempfile
import unittest
from pathlib import Path

from tally_kmer_hits import main


class TallyTests(unittest.TestCase):
    def test_tallies_mapping_rows_per_kmer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            mapping = tmp / "mapped.tsv"
            out = tmp / "counts.tsv"

            mapping.write_text(
                "kmer\tseq_id\tstart\tend\tstrand\tsource_fasta\n"
                "ACGA\tChr1\t4\t7\t+\treference.fa\n"
                "ACGA\tChr1\t8\t11\t-\treference.fa\n"
                "GAAA\tChr1\t13\t16\t+\treference.fa\n",
                encoding="utf-8",
            )

            self.assertEqual(main([str(mapping), str(out)]), 0)

            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(row["kmer"], row["mapped_positions"]) for row in rows],
                [("ACGA", "2"), ("GAAA", "1")],
            )

    def test_includes_zero_counts_from_original_kmer_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            mapping = tmp / "mapped.tsv"
            kmers = tmp / "kmers.txt"
            out = tmp / "counts.tsv"

            mapping.write_text(
                "k-mer\tchr\tpos\tstrand\n"
                "ACGA\tChr1\t4\t+\n"
                "ACGA\tChr1\t8\t-\n",
                encoding="utf-8",
            )
            kmers.write_text("ACGA\nTTTT\n", encoding="utf-8")

            self.assertEqual(
                main([str(mapping), str(out), "--kmers", str(kmers)]), 0
            )

            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [(row["kmer"], row["mapped_positions"]) for row in rows],
                [("ACGA", "2"), ("TTTT", "0")],
            )


if __name__ == "__main__":
    unittest.main()
