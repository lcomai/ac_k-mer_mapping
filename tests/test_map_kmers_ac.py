import csv
import tempfile
import unittest
from pathlib import Path

from map_kmers_ac import main


class MapperTests(unittest.TestCase):
    def test_maps_forward_and_reverse_complement_hits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            kmers = tmp / "kmers.txt"
            fasta = tmp / "reference.fa"
            out = tmp / "mapped.tsv"

            kmers.write_text("ACGA\nGAAA\nTTTT\n", encoding="utf-8")
            fasta.write_text(
                ">Chr1 example chromosome\nTTTACGATCGTGGAAAC\n"
                ">scaffold_2\nCCCCGAAATTTT\n",
                encoding="utf-8",
            )

            self.assertEqual(main([str(kmers), str(fasta), str(out)]), 0)

            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            observed = {
                (
                    row["kmer"],
                    row["seq_id"],
                    int(row["start"]),
                    int(row["end"]),
                    row["strand"],
                )
                for row in rows
            }
            self.assertIn(("ACGA", "Chr1", 4, 7, "+"), observed)
            self.assertIn(("ACGA", "Chr1", 8, 11, "-"), observed)
            self.assertIn(("GAAA", "Chr1", 13, 16, "+"), observed)
            self.assertIn(("GAAA", "scaffold_2", 5, 8, "+"), observed)

            unmatched = (tmp / "mapped.tsv.no_match").read_text(encoding="utf-8")
            self.assertEqual(unmatched, "")

    def test_directory_input_and_full_name_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            genome_dir = tmp / "genomes"
            genome_dir.mkdir()
            kmers = tmp / "kmers.txt"
            out = tmp / "mapped.tsv"

            kmers.write_text("ACGA\nBAD_HEADER\n", encoding="utf-8")
            (genome_dir / "a.fa").write_text(
                ">Chr1 full name\nACGA\n", encoding="utf-8"
            )

            self.assertEqual(
                main([str(kmers), str(genome_dir), str(out), "--name-mode", "full"]),
                0,
            )

            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["seq_id"], "Chr1 full name")

            invalid = (tmp / "mapped.tsv.invalid.tsv").read_text(encoding="utf-8")
            self.assertIn("BAD_HEADER", invalid)


if __name__ == "__main__":
    unittest.main()
