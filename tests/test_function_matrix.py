import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'bin/rmetag_function_matrix.py'


class FunctionMatrixTests(unittest.TestCase):
    def test_fractional_mapping_and_missing_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'map.tsv').write_text('gene_id\tfunction_id\ng1\tK1\ng1\tK2\ng1\tK2\ng2\tK1\n')
            sample = root / 's1'
            sample.mkdir()
            (sample / 'abundance.tsv').write_text('target_id\ttpm\ng1\t60\ng2\t30\ng3\t10\n')
            subprocess.run([sys.executable, str(SCRIPT), '--mapping', str(root / 'map.tsv'), '--quant', str(sample)], cwd=root, check=True)
            with (root / 'function_fractional_tpm.tsv').open() as f:
                rows = list(csv.reader(f, delimiter='\t'))
            self.assertEqual(rows, [['feature_id', 's1'], ['K1', '60'], ['K2', '30']])
            with (root / 'annotation_summary.tsv').open() as f:
                summary = list(csv.reader(f, delimiter='\t'))
            self.assertEqual(summary[1], ['s1', '3', '2', '100.0', '90.0'])

    def test_invalid_mapping_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'bad.tsv').write_text('gene\tannotation\ng1\tK1\n')
            result = subprocess.run([sys.executable, str(SCRIPT), '--mapping', str(root / 'bad.tsv'), '--quant', str(root)], cwd=root, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b'gene_id and function_id', result.stderr)


if __name__ == '__main__':
    unittest.main()
