import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/databases.py'
spec = importlib.util.spec_from_file_location('rmetag_databases', SCRIPT)
db = importlib.util.module_from_spec(spec)
spec.loader.exec_module(db)


class DatabaseHelperTests(unittest.TestCase):
    def test_verified_fetch_and_refuse_overwrite(self):
        payload = b'reference index'
        checksum = 'sha256:' + hashlib.sha256(payload).hexdigest()
        with tempfile.TemporaryDirectory() as tmp, patch.object(db, 'request', side_effect=lambda _: io.BytesIO(payload)):
            dest = Path(tmp) / 'index.syldb'
            db.fetch('https://example.org/index', dest, checksum, 'test')
            self.assertEqual(dest.read_bytes(), payload)
            manifest = json.loads(Path(str(dest) + '.manifest.json').read_text())
            self.assertEqual(manifest['sha256'], checksum.split(':')[1])
            with self.assertRaises(FileExistsError):
                db.fetch('https://example.org/index', dest, checksum, 'test')

    def test_checksum_failure_leaves_no_output(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(db, 'request', side_effect=lambda _: io.BytesIO(b'bad')):
            dest = Path(tmp) / 'index'
            with self.assertRaises(ValueError):
                db.fetch('https://example.org/index', dest, 'sha256:' + '0' * 64, 'test')
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
