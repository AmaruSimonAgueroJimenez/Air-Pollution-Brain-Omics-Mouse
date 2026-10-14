import hashlib
import http.server
import importlib.util
from pathlib import Path
import tempfile
import threading
import unittest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/download_data.py"
SPEC = importlib.util.spec_from_file_location("downloader", MODULE_PATH)
DOWNLOADER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DOWNLOADER)
PAYLOAD = b"Scientific sample bytes\n" * 1000


class RangeHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        start = int(self.headers.get("Range", "bytes=0-").split("=")[1].split("-")[0])
        self.send_response(206 if start else 200)
        if start:
            self.send_header("Content-Range", f"bytes {start}-{len(PAYLOAD)-1}/{len(PAYLOAD)}")
        self.send_header("Content-Length", str(len(PAYLOAD) - start))
        self.end_headers()
        self.wfile.write(PAYLOAD[start:])

    def log_message(self, *args):
        pass


class IntegrityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        DOWNLOADER.ROOT = Path(self.temporary.name)
        (DOWNLOADER.ROOT / "logs").mkdir()
        DOWNLOADER.RESULTS.clear()
        self.row = {"dataset": "fixture", "kind": "processed", "path": "data/fixture/counts.txt",
                    "url": f"http://127.0.0.1:{self.server.server_address[1]}/counts.txt",
                    "bytes": len(PAYLOAD), "md5": hashlib.md5(PAYLOAD).hexdigest()}

    def tearDown(self):
        self.temporary.cleanup()

    def test_resume_then_atomic_verified_file(self):
        destination = DOWNLOADER.ROOT / self.row["path"]
        destination.parent.mkdir(parents=True)
        partial = destination.with_name(destination.name + ".part")
        partial.write_bytes(PAYLOAD[:321])
        self.assertTrue(DOWNLOADER.download(self.row))
        self.assertEqual(destination.read_bytes(), PAYLOAD)
        self.assertFalse(partial.exists())
        self.assertEqual(DOWNLOADER.RESULTS[self.row["path"]]["state"], "verified")

    def test_wrong_source_checksum_cannot_be_marked_complete(self):
        self.row["md5"] = "0" * 32
        self.assertFalse(DOWNLOADER.download(self.row))
        self.assertFalse((DOWNLOADER.ROOT / self.row["path"]).exists())
        self.assertEqual(DOWNLOADER.RESULTS[self.row["path"]]["state"], "failed")

    def test_corrupted_partial_is_quarantined_and_refetched(self):
        destination = DOWNLOADER.ROOT / self.row["path"]
        destination.parent.mkdir(parents=True)
        destination.with_name(destination.name + ".part").write_bytes(b"X" * 321)
        self.assertTrue(DOWNLOADER.download(self.row))
        self.assertEqual(destination.read_bytes(), PAYLOAD)
        self.assertEqual(len(list(destination.parent.glob("*.invalid.*"))), 1)

    def test_invalid_existing_destination_is_recovered(self):
        destination = DOWNLOADER.ROOT / self.row["path"]
        destination.parent.mkdir(parents=True)
        destination.write_bytes(b"invalid existing bytes")
        self.assertTrue(DOWNLOADER.download(self.row))
        self.assertEqual(destination.read_bytes(), PAYLOAD)
        self.assertEqual(len(list(destination.parent.glob("*.invalid.*"))), 1)

    def test_html_is_not_accepted_as_a_dataset(self):
        destination = DOWNLOADER.ROOT / "counts.txt"
        destination.write_bytes(b"<!DOCTYPE html><html>Access error</html>")
        with self.assertRaisesRegex(ValueError, "HTML"):
            DOWNLOADER.validate({"bytes": None, "md5": None, "path": "counts.txt"}, destination)


if __name__ == "__main__":
    unittest.main()
