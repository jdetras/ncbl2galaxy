import tempfile
import unittest
from pathlib import Path

import ncbi_to_galaxy


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, payload):
        self.payload = payload

    def get(self, *args, **kwargs):
        return FakeResponse(self.payload)


class NCBIToGalaxyTests(unittest.TestCase):
    def test_chunked(self):
        data = list(range(7))
        self.assertEqual(ncbi_to_galaxy.chunked(data, 3), [[0, 1, 2], [3, 4, 5], [6]])

    def test_ena_run_record_normalization(self):
        payload = [
            {
                "fastq_ftp": "ftp.sra.ebi.ac.uk/path1.fastq.gz;https://x.test/path2.fastq.gz",
                "sample_accession": "SAMEA1",
                "library_layout": "PAIRED",
            }
        ]
        session = FakeSession(payload)
        record = ncbi_to_galaxy.ena_run_record(session, "SRR123")
        self.assertIsNotNone(record)
        self.assertEqual(record.sample_accession, "SAMEA1")
        self.assertEqual(record.library_layout, "PAIRED")
        self.assertEqual(
            record.fastq_urls,
            [
                "https://ftp.sra.ebi.ac.uk/path1.fastq.gz",
                "https://x.test/path2.fastq.gz",
            ],
        )

    def test_ena_run_record_empty(self):
        session = FakeSession([{}])
        self.assertIsNone(ncbi_to_galaxy.ena_run_record(session, "SRR123"))

    def test_group_runs_by_sample(self):
        records = [
            ncbi_to_galaxy.RunRecord("SRR1", "S1", "SINGLE", ["u1"]),
            ncbi_to_galaxy.RunRecord("SRR2", "S1", "SINGLE", ["u2"]),
            ncbi_to_galaxy.RunRecord("SRR3", "S2", "PAIRED", ["u3", "u4"]),
        ]
        grouped = ncbi_to_galaxy.group_runs_by_sample(records)
        self.assertEqual(sorted(grouped.keys()), ["S1", "S2"])
        self.assertEqual(len(grouped["S1"]), 2)
        self.assertEqual(len(grouped["S2"]), 1)

    def test_state_store_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "state.json"
            store = ncbi_to_galaxy.StateStore(str(state_path))
            store.processed_runs = {"SRR3", "SRR1"}
            store.save()

            loaded = ncbi_to_galaxy.StateStore(str(state_path))
            loaded.load()
            self.assertEqual(loaded.processed_runs, {"SRR1", "SRR3"})


if __name__ == "__main__":
    unittest.main()
