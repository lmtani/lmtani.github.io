import os
import pathlib
import pytest

@pytest.mark.workflow('Run playing-with-bioinformatics example')
def test_ensure_crumble_produces_small_cram(workflow_dir):
    """crumble somethimes can fail and still present exit_code 0
    This test is to ensure the created CRAM has data"""
    crumble = pathlib.Path(workflow_dir, "STEP_4/alignment_nv9_lossnames.cram")
    file_stats = os.stat(crumble)
    size_in_mb = file_stats.st_size / (1024 * 1024)
    assert size_in_mb > 3.0
    assert size_in_mb < 3.5