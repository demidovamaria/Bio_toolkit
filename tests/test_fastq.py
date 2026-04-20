import pytest
import logging
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main import filter_fastq


class TestFastq:

    def test_output_file_created(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nATGC\n+\nIIII\n")

        filter_fastq(str(input_file), str(output_file))

        assert output_file.exists()


    def test_filter_by_length(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nATGC\n+\nIIII\n")

        filter_fastq(
            str(input_file),
            str(output_file),
            length_bounds=(5, 10)
        )

        assert output_file.read_text() == ""


    def test_filter_by_gc(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nAAAA\n+\nIIII\n")

        filter_fastq(
            str(input_file),
            str(output_file),
            gc_bounds=(0.5, 1.0)
        )

        assert output_file.read_text() == ""


    def test_filter_by_quality(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nATGC\n+\n!!!!\n")

        filter_fastq(
            str(input_file),
            str(output_file),
            quality_threshold=30
        )

        assert output_file.read_text() == ""


    def test_valid_read_passes(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nATGC\n+\nIIII\n")

        filter_fastq(
            str(input_file),
            str(output_file),
            gc_bounds=(0.0, 1.0),
            length_bounds=(1, 10),
            quality_threshold=0
        )

        content = output_file.read_text()
        assert "@r1" in content


    def test_all_filtered_empty_file(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("@r1\nAAAA\n+\n!!!!\n")

        filter_fastq(
            str(input_file),
            str(output_file),
            gc_bounds=(0.9, 1.0),
            quality_threshold=40
        )

        assert output_file.exists()
        assert output_file.read_text() == ""


class TestErrors:

    def test_input_file_not_found(self, tmp_path):
        output_file = tmp_path / "output.fastq"

        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent.fastq", str(output_file))


class TestLogging:

    def test_logging_file_created(self, tmp_path):
        log_file = tmp_path / "test.log"

        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            force=True  # важно!
        )

        logging.info("test message")

        assert log_file.exists()
        assert "test message" in log_file.read_text()