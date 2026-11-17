import hipercam as hcam
from hipercam.scripts.checkflags import checkflags


def test_checkflags_log(tmp_path, monkeypatch, capsys):
    log = tmp_path / "run.log"
    log.write_text(
        "# Start of column name definitions:\n"
        "# 1 = CCD MJDok flag_target flag_comp\n"
        "# 2 = CCD flag_target\n"
        "# 3 = CCD flag_empty\n"
        "# End of column name definitions\n"
        "# Start of data type definitions:\n"
        "# 1 = s ? u4 u4\n"
        "# 2 = s u4\n"
        "# 3 = s u4\n"
        "# End of data type definitions\n"
        f"1 1 {int(hcam.NO_SKY | hcam.BAD_TIME)} 0\n"
        f"1 0 {int(hcam.NO_SKY)} 0\n"
        "1 1 2147483648 0\n"
        f"2 {int(hcam.BAD_TIME)}\n"
    )
    monkeypatch.setenv("HIPERCAM_ENV", str(tmp_path))

    checkflags(["checkflags", str(log.with_suffix(""))])

    output = capsys.readouterr().out
    target = output.split("CCD 1, aperture target: 3/3 flagged\n")[1].split("CCD 2")[0]
    assert "  NO_SKY: 2 - zero sky pixels" in target
    assert "  BAD_TIME: 1 - the time was flagged as bad" in target
    assert "  UNKNOWN 0x80000000: 1" in target
    assert "CCD 1, aperture comp: 0/3 flagged\nCCD 1, aperture target:" in output
    assert "CCD 2, aperture target: 1/1 flagged\n  BAD_TIME: 1" in output
    assert "CCD 3, aperture empty: 0/0 flagged" in output
    assert "ANY_FLAG:" not in output
    assert "ALL_OK:" not in output
