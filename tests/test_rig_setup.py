import pytest

from raymondlab.apps import rig_setup


def test_main_exits_with_code_2_and_says_not_built(capsys):
    with pytest.raises(SystemExit) as exc:
        rig_setup.main()
    assert exc.value.code == 2
    out = capsys.readouterr()
    assert "not built yet" in (out.out + out.err)
