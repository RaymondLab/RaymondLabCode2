import pytest

from raymondlab.apps import session_manager


def test_main_exits_with_code_2_and_says_not_built(capsys):
    with pytest.raises(SystemExit) as exc:
        session_manager.main()
    assert exc.value.code == 2
    assert "not built yet" in capsys.readouterr().err
