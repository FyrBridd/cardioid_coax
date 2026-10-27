import pytest

from fakes import DesignTuple, default_design


@pytest.fixture
def design() -> DesignTuple:
    return default_design()
