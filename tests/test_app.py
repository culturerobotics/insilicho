import insilicho as package


class TestSuite:
    def test_version(self) -> None:
        assert package.__version__ == "1.0.0"
