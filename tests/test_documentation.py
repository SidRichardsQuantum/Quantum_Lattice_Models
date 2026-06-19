from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path

from scripts.build_model_docs import main


def test_generated_model_documentation_is_current() -> None:
    main(["--check"])


def test_every_model_source_has_generated_page() -> None:
    sources = sorted(Path("docs/models").glob("*.md"))
    assert len(sources) >= 19
    for source in sources:
        generated = source.with_suffix(".html")
        assert generated.exists()
        text = generated.read_text()
        assert "<h1>" in text
        assert "{{PARAMETERS}}" not in text


def test_generated_model_pages_have_valid_local_links() -> None:
    class References(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.values: list[str] = []

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            for name, value in attrs:
                if name in {"href", "src"} and value is not None:
                    self.values.append(value)

    pages = [Path("docs/theory.html"), *Path("docs/models").glob("*.html")]
    for page in pages:
        parser = References()
        parser.feed(page.read_text())
        for reference in parser.values:
            if reference.startswith(("http://", "https://", "#")):
                continue
            target = (page.parent / reference.split("#", 1)[0]).resolve()
            assert target.exists(), f"{page}: missing {reference}"


def test_schema_and_import_guides_are_published() -> None:
    schema = Path("SCHEMA.md").read_text()
    importing = Path("IMPORTING.md").read_text()

    assert "Compatibility policy" in schema
    assert "Complex-number encoding" in schema
    assert "CSV import" in importing
    assert "NetworkX and GraphML" in importing
