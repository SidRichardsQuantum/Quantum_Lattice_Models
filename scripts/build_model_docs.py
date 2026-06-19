"""Generate model-reference HTML pages and registry parameter tables."""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

from markdown_it import MarkdownIt

from quantum_lattice_models.registry import get_model_info

SOURCE_DIR = Path("docs/models")
BUILDERS_PATTERN = re.compile(r"<!--\s*builders:\s*([^>]+?)\s*-->")
PARAMETERS_TOKEN = "{{PARAMETERS}}"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Fail if generated HTML is stale")
    args = parser.parse_args(argv)
    markdown = MarkdownIt("commonmark", {"html": True, "linkify": True}).enable("table")
    stale = []
    for source in sorted(SOURCE_DIR.glob("*.md")):
        text = source.read_text()
        builders_match = BUILDERS_PATTERN.search(text)
        builders = ()
        if builders_match:
            builders = tuple(name.strip() for name in builders_match.group(1).split(","))
            text = BUILDERS_PATTERN.sub("", text)
        if PARAMETERS_TOKEN in text:
            text = text.replace(PARAMETERS_TOKEN, _parameter_markdown(builders))
        text = _web_links(text, model_page=True)
        body = markdown.render(text)
        output = source.with_suffix(".html")
        rendered = _page(source.stem, _title(text), body, "../styles.css")
        if args.check:
            if not output.exists() or output.read_text() != rendered:
                stale.append(output)
        else:
            output.write_text(rendered)

    theory_source = Path("THEORY.md")
    theory_text = _web_links(theory_source.read_text(), model_page=False)
    theory_output = Path("docs/theory.html")
    theory_rendered = _page(
        "theory", _title(theory_text), markdown.render(theory_text), "styles.css"
    )
    if args.check:
        if not theory_output.exists() or theory_output.read_text() != theory_rendered:
            stale.append(theory_output)
    else:
        theory_output.write_text(theory_rendered)
    if stale:
        paths = ", ".join(str(path) for path in stale)
        raise SystemExit(f"Generated model documentation is stale: {paths}")


def _parameter_markdown(builders: tuple[str, ...]) -> str:
    if not builders:
        return "_No registered builders are associated with this page._"
    lines = [
        "| Builder | Parameter | Type | Default | Constraint |",
        "|---|---|---:|---:|---|",
    ]
    for builder in builders:
        info = get_model_info(builder)
        for parameter in info.parameters:
            constraint = f">= {parameter.minimum:g}" if parameter.minimum is not None else ""
            row = (
                f"| `{builder}` | `{parameter.name}` | "
                f"`{parameter.type.__name__}` | `{parameter.default}` | {constraint} |"
            )
            lines.append(row)
    return "\n".join(lines)


def _title(text: str) -> str:
    for line in text.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return "Model Reference"


def _web_links(text: str, *, model_page: bool) -> str:
    text = re.sub(r"\(([^)]+)\.md\)", r"(\1.html)", text)
    if model_page:
        text = text.replace("../../THEORY.html", "../theory.html")
    else:
        text = text.replace("docs/models/", "models/")
        text = text.replace(
            "ROADMAP.html",
            "https://github.com/SidRichardsQuantum/Quantum_Lattice_Models/blob/main/ROADMAP.md",
        )
    return text


def _page(slug: str, title: str, body: str, stylesheet: str) -> str:
    home = "../index.html" if slug != "theory" else "index.html"
    index = "index.html" if slug != "theory" else "models/index.html"
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="color-scheme" content="light dark">
    <title>{html.escape(title)} — Quantum Lattice Models</title>
    <link rel="stylesheet" href="{stylesheet}">
    <script>
      MathJax = {{tex: {{inlineMath: [['$', '$']], displayMath: [['$$', '$$']]}}}};
    </script>
    <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  </head>
  <body>
    <main class="model-reference">
      <nav class="model-breadcrumbs">
        <a href="{home}">Home</a>
        <span>/</span>
        <a href="{index}">Model reference</a>
        <span>/</span>
        <span>{html.escape(title)}</span>
      </nav>
      <article class="model-article" data-model-page="{html.escape(slug)}">
{body}
      </article>
    </main>
  </body>
</html>
"""


if __name__ == "__main__":
    main()
