"""Build the compact executed-notebook summary used by repository documentation."""

from __future__ import annotations

import json
from pathlib import Path

import nbformat


def main() -> None:
    summaries = []
    for path in sorted(Path("notebooks").glob("*.ipynb")):
        notebook = nbformat.read(path, as_version=4)
        text: list[str] = []
        images = 0
        for cell in notebook.cells:
            if cell.cell_type != "code":
                continue
            for output in cell.get("outputs", []):
                if output.output_type == "stream":
                    text.extend(line for line in output.get("text", "").splitlines() if line)
                elif output.output_type in {"display_data", "execute_result"}:
                    data = output.get("data", {})
                    if "image/png" in data:
                        images += 1
                    plain = data.get("text/plain", "")
                    if plain and not plain.startswith("<Figure"):
                        text.extend(line for line in plain.splitlines() if line)
        summaries.append(
            {
                "notebook": path.name,
                "html": f"docs/notebooks/{path.stem}.html",
                "images": images,
                "text": text[:15],
            }
        )

    output = Path("results/notebook_summary.json")
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    main()
