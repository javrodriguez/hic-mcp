#!/usr/bin/env python
"""Structural checks on the README that a reader would notice before anything else.

1. The client config block appears within the first 200 words, so a screener can
   paste it without scrolling.
2. That block is byte-identical to examples/claude_desktop_config.json, so the two
   cannot drift apart.

Run directly (`python scripts/check_readme.py`) or via the test suite.
"""

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CONFIG_WORD_LIMIT = 200


def readme_text() -> str:
    return (REPO / "README.md").read_text()


def config_block() -> dict:
    """The first ```json fenced block in the README that declares mcpServers."""
    for block in re.findall(r"```json\n(.*?)```", readme_text(), re.DOTALL):
        parsed = json.loads(block)
        if "mcpServers" in parsed:
            return parsed
    raise AssertionError("README has no mcpServers json block")


def check_config_is_early() -> int:
    text = readme_text()
    idx = text.find('"mcpServers"')
    assert idx > 0, "README has no mcpServers block"
    words = len(text[:idx].split())
    assert words <= CONFIG_WORD_LIMIT, (
        f"the client config sits {words} words in; it must appear within the first "
        f"{CONFIG_WORD_LIMIT}"
    )
    return words


def check_config_matches_example() -> None:
    example = json.loads((REPO / "examples" / "claude_desktop_config.json").read_text())
    assert config_block() == example, (
        "the README's config block and examples/claude_desktop_config.json have drifted apart"
    )


def main() -> None:
    words = check_config_is_early()
    check_config_matches_example()
    print(f"README ok: client config at word {words} (limit {CONFIG_WORD_LIMIT}); "
          "matches examples/claude_desktop_config.json")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        sys.exit(f"README check failed: {e}")
