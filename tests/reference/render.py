"""Capture an actual Mol* browser render and its numerical geometry."""

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8765")
    parser.add_argument("--output", type=Path, default=Path(".cache/reference"))
    parser.add_argument("--options", type=Path)
    parser.add_argument("--name", default="molstar")
    parser.add_argument("--browser", help="Optional Chromium executable")
    parser.add_argument(
        "--suite", type=Path, help="JSON array of named comparison cases"
    )
    args = parser.parse_args()
    options = json.loads(args.options.read_text()) if args.options else {}
    cases = (
        json.loads(args.suite.read_text())
        if args.suite
        else [{"name": args.name, "options": options}]
    )
    args.output.mkdir(parents=True, exist_ok=True)
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=args.browser,
            args=[
                "--use-gl=angle",
                "--use-angle=swiftshader",
                "--enable-unsafe-swiftshader",
            ],
        )
        page = browser.new_page(
            viewport={"width": 1200, "height": 900}, device_scale_factor=1
        )
        page.on(
            "pageerror",
            lambda error: errors.append(str(error.stack)) if len(errors) < 3 else None,
        )
        page.goto(args.url)
        page.wait_for_load_state("networkidle")
        for case in cases:
            name, options = case["name"], case["options"]
            page.evaluate("options => reference.load(options)", options)
            page.wait_for_timeout(500)
            page.screenshot(path=str(args.output / f"{name}.png"))
            dump = page.evaluate("reference.dump()")
            dump["options"] = options
            (args.output / f"{name}.json").write_text(json.dumps(dump) + "\n")
            print(
                json.dumps(
                    {
                        "errors": errors,
                        "segments": len(dump["segments"]),
                        "meshes": len(dump["meshes"]),
                    }
                )
            )
        browser.close()
    if errors:
        raise RuntimeError("; ".join(errors))


if __name__ == "__main__":
    main()
