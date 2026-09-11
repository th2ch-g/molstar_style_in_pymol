# Development

- Reply in concise Japanese; write code comments and commit messages in English.
- Run Python through uv with the pixi interpreter.
- Keep runtime independent of Node.js, Mol*, CueMol, and network services.
- Read docs/guide.md before changing the command interface or renderer.
- Preserve source atoms and unrelated views; never patch PyMOL commands.
- Validate real GPU and ray output in a separate PyMOL process.
- Keep caches, environments, downloaded inputs, and temporary renders ignored.
- Version published gallery PNGs in docs/gallery; render them locally, without CI.
- Maintain English and Japanese guides together. Keep paths portable.
- Scope searches to relevant files. Preserve third-party attribution in NOTICE.
- Append a blank line and Co-authored-by: Codex <noreply@openai.com> to commits.
