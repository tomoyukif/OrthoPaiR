---
description: Update DESCRIPTION and inst/NEWS.md with a version bump
---

Update `DESCRIPTION` and `inst/NEWS.md`.

Rules:
- If the user explicitly provides a target version, use that version.
- If no version is specified, increment the 3rd segment of the current version by 1
  (e.g. `0.7.4` -> `0.7.5`).
- Update `DESCRIPTION`:
  - `Version: <new_version>`
  - Keep `Date` as today's date.
- Update `inst/NEWS.md`:
  - Add a new top entry:
    - `Changes in version <new_version> (<today>)`
    - A short bullet list summarizing the current changes.
- Do not modify unrelated files unless needed to complete these updates.
- After edits, show the updated version line and the newly added NEWS header.
