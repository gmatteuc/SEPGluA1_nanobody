# History

Records of the refactor of September and October 2026, kept for the reasons
behind decisions and for the numbers the checks gave. They describe the code
as it was then and are not updated; what is still open is in
[ROADMAP.md](../ROADMAP.md). Paths on G: in them are working folders of the
refactor and may have been cleared since.

| file | what |
|---|---|
| [REFACTOR_PLAN.md](REFACTOR_PLAN.md) | the plan: rules, layout, decisions L1 to L7, Y1 to Y7 and S1 to S6, the verification design, the bug list and the progress, 30 September to 5 October |
| [production_settings.md](production_settings.md) | the settings, driver by driver, of the reference run of step 3 and of every check after it (30 September) |
| [HAND_CHECK.md](HAND_CHECK.md) | the hand check of the GUIs after step 6 (2 October) |
| [FIXES_STEP8.md](FIXES_STEP8.md) | every item of the step 8 list: applied, held for a decision, or skipped, with the reason (2 and 3 October) |
| [MORNING_REPORT.md](MORNING_REPORT.md) | the held fixes and Giulio's answers to them (3 October) |
| [STEP8_REPORT.md](STEP8_REPORT.md) | the fixes of step 8 as applied, what each changed and how it was checked (4 October) |
| [DECISIONS_REPORT.md](DECISIONS_REPORT.md) | the four decisions of 4 October on the plasticity comparison and the Python route, and their measurements |

The old and new names of every script are in
[refactor_name_map.csv](../refactor_name_map.csv), which stays in `docs/`:
`tools/check_code_identity.py --map` reads it, and the warning of
`sep_setup_paths` points to it.
