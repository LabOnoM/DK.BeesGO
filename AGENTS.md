This repository hosts a knowledge base website built with Quartz. If you change code in the `source` directory, run `npm run check` and `npm test` inside that folder before committing. You can preview the generated site locally with `npx quartz build --serve`.

For documentation or Markdown content only, tests are not required. Follow the community Code of Conduct in `source/CODE_OF_CONDUCT.md` when contributing.

## Agent Governance & Workspace Policies

To maintain auditability, knowledge grounding, and workspace cleanliness, all AI agents operating in this repository MUST strictly follow these rules:

### 1. Wiki-First Resolution & Grounding
- Before proposing any design or answering domain-specific scientific/bioinformatics questions, query the local LLM-Wiki (`.wiki/` or symlink `Wiki`).
- State clearly if information is not found in the wiki. Use double-bracket links `[[PageName]]` for citing wiki pages.

### 2. Workspace Hygiene (Scratch Files)
- Ad-hoc scripts, debug code, or temporary files **MUST NOT** be saved in the root or source folders.
- All scratch files must be placed in the designated workspace scratch directory: `~/.gemini/antigravity/brain/<conversation-id>/scratch/`.

### 3. Large Artifact Generation Policy
- Direct generation of massive text documents, manuscripts, or reports (exceeding 3000 words) inside the LLM context is prohibited due to token constraints and quality issues.
- Agents must instead write Python scripts or programmatic compilers that combine structured sections from `.wiki/` or text chunks into the final document.

### 4. Gated File I/O & Audit Layer
- All agent actions must run under the `re_gent` VCS layer.
- Modifications to existing files must be completed via the `/lab-commit` workflow.
