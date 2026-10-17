# Target venue and submission requirements

**Provenance warning.** The sandbox this study was run in cannot open journal
websites (DNS/egress is blocked for `emerginginvestigators.org`, `peerj.com`,
`pmc.ncbi.nlm.nih.gov`, etc.), so every requirement below was gathered from
web-search result summaries of the journals' own pages, **not** from the pages
themselves. Items marked **[VERIFY]** could not be confirmed at all. Open the
live pages (linked) and the official manuscript template before submitting;
journals change these rules often (one JEI page we saw was dated Aug 2025,
another Jan 2022, and the submission fee changed between pages).

## Recommendation

| Role | Venue | Why |
|---|---|---|
| **Primary** | [Journal of Emerging Investigators (JEI)](https://emerginginvestigators.org/) | Peer-reviewed, free to read, run for middle/high-school authors; the pasted idea list explicitly targets USABO/ISEF-level students. Fit is conditional on the hypothesis-driven rule below, so the paper is written hypothesis-first. |
| **Backup (peer-reviewed, any age)** | [PeerJ](https://peerj.com/about/editorial-criteria) (Life & Environment / Biochemistry & molecular biology sections) | Judges *soundness, not novelty or importance*; accepts computational research articles. Open-access fee applies **[VERIFY amount]**. |
| **Preprint (no review, free)** | bioRxiv, or JEI's own student preprint server eiRxiv **[VERIFY eligibility]** | Establishes a public timestamp while review is pending. |

## JEI requirements (from search summaries of JEI's own pages)

**Eligibility and process**
- Authors are middle- or high-school students who have not started college;
  students must be 13+ to submit. An adult must submit on the student's behalf
  (students may not submit their own manuscript); a senior mentor is listed
  last; undergraduates cannot be mentors.
  ([FAQ](https://emerginginvestigators.org/submissions/faq),
  [Pay and Submit](https://emerginginvestigators.org/submissions/checklist))
- Fee: **US$49** per submission on the current pay page (an older blog post
  said $35); nonrefundable; waivers/scholarships exist; one submission per
  student at a time.
  ([Pay the fee](https://emerginginvestigators.org/submissions/pay-and-submit/pay))
- Process: managing-editor formatting check (~2 weeks), then scientific review
  by 3-4 reviewers, usually several rounds.
- Submission is through Editorial Manager; the mandatory template is on the
  [author resources page](https://emerginginvestigators.org/documents).
  Manuscripts not on the template are returned.

**Scope (the biggest risk for an in-silico paper)**
- [Topics We Accept](https://emerginginvestigators.org/submissions/hypothesis-driven-research):
  *all* JEI manuscripts must be **hypothesis-driven**; the hypothesis must be
  specific and measurable ("We hypothesized ..."). Mathematical modelling is
  acceptable *if* a hypothesis is tested.
- [Guidelines for engineering and machine-learning projects](https://emerginginvestigators.org/submissions/engineering-and-machine-learning-based-projects)
  (updated 15 Aug 2025): tool-building, method-optimisation and "I built a
  model" papers are rejected; the page says JEI currently lacks the expertise
  to evaluate the growing number of computational manuscripts.
  **Implication:** this paper is framed as tests of stated hypotheses about
  binding, with a retrospective benchmark and pre-specified success criteria,
  not as a tool or a model. Email the editors
  (submissions@emerginginvestigators.org) to ask whether a molecular-modelling
  study fits **before** paying the fee.
- Ethics: ISEF-style approval is needed only for human subjects / vertebrate
  animals (not applicable here: no human or animal data).
- [Academic Honesty and AI](https://emerginginvestigators.org/submissions/academic-honesty-and-ai):
  JEI has a dedicated policy page; violations can lead to pre-review failure,
  public withdrawal, or a ban from JEI and eiRxiv. **[VERIFY: we could not read
  the rule text. Read it before submitting; this manuscript was drafted with an
  AI assistant and the student author must be able to explain and defend every
  step, and must comply with whatever disclosure/prohibition the page states.]**

**Format** ([Manuscript Format & Content](https://emerginginvestigators.org/submissions/manuscript-content),
[Reference Formatting](https://emerginginvestigators.org/submissions/references),
[Preparing for Submission](https://emerginginvestigators.org/submissions/preparing-for-submission))
- Section order, exactly: **Summary (abstract, < 250 words) -> Introduction ->
  Results -> Discussion -> Materials and Methods -> References.** No extra
  sections; subheadings allowed only in Results and Methods.
- Body max **10 pages** (excludes title page, references, figures, captions);
  Times New Roman 11 pt, 1.5 line spacing, 1-inch margins.
- Title <= **110 characters** including spaces.
- Figures go **at the end** of the Word file, each above its caption; tables are
  editable Word tables numbered separately from figures.
- References: *modified MLA8*, numbered in order of first citation, no hanging
  indent, no publisher for journal articles, "first author et al." for 3+
  authors. Every reference must be cited in the text.
- No Word equation objects in-line; use Insert > Symbol.
- Outside data and software must be cited.
- Word limits for Introduction/Discussion/figure legends: **[VERIFY]**.

## PeerJ (backup) -- what we could confirm
- Reviewers judge *soundness, not importance or novelty*
  ([editorial criteria](https://peerj.com/about/editorial-criteria)), using
  Basic Reporting / Experimental Design / Validity of the Findings.
- Research Articles and Bioinformatics Software Tools are accepted article
  types. Article-processing charge, data/code policy, and formatting:
  **[VERIFY]**. (For comparison, PLOS Computational Biology requires all
  study code to be public at publication; releasing this repository meets that
  bar.)

## What this repository supplies to meet those requirements
| Requirement | Where |
|---|---|
| Hypothesis-driven, pre-specified tests | `manuscript/` Introduction + `docs/02_study_design.md` |
| Reproducible methods | `scripts/`, `requirements.txt`, `data/` |
| Raw data and figures | `data/results/`, `figures/` |
| Data/code availability statement | Methods (repository URL) |
| JEI section order, <= 10 pp body, 110-char title | `manuscript/` |
