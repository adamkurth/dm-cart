# WRITE-UP — translating the thesis into `tex/article/main.tex`

Companion documents: [`notes.md`](notes.md) (your own interpretation of the code
revisions) and [`REVISIONS.md`](REVISIONS.md) (the full revision log, including its own
["Write-up guide"](REVISIONS.md#write-up-guide--what-to-say-where-and-how-to-restructure)
section, which this file assumes and does not repeat).

Those two files tell you **what changed in the code and why**. This file is scoped
differently: it is about turning `tex/thesis/` (four-chapter ASU dissertation, 2020/2021
cohort, ~10,200 words of chapter prose) into `tex/article/main.tex` (single-file,
twocolumn JCSA journal article). Two things have to happen at once, and they are easy to
tangle together if you don't separate them:

1. **STRIP** — cut content that a dissertation needs and a journal article does not.
   Nothing here is wrong, it just does not fit the format.
2. **RESPECIFY** — fix content that is factually wrong in the thesis, using the 2023/2024
   revision (`REVISIONS.md`). This is not optional shortening; the old numbers are
   incorrect and must not be copied forward even in a condensed form.

Do the STRIP pass and the RESPECIFY pass separately. Cutting a wrong sentence for length
and cutting a right sentence for length look the same in a diff, and it is easy to
"condense" a paragraph and accidentally preserve the error in fewer words. Fix first,
then trim.

---

## 1. What's wrong in the thesis (must not be copied forward, condensed or not)

These are corrections, not cuts. Every one of these appears in the thesis text as
currently written and is factually superseded. Full detail and exact numbers are in
[`REVISIONS.md` §"Write-up guide" §1–4](REVISIONS.md#1-the-three-claims-the-thesis-currently-makes-that-are-no-longer-true).

| # | Thesis says | Where (thesis) | Now |
|---|---|---|---|
| 1 | Cohort is 2020/2021 | `abstract.tex:2`, `chapter2.tex:33,35,40,76`, `chapter3.tex:8,15,19,43,61,156` (every `2020`/`2021`) | Cohort is **2023/2024**: cutpoints + α from 2023, counts from 2024. These were silently the *same* year before (R1) — do not restate the old claim even loosely |
| 2 | `I = 128` classes, `X ∈ {0,1}^{I×7}`, all 7 predictors binary | `chapter2.tex:40,71,76,129,153,255,269,274,277` | `I = 672`; `dmar` is 3-level, `race` is 7-level. Restate as `xᵢ ∈ ∏ⱼ 𝓛ⱼ`, `|𝓛| = (2,3,2,2,2,2,7)`, `I = ∏|𝓛ⱼ| = 672`. The DM likelihood and the ΔL splitting rule are **unchanged** — only the index set grows. Say so explicitly |
| 3 | Race binary split (`mrace15`) justified by national population share, flagged `entirely arbitrary. (NEED TO REVISIT)` | `chapter2.tex:45` | Full 7-level `mracehisp` (White, Black, Hispanic, Asian, AIAN, NHOPI, Multiracial). Do not amend the old paragraph — replace it. The new argument is methodological (sparse cells are what the informed Dirichlet prior is for), not demographic |
| 4 | `meduc == 3` labelled "HS completed" | Table 2.1, `chapter2.tex` predictor table | `meduc >= 3`. Old rule flagged 28.5%, correct rule flags 88.9%. A bachelor's-degree mother was being recorded as *not* having finished high school |
| 5 | `precare5` glossed as "adequate prenatal care" | Table 2.1 | Renamed `precare1st`, glossed as "prenatal care began in the 1st trimester." No numbers change, only the label — adequacy is the APNCU/Kotelchuck index, which is not in this data extract |
| 6 | High/low-risk subgroups defined by reading covariates (`i=69`/`i=28` in the old 128-class design) | `chapter3.tex` §3.3.5-equivalent | Subgroups are **found from the data**: rank all classes by bootstrap posterior risk, threshold at `N ≥ 1,000` births, report the extremes with tie structure. The assumed profiles are kept but reported as *where they landed* (rank 13 and rank 25-from-bottom), not as the headline result |
| 7 | Root split on race described as `mrace15 = 1` (Black) vs `mrace15 = 0` (everyone else) | `chapter3.tex:48` | Root split is a genuine multi-level partition the old binary encoding could not represent at all. Report the actual partition (see `REVISIONS.md` R15/§3 of the write-up guide) |

**Do not hand-fix these inline while translating.** Rewrite each paragraph from the
correct numbers in `REVISIONS.md`, the way §2–4 of its write-up guide already lay out
sentence-by-sentence. Treat the thesis prose as a *structure* to reuse, not text to
patch.

---

## 2. What to strip (format, not correctness)

The article template (`tex/article/main.tex`) is a single twocolumn, 10pt research
article. JCSA's own template caps the abstract at 150–250 words. There is no stated hard
page limit, but the structural signal is clear: this is a ~4,000–6,000 word paper, not a
~10,000-word, four-chapter document plus front matter. Cut in this order — each item
below is safe to drop with no loss of the argument the paper is actually making:

1. **All ASU dissertation scaffolding.** None of this has a place in a journal article —
   don't translate it, just don't carry it over:
   - `cover.tex`, `vita.tex`, `notation.tex`, `ack.tex`, the copyright page macro in
     `dis.tex`
   - List of Tables / List of Figures / Table of Contents machinery
   - `asudis.sty`, `tocloft`, `afterpage`, `pgffor` (the chapter-loop `\foreach`) — see the
     earlier package cleanup in `tex/article/main.tex`, these were deliberately not
     carried over
   - `appendix1.tex` in full — pull out only the two or three numbers a reviewer would
     actually ask for (missingness table, split-partition table) and inline them as a
     single compact table or drop them to supplementary material
2. **Abstract: 421 words → 150–250.** Current abstract (`abstract.tex`) spends sentences
   on background (LBW as a public-health indicator) that the Introduction also covers.
   Journal abstracts do not restate general background — state the method, the data, and
   the one or two headline numbers (the empirical high/low-risk contrast, the class
   count) and stop.
3. **Chapter 1 (Introduction, 1,137 words): compress to ~2–3 paragraphs.** The literature
   review of individual LBW risk factors (Kramer 1987, CART precedents, Bayesian
   nonparametric alternatives) can be one paragraph of citations, not four. Keep only
   what motivates *this* method choice: (a) LBW risk factors interact, not add — motivates
   tree-based partitioning; (b) existing CART-on-LBW work (Kitsantas 2006) has no
   uncertainty quantification — motivates the Bayesian/DM layer; (c) the framing sentence
   from `REVISIONS.md` §6: *a DM tree lets you disaggregate a categorical risk factor to a
   resolution a plain multinomial CART cannot support, and read risk profiles off the
   result.* That's the paper's contribution in one sentence — lead with it, don't bury it
   in paragraph 4.
4. **Chapter 2 (Methods, 4,548 words): keep the DM derivation, compress the exposition
   around it.** The derivation itself (Dirichlet-Multinomial conjugacy, the marginal
   likelihood, the ΔL splitting criterion) is the paper's methodological core and should
   stay close to full length — it's what makes this publishable as methods, not just an
   applied result. What compresses:
   - The general dataset description (NCHS/NBER background, all available variable
     domains) — two sentences, not a paragraph. The reader needs to know it's
     administrative birth-certificate data with ~3.6M annual records; the demographic/
     health/geographic variable taxonomy is not needed for the modeling section
   - Restate preprocessing as three tight steps (recoding → missingness → class
     aggregation) rather than a walkthrough
5. **Chapter 3 (Results, 3,571 words): cut the figure count roughly in half.** The thesis
   has depth-2/3/4/5 comparison trees as four separate full-page figures for *each* model
   arm (8 figures total) — for a journal article, one depth-comparison panel figure
   (small multiples) replaces all eight. Same for the tree renders: one annotated full
   tree per arm, not a full tree + a display tree + four depth variants. See the
   figure-mapping table in `REVISIONS.md` §4.5 for which output files are the ones worth
   keeping as figures versus which are diagnostic-only and belong in supplementary
   material at most.
6. **Chapter 4 (Discussion/Future Work, 517 words): fold into the article's Conclusion
   section, do not keep as a separate section.** The JCSA template has one `Conclusion and
   Future Work` section — this chapter's content (richer predictor sets, environmental/
   contextual variables, temporal trends) becomes the last paragraph, compressed to ~4–5
   sentences of concrete next steps, not a full section of its own.
7. **Citation style leftovers.** The thesis uses `natbib`-style `\citep`/`\citet`/
   `\citeyear` throughout; `tex/article/main.tex` now accepts these directly
   (`natbib=true` was added to the `biblatex` call), so citation commands do not need to
   be rewritten — but do check that every `\cite` key used in the translated text
   actually exists in `tex/article/references.bib` (the article's bib file is currently
   the JCSA sample file, not `tex/thesis/dis.bib` — the two bibliographies have not been
   merged yet, see open items below).

---

## 3. What to respecify (correctness fixes that change the shape of the argument, not just the numbers)

These aren't sentence-level corrections — they change what the results section is
*about*, so translate the restructured version, not the old section with new numbers
substituted in. `REVISIONS.md` §4 has the full sentence-level guide; this is the
shape-level summary so you don't paste old structure by habit.

- **The results chapter's center of gravity moves.** Old order: tree structure → depth
  comparison → variable importance → (assumed) subgroup contrast, as an afterthought. New
  order: tree structure → depth comparison → variable importance/stability →
  **risk-profile ranking** → extreme-ten and shared features. The ranking and the
  extreme-ten analysis are the new headline result, not a validation check tacked onto
  the end — restructure the section so the paper visibly builds toward it.
- **New subsection with no thesis precedent at all: the extreme ten and their shared
  features.** This is the strongest new result (`REVISIONS.md` §4.4) — unmarried status
  is the only unanimous marker among the ten highest-risk profiles (race 8/10, smoking
  6/10), while the ten lowest-risk profiles share three unanimous markers (male, white,
  non-smoker). That asymmetry — many routes to high risk, one narrow route to low risk —
  is a finding that literally does not exist in the thesis and needs its own paragraph,
  not a footnote.
- **New paragraph with no thesis precedent: the categorical split search itself is a
  methods contribution**, not an implementation detail (`REVISIONS.md` §3 of the write-up
  guide, and Chapter 2 §2.4.4a). The thesis has zero text on this. One paragraph: rpart's
  user-split protocol only reaches contiguous orderings; exhaustive search reaches all
  `2^(k-1)-1` bipartitions; report that this changed the fitted root split on real data,
  not just a synthetic demonstration.
- **`dmar`'s `Unknown` level is a methods contribution, argued, not mentioned.** The
  thesis doesn't have this predictor at all in its old 2-level form. Make the argument in
  full (geography is unrecoverable from the public file → imputation is not credible →
  explicit `Unknown` level assumes nothing → and it pays off empirically, the model groups
  `Unknown` with `married`) rather than a one-line footnote. `REVISIONS.md` §2.2 has the
  exact five-point argument to use.
- **A limitations paragraph the thesis does not have at all.** Four short points, given in
  full in `REVISIONS.md` §5 — `dmar` unknown-as-level, AIAN/NHOPI sparse-cell interval
  width, the `N ≥ 1,000` support threshold's effect on which classes qualify, and
  `precare1st` as timing-not-adequacy. A journal reviewer will ask about at least two of
  these if they are not pre-empted.

---

## 4. Section-by-section translation map

| Article section (`tex/article/main.tex`) | Thesis source | Strip | Respecify |
|---|---|---|---|
| Abstract | `abstract.tex` | Cut to 150–250 words; drop general LBW background | Cohort year, `I=672`, headline risk contrast numbers |
| Introduction | `chapter1/chapter1.tex` | Literature review to one paragraph | Lead with the "disaggregate a categorical risk factor" framing sentence |
| Literature Review | `chapter1/chapter1.tex` (back half) | Merge into Introduction if space is tight — JCSA template treats these as separate sections but a short paper can combine them | — |
| Methodology and Data | `chapter2/chapter2.tex` | Dataset-domain taxonomy, verbose preprocessing walkthrough | `I=672` notation, missingness subsection (new), `dmar` argument (new), race rewrite, Table 2.1 rebuild, split-search paragraph (new) |
| Results and Discussion | `chapter3/chapter3.tex` | Depth-comparison figures (8→1), duplicate tree renders | Cohort year throughout, subgroup section replaced wholesale, extreme-ten subsection (new), tie/stability caveats (new) |
| Conclusion and Future Work | `chapter4/chapter4.tex` | Compress to ~5 sentences | Add the limitations paragraph (new — was not a limitations section before, was future-work only) |
| Declaration | — (thesis has no equivalent) | — | Fill in conflict-of-interest / author-contribution placeholders |

---

## 5. Open items before this is submission-ready

- `tex/article/references.bib` is still the JCSA sample bibliography, not
  `tex/thesis/dis.bib`. Merge the citation keys you actually use once the text is
  translated — don't merge the whole file wholesale, the thesis bib has entries the
  condensed article won't cite.
- Every number in `REVISIONS.md` is from a `B = 150–200` smoke run, not the publication
  `B = 10,000` bootstrap. Do not typeset final tables/figures until `REVISIONS.md` §7's
  regeneration commands have been run:
  ```bash
  Rscript quantile.R 2023-2024
  Rscript data.R     2023-2024
  DMCART_B=10000 Rscript dm-cart.R 2023-2024
  Rscript visual.R   2023-2024
  Rscript visual.R   2023-2024 lbw
  ```
  The extreme classes and unanimity counts are stable at `B=200`; only interval endpoints
  move in the third decimal, but a reviewer can ask for the exact run.
- `visual.R` had a parse error as of the last `REVISIONS.md` pass (stray `»` character on
  line 18) — confirm it currently runs before pulling any figure from it.
- The `"binary"` and `"collapsed5"` race arms exist as one-line sensitivity checks but
  have not been run. Decide whether the article needs a sensitivity-analysis mention at
  all given the space budget — the thesis didn't have one, so this is a genuine addition,
  not a restoration.
