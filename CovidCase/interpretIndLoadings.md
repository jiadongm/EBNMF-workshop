# Interpretation of DIVAS Individual Component Loadings

> Source: `08c_divas_individual_loadings.R` and `results/divas_2block.rds`
> DIVAS 2-block proteomics analysis: T1 (admission) vs T2 (follow-up), 481 Olink proteins, 120 patients.
> **Centering:** `colCent=TRUE, rowCent=FALSE` (proteins centred within each block).

---

## Overview

DIVAS decomposed the two-time-point proteomics data into 24 components:
10 shared (T1+T2), 8 T1-individual, and 6 T2-individual. The **individual
components** capture proteomics variation unique to one time point -- axes
present at admission but not follow-up, or vice versa.

Unlike shared components (whose scores correlate strongly with severity --
e.g. T1+T2-1 rho=−0.721), individual components have uniformly weak
severity associations (all FDR > 0.3). This is expected: severity-related
variation is predominantly time-stable. Individual components instead
capture **time-specific biological processes** -- e.g. acute-phase responses
at admission that resolve by follow-up, or recovery/compensation programs
that emerge only later.

### Recurring proteins across individual components

Several proteins appear in the top 20 of 4+ components, indicating they are
major axes of individual variation:

| Protein | Appearances (of 14) | Biology |
|---------|---------------------|---------|
| Growth hormone 1 (GH1) | 8–9 | Somatotropic axis; stress/catabolism marker |
| Troponin I (cardiac) | 7–8 | Myocardial injury marker |
| IFN-gamma | 5–6 | Th1 cytokine; antiviral response |
| IL-10 (4 isoforms) | 5–6 | Anti-inflammatory cytokine |
| HAVCR1/TIM-1 (4 isoforms) | 5–6 | Kidney injury / hepatitis A receptor |
| SERPINA12 (vaspin) | 5 | Adipokine; insulin-sensitising |
| FGF21 (multiple isoforms) | 5 | Hepatokine; metabolic/mitochondrial stress |
| HAO1 | 4–5 | Peroxisomal oxidase; liver/metabolic |
| BNP (natriuretic peptide B) | 4 | Cardiac stress marker |
| NCF2/p67phox | 4 | NADPH oxidase component; neutrophil ROS |
| IL-5 (4 isoforms) | 4–5 | Eosinophil-activating cytokine |
| FBP1 | 4 | Gluconeogenesis enzyme |
| CA5A | 4 | Mitochondrial carbonic anhydrase; metabolic |

These proteins define recurring biological themes: cardiac/renal stress,
hepatic/metabolic function, Th1/Th2 immune balance, and neutrophil activation.

---

## T1-Individual Components (admission-specific)

### T1-Individual-1: Troponin I / IFN-gamma vs. IL-10 axis

**Top positive:** TnI (+0.188), parvalbumin, FKBP4, CGRP, FGF21

**Top negative:** IFN-gamma (−0.178), FABP2, HAO1, IL-10

**Interpretation:** A bipolar axis at admission separating a cardiac stress /
immunophilin profile (TnI, FKBP4 -- involved in steroid receptor signalling)
from an antiviral/anti-inflammatory response (IFN-gamma, IL-10). FABP2
(intestinal fatty acid binding) and HAO1 (hepatic) on the negative pole suggest
gut/hepatic involvement in the IFN-gamma response. CGRP (calcitonin gene-related
peptide, a neuropeptide and vasodilator) co-loading with TnI may reflect
neurogenic cardiac involvement. This axis is admission-specific -- by follow-up,
the IFN-gamma response and cardiac markers diverge onto different axes.

### T1-Individual-2: Troponin I / HAVCR1 cardiac-renal damage

**Top positive:** TnI (+0.262), HAVCR1/TIM-1 × 4 isoforms, meprin A

**Top negative:** IL-5 × 2 isoforms (−0.178), GH1, FGF21 × 3 isoforms

**Interpretation:** Dominated by organ damage markers: cardiac troponin I
(myocardial injury) and HAVCR1/TIM-1 (kidney injury molecule, also the
hepatitis A virus receptor). Meprin A (metalloprotease) adds active tissue
remodelling. The negative pole (IL-5, GH1, FGF21) represents a cytokine/
metabolic stress axis. This captures patients at admission with acute cardiac
and renal damage. The four HAVCR1 isoforms loading together indicate that
this is a genuine kidney injury signal, not noise.

### T1-Individual-3: HAO1 / FGF21 hepatic metabolic stress

**Top positive:** RNase 3/ECP (+0.147), GH1, IL-27, LRIG1, CGRP, PTX3

**Top negative:** HAO1 (−0.313), CA5A, FGF21 × 4 isoforms

**Interpretation:** The strongest loading is HAO1 (liver peroxisomal
hydroxyacid oxidase) in the negative direction, co-loading with CA5A
(mitochondrial carbonic anhydrase) and FGF21 (hepatokine induced by
metabolic stress and fasting). The negative pole is a concentrated
hepatic/mitochondrial stress signature. PTX3 (pentraxin 3, acute-phase
protein) and RNase 3 on the positive pole represent innate immune activation
opposed by hepatic dysfunction. Distinct from T1-Indiv-2 because it
emphasises metabolic/mitochondrial liver stress rather than cardiac/renal
injury.

### T1-Individual-4: NAD kinase / TRIM21 vs. HAVCR1 axis

**Top positive:** NAD kinase (+0.129), TRIM21, thymidine phosphorylase,
IFN-gamma

**Top negative:** HAVCR1 × 2 isoforms (−0.133), FGF21 × 2 isoforms

**Interpretation:** NAD kinase (metabolic enzyme for NAD+ homeostasis) and
TRIM21 (E3 ubiquitin ligase; innate immune sensor for intracellular antibody-
coated pathogens) on the positive pole suggest a coordinated metabolic and
innate immune sensing response. IFN-gamma co-loading on the positive side
links this to antiviral immunity. The negative pole (HAVCR1, FGF21) represents
kidney injury / hepatic stress. This axis captures patients with active
innate immune metabolic reprogramming at admission, distinct from the
hepatic/cardiac patterns of other components.

### T1-Individual-5: Meprin A / BNP vs. IFN-gamma axis

**Top positive:** Meprin A (+0.165), TnI, CRH (corticotropin-releasing
hormone), BNP, renin

**Top negative:** IFN-gamma (−0.269), HAO1, CA5A, CXCL9

**Interpretation:** The positive pole is an acute cardiovascular/
neuroendocrine stress axis: meprin A (tissue-remodelling metalloprotease),
troponin I, BNP (heart failure marker), renin (renin-angiotensin system),
and CRH (HPA axis activation). The negative pole is a coordinated IFN-gamma
signalling module (IFN-gamma, CXCL9 -- IFN-gamma-induced chemokine). This
is one of the clearest Th1 antiviral vs. cardiovascular/endocrine stress
axes. At follow-up, this specific pairing resolves: the antiviral response
and cardiovascular stress either normalise or separate onto shared axes.

### T1-Individual-6: SERPINA12 / RNase 3 / troponin cardiovascular axis

**Top positive:** SERPINA12/vaspin (+0.221), RNase 3/ECP, TnI, NCF2,
CXCL5, BNP

**Top negative:** KRT19 (−0.162), AGR2, GH1, CA5A

**Interpretation:** The positive pole combines cardiovascular markers (TnI,
BNP, vaspin -- an adipokine linked to insulin resistance and cardiovascular
risk) with innate immune markers (RNase 3/ECP from eosinophils, NCF2/p67phox
from neutrophils, CXCL5 -- neutrophil chemokine). This co-loading suggests
concurrent cardiac injury and granulocyte activation at admission in a
patient subgroup. KRT19 (keratin-19, epithelial marker) and AGR2 (anterior
gradient 2, mucus-secreting cells) on the negative pole suggest mucosal/
epithelial compartment biology. This axis is admission-specific.

### T1-Individual-7: IL-10 / pleiotrophin vs. GH1 axis

**Top positive:** IL-10 × 4 isoforms (+0.166), pleiotrophin, GPR56

**Top negative:** GH1 (−0.245), TnI, SERPINA12, BNP, parvalbumin

**Interpretation:** IL-10 (all four isoforms loading together) dominates the
positive pole, indicating a robust anti-inflammatory signal. Pleiotrophin
(heparin-binding growth factor with anti-inflammatory effects) and GPR56
(adhesion GPCR expressed on NK cells and microglia) support an
immunomodulatory program. The negative pole is a cardiac/somatotropic cluster
(GH1, TnI, BNP, SERPINA12). This axis captures patients at admission where
higher anti-inflammatory (IL-10) signalling is opposed by cardiac injury and
metabolic stress markers. By follow-up, this particular pairing is absent.

### T1-Individual-8: NCF2 / nibrin neutrophil oxidative stress axis

**Top positive:** NCF2/p67phox (+0.191), PCBP2, nibrin/NBN, SH2D1A,
paxillin

**Top negative:** SERPINA12 (−0.209), TnI, BNP, IFN-gamma, FBP1

**Interpretation:** The positive pole is a DNA damage / oxidative burst
module: NCF2/p67phox (NADPH oxidase subunit -- essential for the neutrophil
respiratory burst that generates ROS), nibrin/NBN (a key DNA double-strand
break repair sensor -- activated by ROS-induced genotoxic stress), and
PCBP2 (RNA-binding protein involved in iron homeostasis and viral response).
SH2D1A (SAP -- required for NKT cell and CD8 T cell function) and paxillin
(cell adhesion/cytoskeletal) suggest immune cell migration and degranulation.
The negative pole (SERPINA12, TnI, BNP, FBP1) represents cardiovascular/
metabolic markers. This admission-specific axis captures patients with active
neutrophil-driven oxidative stress, potentially causing secondary genotoxic
damage.

---

## T2-Individual Components (follow-up-specific)

### T2-Individual-1: FGF21 / CCL3 vs. HAO1 / ARNT metabolic axis

**Top positive:** FGF21 × 4 isoforms (+0.152), CCL3 × 2 isoforms

**Top negative:** ARNT (−0.150, hypoxia-inducible factor 1β), HAO1, KPNA1,
RASA1, meprin A

**Interpretation:** FGF21 (hepatic stress / fasting response) dominates the
positive pole at follow-up, opposed by ARNT (the obligate partner for HIF-1α
in the hypoxia response) and HAO1 (hepatic peroxisomal). CCL3/MIP-1α (T cell
and macrophage chemokine) co-loading with FGF21 links hepatic metabolic stress
to immune cell recruitment. RASA1 (RAS GTPase activating protein; negative
regulator of RAS-MAPK signalling) on the negative pole suggests active growth
factor signalling during recovery. This is a follow-up-specific hepatic
metabolism axis: FGF21 elevation at T2 (metabolic recovery) vs. ARNT/HAO1
(residual hypoxic/peroxisomal stress).

### T2-Individual-2: IL-6 vs. troponin I / HAVCR1 inflammation axis

**Top positive:** IL-6 × 6 isoforms (+0.129)

**Top negative:** TnI (−0.245), HAVCR1 × 4 isoforms, CGRP

**Interpretation:** IL-6 (all six panel isoforms loading together as a
coherent signal) defines the positive pole: persistent IL-6 inflammation at
follow-up. This is opposed by troponin I (cardiac injury) and all four
HAVCR1 isoforms (kidney injury). The negative-positive pairing indicates
that IL-6-high patients at T2 have low cardiac/renal damage, and vice versa.
CGRP on the negative side suggests neurogenic vasodilation in patients with
cardiac injury. This axis is T2-specific: IL-6 variation at admission is
captured differently (through shared components).

### T2-Individual-3: GH1 / IL-5 vs. HAO1 / meprin hepatorenal axis

**Top positive:** GH1 (+0.266), IL-5 × 4 isoforms, CD5

**Top negative:** HAO1 (−0.236), meprin A, FBP1, HNMT, TMPRSS15, CA5A

**Interpretation:** GH1 dominates the positive pole with IL-5 isoforms
(eosinophil activation), suggesting somatotropic recovery with concurrent
eosinophil-mediated immunity at follow-up. CD5 (T cell / B-1 cell marker)
adds an adaptive immune component. The negative pole is a concentrated
hepatorenal damage signature: HAO1 (liver peroxisomal), FBP1 (gluconeogenesis),
CA5A (mitochondrial), HNMT (histamine N-methyltransferase, liver enzyme),
meprin A (metalloprotease), and TMPRSS15 (enterokinase, GI). This mirrors
T1-Individual-3 but at T2 -- persistent organ damage remaining as a distinct
axis at follow-up.

### T2-Individual-4: GH1 / IL-10 vs. HAVCR1 / SERPINA12 axis

**Top positive:** GH1 (+0.156), DDX58/RIG-I, IL-10 × 2 isoforms, PCBP2

**Top negative:** HAVCR1 × 4 isoforms (−0.175), SERPINA12, CXCL5

**Interpretation:** GH1 and IL-10 (anti-inflammatory) co-load with DDX58/
RIG-I (cytosolic viral RNA sensor -- innate immune activation against
remaining viral RNA at follow-up) and PCBP2 (iron/viral response regulator).
The positive pole thus captures a mixed recovery/innate-antiviral program.
The negative pole (HAVCR1 isoforms, SERPINA12, CXCL5) represents ongoing
kidney injury and adipokine/neutrophil signalling. This axis separates patients
with active innate antiviral + metabolic recovery at T2 from those with
persistent kidney injury.

### T2-Individual-5: GH1 / HAVCR1 vs. TMPRSS15 / EPO axis

**Top positive:** GH1 (+0.305), PAPPA (pregnancy-associated plasma protein
A; metalloprotease; tissue remodelling), HAVCR1 × 4 isoforms, RNase 3/ECP

**Top negative:** TMPRSS15 (−0.136, enterokinase/GI marker), EPO
(erythropoietin), IL-5 × 2 isoforms, PKC theta

**Interpretation:** GH1 and PAPPA on the positive pole suggest active tissue
remodelling / growth factor signalling at follow-up. HAVCR1 isoforms and
RNase 3/ECP indicate concurrent kidney injury and eosinophil activation.
The negative pole features TMPRSS15 (GI serine protease), EPO (erythropoietin
-- kidney-produced; depletion may indicate renal failure), and IL-5 (eosinophil
growth factor). This follow-up-specific axis captures patients with active
tissue remodelling (GH1, PAPPA) and kidney involvement (HAVCR1, EPO depletion).

### T2-Individual-6: GH1 / SERPINA12 vs. IL-10 / CCL20 axis

**Top positive:** GH1 (+0.294), SERPINA12/vaspin, FBP1, WAS (Wiskott-Aldrich
syndrome protein; cytoskeletal regulation in immune cells), pleiotrophin

**Top negative:** IL-10 × 4 isoforms (−0.194), CCL20, TMPRSS15

**Interpretation:** GH1 and SERPINA12 (adipokine / insulin-sensitising)
on the positive pole, combined with FBP1 (gluconeogenesis) and pleiotrophin
(growth factor), represent a metabolic/endocrine recovery axis. WAS protein
(required for actin polymerisation in immune cells; regulates T cell and NK
cell function) co-loading suggests cytoskeletal immune activity alongside
metabolic recovery. The negative pole (IL-10 isoforms and CCL20) captures
anti-inflammatory and mucosal immune signalling. This mirrors T1-Indiv-7
(IL-10 / GH1 pairing) but at T2 and with reversed polarity: at T2, the
metabolic recovery (GH1) and anti-inflammatory (IL-10) programs are more
clearly separated into individual axes.

---

## Cross-Component Themes

### 1. Troponin I and cardiac markers are prominent at admission

TnI appears in the top 20 of 5–6 T1-individual components, BNP in 3–4.
T1-Indiv-2 (troponin + HAVCR1) and T1-Indiv-5 (troponin + BNP + meprin)
capture distinct cardiac injury patterns at admission. These are largely
absent as individual axes at T2, suggesting cardiac injury either resolves
or merges into shared severity components by follow-up.

### 2. GH1 is ubiquitous at follow-up

Growth hormone 1 appears in 5–6 of the 6 T2-individual components, always
as a major loading. GH1 is an acute-phase protein whose direction reflects
different biological contexts (catabolic stress, pituitary suppression by
cytokines, or somatotropic recovery). Its dominance at T2 suggests that
somatotropic axis variation is a major source of inter-patient heterogeneity
during recovery, while at admission, cardiac and immune markers account for
more variation.

### 3. HAVCR1/TIM-1 isoforms cluster together

The four HAVCR1 isoforms consistently load together across multiple
components (T1-Indiv-2, T1-Indiv-4, T2-Indiv-2, T2-Indiv-4, T2-Indiv-5).
This coherent co-loading confirms HAVCR1 as a genuine biological signal
(kidney injury / hepatitis A receptor) rather than assay noise.

### 4. IFN-gamma and IL-10 reflect admission-specific immune state

IFN-gamma is most prominent in T1-individual components (T1-Indiv-1, T1-4,
T1-5), consistent with the acute antiviral response at admission. IL-10
spans both admission and follow-up individual components (T1-Indiv-7,
T2-Indiv-4, T2-Indiv-6), reflecting its role as a persistent regulatory
brake that tracks with different opposing signals at each time point.

### 5. FGF21 / HAO1 hepatic stress is admission-prominent

FGF21 (multiple isoforms) and HAO1 are prominent in T1-individual
components (T1-2, T1-3, T1-4) reflecting hepatic metabolic derangement at
admission. At T2, FGF21 re-emerges in T2-Indiv-1 (opposing ARNT/HAO1)
but with lower magnitude, suggesting partial hepatic recovery.

---

## Relationship to Shared Components

The individual components capture variation that is **not** present in the
10 shared components (T1+T2-1 through T1+T2-10). By definition, individual
variation at one time point is orthogonal to shared variation. Key contrasts:

- **T1+T2-1** (rho=−0.721 with severity) captures time-stable severity --
  the same protein axis distinguishes mild/severe at both time points.
  Individual components capture what is **different** about admission vs
  follow-up, not what is the same.

- **T1+T2-6** (rho=−0.334 with delta) captures shared proteomics change
  associated with improvement. Individual components have no equivalent
  signal (all FDR > 0.3), confirming improvement-associated variation is
  also time-stable.

- The absence of strong severity or delta associations in individual
  components (all FDR > 0.3) confirms that COVID-19 severity and clinical
  trajectory are predominantly **time-stable** proteomics phenomena in
  this cohort.

---

## Summary Table

| Component | Dominant theme | Key proteins | Direction |
|-----------|---------------|--------------|-----------|
| T1-Indiv-1 | TnI / IFN-gamma vs. IL-10 | TnI, FKBP4, CGRP vs. IFN-γ, IL-10 | Positive = cardiac/immunophilin |
| T1-Indiv-2 | Troponin / kidney damage | TnI, HAVCR1 × 4 vs. IL-5, GH1 | Positive = organ damage |
| T1-Indiv-3 | Hepatic metabolic stress | HAO1, CA5A, FGF21 vs. RNase 3, GH1 | Negative = hepatic stress |
| T1-Indiv-4 | NAD kinase / innate immunity | NAD kinase, TRIM21, IFN-γ vs. HAVCR1, FGF21 | Positive = metabolic-immune |
| T1-Indiv-5 | Meprin / cardiovascular vs. IFN-γ | Meprin, TnI, BNP, renin vs. IFN-γ, CXCL9 | Positive = cardiovascular |
| T1-Indiv-6 | SERPINA12 / cardiac-innate | SERPINA12, RNase 3, TnI, NCF2 vs. KRT19, GH1 | Positive = cardiac-innate |
| T1-Indiv-7 | IL-10 / pleiotrophin vs. GH1 | IL-10 × 4 vs. GH1, TnI, SERPINA12 | Positive = anti-inflammatory |
| T1-Indiv-8 | NCF2 / nibrin neutrophil ROS | NCF2, nibrin, SH2D1A vs. SERPINA12, TnI | Positive = neutrophil oxidative |
| T2-Indiv-1 | FGF21 / CCL3 vs. ARNT / HAO1 | FGF21 × 4, CCL3 vs. ARNT, HAO1 | Positive = FGF21/metabolic |
| T2-Indiv-2 | IL-6 vs. troponin / HAVCR1 | IL-6 × 6 vs. TnI, HAVCR1 × 4 | Positive = IL-6 persistence |
| T2-Indiv-3 | GH1 / IL-5 vs. hepatorenal | GH1, IL-5 × 4 vs. HAO1, FBP1, CA5A | Positive = somatotropic |
| T2-Indiv-4 | GH1 / RIG-I vs. HAVCR1 | GH1, DDX58, IL-10 vs. HAVCR1 × 4 | Positive = recovery-antiviral |
| T2-Indiv-5 | GH1 / HAVCR1 vs. TMPRSS15 / EPO | GH1, PAPPA, HAVCR1 vs. TMPRSS15, EPO | Positive = GH1/remodelling |
| T2-Indiv-6 | GH1 / SERPINA12 vs. IL-10 | GH1, SERPINA12, FBP1 vs. IL-10 × 4 | Positive = metabolic recovery |
