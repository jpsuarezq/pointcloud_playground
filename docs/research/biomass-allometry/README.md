# Above-ground biomass (AGB) allometry for canopy trees in a cacao agroforestry system

**Site:** cacao (*Theobroma cacao* L.) agroforestry system, Santander, Colombia — approx. 6.86 °N, −73.41 °W, ~900 m a.s.l.
**Available per-tree inputs:** DBH (diameter at breast height, 1.3 m) in **cm**; total height **H** in **m**.
**Not available:** species identity; any field (destructive) measurements; wood density measurements.

This document is written to be used directly as a **Methods citation**: each equation is given with its exact
form, parameter values, units, validity range, and primary source. It ends with a single **recommended default**
(equation + wood specific gravity + carbon fraction) and a DBH-only fallback.

---

## 0. Executive summary

| Step | Recommended choice | Formula | Source |
|---|---|---|---|
| AGB (primary, DBH + H) | **Chave et al. 2014, Eq. 4** | `AGB = 0.0673 × (ρ·D²·H)^0.976` | Chave et al. 2014, *Glob. Change Biol.* |
| Wood specific gravity (ρ) | **0.60 g cm⁻³** (generic tropical default); regional alt. 0.63 | — | IPCC/FAO Tier-1 practice; Zanne et al. 2009 |
| AGB (DBH-only fallback) | **Chave et al. 2014, Eq. 7** (corrected by Réjou-Méchain et al. 2017) | `AGB = exp(−2.024 − 0.896E + 0.920 ln ρ + 2.795 ln D − 0.0461 (ln D)²)` | Chave et al. 2014; Réjou-Méchain et al. 2017 |
| Carbon fraction (CF) | **0.47** t C (t d.m.)⁻¹ | `C = AGB × 0.47` | IPCC 2006 GL, Vol. 4, Ch. 4, Table 4.3 |
| CO₂ equivalent | **44/12 = 3.6667** | `CO₂e = C × 44/12` | IPCC; molar-mass ratio |

Net result for one tree: `CO₂e (kg) = AGB × 0.47 × 3.6667 = AGB × 1.7233`.

---

## 1. Primary equation — pantropical, uses DBH **and** height (Chave et al. 2014)

### 1.1 The equation (Chave et al. 2014, Eq. 4)

```
AGB = 0.0673 × (ρ × D² × H)^0.976
```

| Symbol | Meaning | Unit |
|---|---|---|
| `AGB` | above-ground biomass, **oven-dry** | **kg** |
| `D` | DBH (diameter at breast height, 1.3 m) | **cm** |
| `H` | total tree height | **m** |
| `ρ` (rho, wood specific gravity) | oven-dry mass / green volume | **g cm⁻³** |

- **Fitted on:** 4,004 destructively harvested trees ≥ 5 cm DBH, from 58 sites across the tropics.
- **Fit statistics:** residual standard error RSE = 0.357 (in ln units), df = 4002; the authors report a
  mean tree-level CV of **56.5 %** and a mean bias of **5.31 %**.
- **Key property:** when D, H and ρ are all included, *a single model holds across tropical vegetation types*
  — no detectable effect of region or environment was found. This is why it is the standard "generic" choice.
- **Validity:** trees ≥ 5 cm DBH. Do not extrapolate to very large stems far beyond the calibration range.

> Chave J, Réjou-Méchain M, Búrquez A, et al. (2014). *Improved allometric models to estimate the aboveground
> biomass of tropical trees.* **Global Change Biology 20(10): 3177–3190.** doi:10.1111/gcb.12629
> URL: https://onlinelibrary.wiley.com/doi/10.1111/gcb.12629 · PDF: https://jeromechave.github.io/pubs/chave-gcb14.pdf

### 1.2 Chave et al. 2005 (Oecologia) — the older pantropical equations

Chave et al. (2005) fitted forest-type-specific equations of the form `ln(AGB) = a + b·ln(ρ·D²·H)`.
They remain widely used and are the fallback when the environmental stress variable `E` (see §2) is not
available. The canonical coefficients (D in cm, H in m, ρ in g cm⁻³, AGB in kg):

| Forest type | Formula | Equivalent power form |
|---|---|---|
| **Moist** (1–4 dry months; ~1500–3500 mm yr⁻¹) | `AGB = exp(−2.977 + ln(ρ·D²·H))` | `AGB = 0.0509 × ρ·D²·H` |
| **Wet** (>3500 mm yr⁻¹, no seasonality; cloud forest) | `AGB = exp(−2.557 + 0.940·ln(ρ·D²·H))` | `AGB = 0.0776 × (ρ·D²·H)^0.940` |
| **Dry** (<1500 mm yr⁻¹, >5 dry months) | `AGB = exp(−2.187 + 0.916·ln(ρ·D²·H))` | `AGB = 0.112 × (ρ·D²·H)^0.916` |
| **All types** (model I.4) | `AGB = exp(−2.922 + 0.990·ln(ρ·D²·H))` | `AGB = 0.0538 × (ρ·D²·H)^0.990` |

- **Validity:** D from **5 to 156 cm**, and ρ·D²·H from 50 to 1,000,000. Broadleaf trees only
  (different models are needed for conifers, palms and lianas).
- A cacao agroforest at ~900 m in Santander (perhumid/moist, Andean foothill) is best matched by the
  **moist** equation; the 2014 Eq. 4 above is nevertheless preferred because it no longer requires choosing
  a forest type.

> Chave J, Andalo C, Brown S, et al. (2005). *Tree allometry and improved estimation of carbon stocks and
> balance in tropical forests.* **Oecologia 145(1): 87–99.** doi:10.1007/s00442-005-0100-x
> URL: https://link.springer.com/article/10.1007/s00442-005-0100-x

**Agreement check:** for a tree with D = 30 cm, H = 25 m, ρ = 0.6 g cm⁻³, Chave 2014 Eq. 4 gives 723 kg,
Chave 2005 moist gives 687 kg — a ~5 % difference, i.e. the two are practically interchangeable at this size.

---

## 2. DBH-only fallback (when H is missing) — Chave et al. 2014, Eq. 7

When total height is unavailable for a tree, AGB can still be inferred from DBH and ρ, but the model must
absorb the (unknown) diameter–height relationship via the **environmental stress variable `E`**:

### 2.1 As originally published (Chave et al. 2014, Eq. 7)

```
AGB = exp[ −1.803 − 0.976·E + 0.976·ln(ρ) + 2.673·ln(D) − 0.0299·(ln D)² ]
```

with

```
E = ( 0.178·TS − 0.938·CWD − 6.61·PS ) × 10⁻³
```

where `TS` = temperature seasonality (WorldClim bioclim 4), `CWD` = climatic water deficit (mm yr⁻¹),
`PS` = precipitation seasonality (WorldClim bioclim 15). `E` is a dimensionless index of environmental
stress (higher E = more stress = shorter trees for a given diameter).

- **Fit statistics:** RSE = 0.431; mean CV = **71.5 %**; mean bias = **9.71 %** — clearly worse than Eq. 4.
- Chave et al. explicitly recommend developing **locally derived H–D relationships** rather than relying on
  this model whenever possible.

> Same source as §1.1. A ready-to-use gridded raster of `E` (2.5 arc-min) is distributed at
> http://chave.ups-tlse.fr/pantropical_allometry.htm (mirror: https://jeromechave.github.io/pantropical_allometry.htm).

### 2.2 Corrected version (use this one if computing E yourself) — Réjou-Méchain et al. 2017

The authors of the `BIOMASS` R package found that Eq. 7 in the 2014 paper was **not** fitted directly to the
data (it was obtained by combining Eq. 4 and the H–D equation 6) and that it omits the Baskerville (1972)
log-normal bias correction. The directly fitted, corrected equation is:

```
AGB = exp( −2.024 − 0.896·E + 0.920·ln(ρ) + 2.795·ln(D) − 0.0461·(ln D)² )
```

This is the equation implemented in `BIOMASS::computeAGB()` when site coordinates are supplied.

> Réjou-Méchain M, Tanguy A, Piponiot C, Chave J, Hérault B (2017). *biomass: an R package for estimating
> above-ground biomass and its uncertainty in tropical forests.* **Methods in Ecology and Evolution 8(9):
> 1163–1167.** doi:10.1111/2041-210X.12753 · R docs: https://search.r-project.org/CRAN/refmans/BIOMASS/html/computeAGB.html

### 2.3 Chave et al. 2005 DBH-only model (with ρ) — use with caution in Colombia

```
AGB = ρ × exp( −1.499 + 2.148·ln(D) + 0.207·(ln D)² − 0.0281·(ln D)³ )
```

(D in cm, ρ in g cm⁻³, AGB in kg.) This is the widely used "Chave 2005 moist, no height" model.
**Important caveat:** Álvarez et al. (2012) tested the Chave (2005) Type-II (no-height) models on 631
destructively harvested trees across Colombia and found they **overestimated AGB by 54.7 ± 135.7 %**.
For Colombian conditions prefer §2.1/§2.2 or the Colombia-specific §6 model, and treat this one as a last
resort. (Implementation reference: ForestPlots `BiomasaFP` `AGBChv05M`.)

---

## 3. Wood specific gravity (ρ) when species is unknown

ρ enters the primary equation as `ρ^0.976` (Chave 2014 Eq. 4), i.e. AGB scales **almost linearly** with ρ
(a ±10 % error in ρ produces a ±~10 % error in AGB). It is therefore one of the most influential inputs.

### 3.1 Recommended default

| Value | Basis | Source |
|---|---|---|
| **0.60 g cm⁻³** | Generic tropical Tier-1 default widely applied in carbon accounting when taxa are unknown; also the rounded stand-level mean used in most pantropical/cacao studies | IPCC 2006 GL practice (Tier-1 default logic); Winrock Allometric Equation Guidance (2016) |
| 0.632 ± 0.178 g cm⁻³ | Regional mean, **tropical South America** (n = 4,191 taxa) | Zanne et al. 2009, *Ecol. Lett.* 12:351–366, Table 2 (Global Wood Density Database) |
| 0.645 g cm⁻³ | Neotropical mean (2,456 species) | Chave et al. 2006, *Ecol. Appl.* 16:2356–2367 |

**Recommendation:** use **ρ = 0.60 g cm⁻³** as the single defensible default for this site; if a
region-specific mean is preferred, **0.63 g cm⁻³** (tropical South America, Zanne et al. 2009) is the better
justified number and raises AGB by only ~5 %. A locally measured stand mean is always preferable if it can
be obtained.

### 3.2 Sensitivity and rules of thumb

- Chave et al. (2005) advise: if no species-level information exists, use a **stand-level average** of ρ;
  if even that is unavailable, use a regional/global database mean.
- Because the exponent is 0.976 (Chave 2014), **ΔAGB/AGB ≈ 0.976 × Δρ/ρ**. Example: switching from 0.60 to
  0.65 g cm⁻³ (+8.3 %) raises AGB by ~8 %.
- FAO/IPCC give the uncertainty of basic wood density as **±10 % to ±40 %** (IPCC 2006 GL, Vol. 4, Ch. 4,
  uncertainty section) — this dominates the error budget when species are unknown.
- The **Global Wood Density Database** (GWDD; Zanne et al. 2009, Dryad doi:10.5061/dryad.234; v2 Fischer et al.
  2026, doi:10.5281/zenodo.16919509) is the primary source to look up **genus- or family-level means** if any
  taxonomic information becomes available (genus means explain ~74 % of species-level variance; Chave et al.
  2006). Any lookup requires identifying the taxa, which is not currently possible here.

> Zanne AE, Lopez-Gonzalez G, Coomes DA, et al. (2009). *Towards a worldwide wood economics spectrum.*
> **Ecology Letters 12(4): 351–366.** doi:10.1111/j.1461-0248.2009.01285.x · Data: https://doi.org/10.5061/dryad.234
> Chave J, Muller-Landau HC, Baker TR, et al. (2006). *Regional and phylogenetic variation of wood density
> across 2456 Neotropical tree species.* **Ecological Applications 16(6): 2356–2367.**
> doi:10.1890/1051-0761(2006)016[2356:RAPVOW]2.0.CO;2

---

## 4. Biomass → carbon → CO₂ equivalent

### 4.1 Carbon fraction

```
C = AGB × CF ,   CF = 0.47 t C (t d.m.)⁻¹
```

The **IPCC default carbon fraction for above-ground forest biomass is 0.47** (2006 IPCC Guidelines for
National Greenhouse Gas Inventories, Vol. 4 [AFOLU], Ch. 4 [Forest Land], **Table 4.3**). The same 0.47 factor
is retained in the 2019 Refinement. (Older guidance and some pools use 0.5; the current IPCC default for
forest woody biomass is 0.47.)

> IPCC (2006). *2006 IPCC Guidelines for National Greenhouse Gas Inventories, Vol. 4: Agriculture, Forestry
> and Other Land Use, Ch. 4: Forest Land.* Table 4.3.
> URL: https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/4_Volume4/V4_04_Ch4_Forest_Land.pdf

### 4.2 CO₂ equivalent

```
CO₂e = C × (44/12) = C × 3.6667
```

44/12 is the molar-mass ratio of CO₂ to C (44 g mol⁻¹ / 12 g mol⁻¹), the standard IPCC/UNFCCC conversion.
Combining with CF:

```
CO₂e = AGB × 0.47 × 3.6667 = AGB × 1.7233
```

(A Colombian cacao-agroforest study applied exactly these two factors — CF = 0.47 and the 3.67 stoichiometric
CO₂:C ratio — to the same kind of inventory: Díaz-Cháux et al. 2025, doi:10.1371/journal.pone.0325477.)

> IPCC (2006), Vol. 4, Ch. 4, and GPG-LULUCF: "To convert tonnes C to Gg CO₂, multiply the value by 44/12".
> URL: https://www.ipcc-nggip.iges.or.jp/public/gpglulucf/gpglulucf_files/Chp3/Chp3_2_Forest_Land.pdf

---

## 5. Cacao / *Theobroma cacao* and cacao-agroforestry-specific equations

These equations estimate the biomass of the **cacao shrub itself** (or specific companion species), **not** the
generic canopy/shade trees. They are provided because you asked, but they should be applied **only** to cacao
plants (and identified species), never to unidentified canopy trees — for those use the pantropical equation
in §1. Cacao is normally measured as **d30** (stem diameter at 30 cm above ground) rather than DBH.

| Target | Formula | Units | Source |
|---|---|---|---|
| *T. cacao* (Costa Rica) | `AGB = 10^(−1.625 + 2.63·log₁₀(d30))`  ≡ `AGB = 0.0237 × d30^2.63` | AGB kg; d30 cm | Andrade et al. (2008), as compiled in Cardozo et al. (2022) and Díaz-Cháux et al. (2025) |
| *T. cacao* (Cameroon) | `ln(AGB) = −1.613 + 1.83·ln(D)`  (site-species model M1) | AGB kg; D = DBH cm; n = 50, 5–27 cm | Fomekong Tane et al. (2024) |
| *T. cacao* (Ghana) | `AGB = 0.7217·ρ·(D²)^0.921` | AGB kg; D cm; ρ g cm⁻³ | Donkor et al. (2023), as compiled in the Côte d'Ivoire comparison (2024) |
| *T. cacao* (Bolivia/Costa Rica) | `AGB = −0.0376 + 0.133·BA` | AGB kg; BA = basal area cm² | Beer et al., as applied in Niether et al. 2016 (PLOS ONE) |

**Key findings relevant to this project:**

- In cacao agroforests, **the cacao shrubs themselves contribute only a small share of stand C** (e.g.
  ~2–12 % of live-tree C in Cameroon, Saj et al. 2013; the canopy/shade trees dominate). So the choice of
  *tree* equation matters far more than the cacao equation.
- Cacao-agroforest studies almost universally estimate **shade trees with pantropical equations** — most
  commonly **Chave et al. 2005 moist** (e.g. Cameroon: Saj et al. 2013; Kotto-Same; Indonesia) or
  **Chave et al. 2014** (e.g. Niether et al. 2016). A recent Colombian Amazon cacao-agroforestry study
  (Díaz-Cháux et al. 2025) used the **Colombia-specific Álvarez et al. (2012)** multispecies tree equation
  for its trees and the Andrade cacao equation for the shrubs — see §6.
- Species-specific shade-tree allometries exist (e.g. Saj et al. 2022, Talamanca, Costa Rica,
  doi:10.1007/s10457-022-00741-y; Mensah et al. 2023, Ghana, doi:10.1038/s41598-023-42219-6) but all require
  species identity, which you do not have. **A generic pantropical equation is acceptable for an
  unidentified shade-tree inventory** and is the current standard practice; its caveats are in §8.

> Fomekong Tane AB, et al. (2024). *Site-Species Allometry Equation for Theobroma cacao L. Biomass Estimation
> in Agroforestry Systems of Cameroon.* **International Journal of Ecology 2024: 7096854.**
> doi:10.1155/2024/7096854
> Saj S, Jagoret P, Todem Ngnogue H, et al. (2013). *Carbon storage and density dynamics of associated trees in
> three contrasting Theobroma cacao agroforests of Central Cameroon.* **Agroforestry Systems 87(6):1309–1320.**
> doi:10.1007/s10457-013-9639-4

---

## 6. Colombia / neotropical-specific tree equations (Álvarez et al. 2012)

Álvarez et al. (2012) built tree allometries from **631 destructively harvested trees across Colombia**
(DBH ≥ 10 cm), following the Chave (2005) approach but classifying by **Holdridge life zone**. These are the
closest thing to a nationally calibrated tree equation for Colombia.

**Type I models (use D, H and ρ — directly comparable to Chave 2014):**

```
All types:   ln(AGB) = −2.277 + 2.017·ln(D) + 0.715·ln(H) + 0.718·ln(ρ)
All types:   ln(AGB) = −2.231 + 0.933·ln(D²·H·ρ)                       (compound form)
```

**Type II models (D and ρ only — no height):**

```
All types:   ln(AGB) = −3.103 + 1.794·ln(D) + 1.290·(ln D)² − 0.128·(ln D)³ + 0.819·ln(ρ)
```

**Life-zone-specific Type II.1 coefficients** (useful for a DBH-only fallback, no E needed).
Form: `ln(AGB) = a + b1·ln(D) + 1.169·(ln D)² − 0.122·(ln D)³ + d·ln(ρ)`:

| Life zone | a | b1 | d (ρ) |
|---|---|---|---|
| Tropical moist (Tm) | −2.406 | 1.289 | 0.445 |
| Premontane moist (PMm) | −1.960 | 1.098 | 1.061 |

- Units: D cm, H m, ρ g cm⁻³, AGB kg. Life zones from Holdridge et al. (1971).
- **Life-zone matching caveat:** a cacao agroforest at ~900 m in Santander sits in the **premontane**
  belt (PMm / premontane wet), which is *between* the Tm and PMm calibration sets. Álvarez's dataset also
  excludes premontane-wet (Pw) specifically. Pick the zone whose climate best matches, and document the
  choice.
- Álvarez et al. also found the **Chave (2005) no-height (Type II) models overestimate Colombian AGB by
  54.7 ± 135.7 %** — a strong reason to prefer the height-based equations or the Colombia-specific models.

> Álvarez E, Duque A, Saldarriaga J, et al. (2012). *Tree above-ground biomass allometries for carbon stocks
> estimation in the natural forests of Colombia.* **Forest Ecology and Management 267: 297–308.**
> doi:10.1016/j.foreco.2011.12.013 · PDF: https://ctfs.si.edu/Public/pdfs/ToDelete/Alvarez%20et%20al%202012.pdf

A recent **Colombian Amazon cacao-agroforestry** study (Caquetá, 78 plots) applied exactly this combination —
the Álvarez Tm tree equation for trees, the Andrade cacao equation for shrubs, and CF = 0.47 with the 3.67
CO₂:C ratio — and is a good methodological template:

> Díaz-Cháux JT, Velasquez-Valencia A, Casanoves F (2025). *Influence of landscape structure on carbon storage
> in agroforestry systems with cacao and silvopastoral systems in the Colombian Amazon.* **PLoS ONE 20(6):
> e0325477.** doi:10.1371/journal.pone.0325477

---

## 7. Units cheat-sheet (read carefully)

| Quantity | Symbol | Unit to use in the equations | Common pitfall |
|---|---|---|---|
| Diameter | D / DBH | **cm** (not m, not mm) | Some code divides by 10 (mm→cm); Chave Eq. 4 is defined in **cm** |
| Height | H | **m** | total height, not commercial/bole height |
| Wood specific gravity | ρ | **g cm⁻³** (= t m⁻³ = Mg m⁻³) | ≠ kg m⁻³; if you have kg m⁻³ divide by 1000 |
| AGB output | AGB | **kg** per tree | divide by 1000 for Mg (t) |
| Carbon fraction | CF | dimensionless **0.47** | not 0.5 in current IPCC default |
| CO₂ conversion | — | **44/12 = 3.6667** | not 44/44 or 12/44 |
| ρ·D²·H | — | cm²·m·g cm⁻³ | Chave 2005 valid range 50–1,000,000 |

ρ is defined as **oven-dry mass / green volume** (basic specific gravity), *not* air-dry density
(Chave et al. 2005).

---

## 8. Uncertainty and caveats

1. **Equation error (irreducible).** Chave 2014 Eq. 4: RSE = 0.357 (ln units), mean tree-level CV ≈ **56.5 %**,
   mean bias ≈ **5.31 %**. The no-height Eq. 7 is much worse (CV ≈ 71.5 %, bias ≈ 9.7 %). Per-tree AGB from
   any pantropical equation therefore carries a **~±40–60 %** uncertainty at the 1-σ level, even with perfect
   inputs. Aggregating many trees reduces the *stand-mean* error (errors partly cancel) but not the individual
   error.

2. **Wood-density uncertainty.** With species unknown, ρ uncertainty is **±10–40 %** (IPCC/FAO). Because
   AGB ∝ ρ^0.976, this propagates almost 1:1 into AGB. This is often the single largest controllable error.

3. **Destructive-sampling origin.** All these equations are calibrated on **felled, oven-dried** trees. They
   encode the average architecture of the calibration trees. The largest/rarest trees are under-represented;
   extrapolation to very large DBH is uncertain.

4. **Agroforestry vs. natural forest structure.** Pantropical equations were fitted mostly on closed-canopy
   natural forest trees. Open-grown agroforestry shade trees are often **shorter and broader-crowned for a
   given DBH** (more light, less competition), so height-based equations can be biased; conversely, the
   diameter–height relationship differs from forest. Documented cases: Chave-type equations **overestimated AGB
   by ~40 % in Central African forests** (Fayolle et al. 2013, doi:10.1016/j.foreco.2013.05.036). The effect
   is system-specific — there is no single correction factor, hence the recommendation to validate locally.

5. **Photogrammetric inputs propagate measurement error.** If DBH and/or H are measured from point clouds /
   photogrammetry (rather than calipers and clinometer), add their errors: AGB ∝ D^~2 (in Eq. 4, D enters as
   (D²)^0.976 = D^1.952), so a **+5 % DBH error → ~+10 % AGB error**; height errors propagate roughly linearly
   (H^0.976). Tree-top/occlusion errors in photogrammetric height are common and often systematic (under-
   or over-estimation), which biases the whole inventory rather than just adding noise.

6. **Range of validity.** Chave 2014: D ≥ 5 cm. Chave 2005: D 5–156 cm and ρD²H 50–1,000,000. Do not apply
   outside these ranges.

7. **Log-normal back-transformation bias.** For any log–log equation reported as `ln(AGB) = a + …`, the
   arithmetic-mean prediction is biased low unless the Baskerville (1972) correction `exp(σ²/2)` is applied.
   Chave 2014 Eq. 4 is already a power law (no correction needed); Eq. 7/§2.2 and the Chave 2005 models are in
   log form — apply the correction or use the corrected coefficients in §2.2.

8. **Lianas, palms, conifers and tree ferns are excluded.** The equations are for broadleaf trees. If the
   inventory contains palms or other non-dicot growth forms, they need separate equations.

---

## 9. Recommended default (copy-paste for Methods)

**AGB per tree (primary — use because you have both DBH and H):**

```
AGB_kg = 0.0673 × (ρ × D² × H)^0.976
```
- D = DBH in **cm**; H = total height in **m**; ρ = **0.60 g cm⁻³** (generic tropical default; regional
  alternative 0.63 g cm⁻³, Zanne et al. 2009); AGB in **kg**.
- Source: Chave et al. (2014), *Global Change Biology* 20:3177–3190, Eq. 4, doi:10.1111/gcb.12629.

**DBH-only fallback (when H is missing for a tree):**

```
AGB_kg = exp( −2.024 − 0.896·E + 0.920·ln(ρ) + 2.795·ln(D) − 0.0461·(ln D)² )
E = (0.178·TS − 0.938·CWD − 6.61·PS) × 10⁻³
```
- Source: Chave et al. (2014) Eq. 7, corrected per Réjou-Méchain et al. (2017), *Methods Ecol. Evol.*
  8:1163–1167, doi:10.1111/2041-210X.12753. Obtain `E` for 6.86 °N, −73.41 °W from the Chave `E` raster
  (http://chave.ups-tlse.fr/pantropical_allometry.htm) or compute from WorldClim TS/CWD/PS.
- If `E` is also unavailable, prefer a locally fitted H–D relationship feeding §9 primary, or the
  Colombia-specific Álvarez et al. (2012) Type-II equation (§6) — **not** the Chave (2005) no-height model
  (large overestimation in Colombia).

**Carbon and CO₂ equivalent:**

```
C_kg      = AGB_kg × 0.47
CO2e_kg   = C_kg × (44/12) = AGB_kg × 1.7233
```
- CF = 0.47: IPCC (2006) Guidelines, Vol. 4, Ch. 4, Table 4.3.
- 44/12: IPCC/UNFCCC molar-mass ratio.

**One-line summary for the Methods section:**
> "Above-ground biomass of each tree was estimated with the pantropical equation of Chave et al. (2014,
> Eq. 4), AGB = 0.0673·(ρ·D²·H)^0.976 (D in cm, H in m, ρ in g cm⁻³, AGB in kg), using a wood specific
> gravity of 0.60 g cm⁻³ (species unknown). Where height was unavailable, the no-height model of Chave et al.
> (2014, Eq. 7, as corrected by Réjou-Méchain et al. 2017) was used. Biomass was converted to carbon with the
> IPCC default carbon fraction of 0.47, and to CO₂-equivalent using the 44/12 molar ratio."

---

## References

1. Chave J, Réjou-Méchain M, Búrquez A, et al. (2014). Improved allometric models to estimate the aboveground
   biomass of tropical trees. *Global Change Biology* 20(10): 3177–3190. doi:10.1111/gcb.12629.
   https://onlinelibrary.wiley.com/doi/10.1111/gcb.12629
2. Chave J, Andalo C, Brown S, et al. (2005). Tree allometry and improved estimation of carbon stocks and
   balance in tropical forests. *Oecologia* 145(1): 87–99. doi:10.1007/s00442-005-0100-x.
3. Réjou-Méchain M, Tanguy A, Piponiot C, Chave J, Hérault B (2017). biomass: an R package for estimating
   above-ground biomass and its uncertainty in tropical forests. *Methods in Ecology and Evolution* 8(9):
   1163–1167. doi:10.1111/2041-210X.12753.
4. Zanne AE, Lopez-Gonzalez G, Coomes DA, et al. (2009). Towards a worldwide wood economics spectrum.
   *Ecology Letters* 12(4): 351–366. doi:10.1111/j.1461-0248.2009.01285.x. Data: Global Wood Density Database,
   Dryad doi:10.5061/dryad.234 (v2: Fischer et al. 2026, doi:10.5281/zenodo.16919509).
5. Chave J, Muller-Landau HC, Baker TR, Easdale TA, ter Steege H, Webb CO (2006). Regional and phylogenetic
   variation of wood density across 2456 Neotropical tree species. *Ecological Applications* 16(6):
   2356–2367. doi:10.1890/1051-0761(2006)016[2356:RAPVOW]2.0.CO;2.
6. Álvarez E, Duque A, Saldarriaga J, et al. (2012). Tree above-ground biomass allometries for carbon stocks
   estimation in the natural forests of Colombia. *Forest Ecology and Management* 267: 297–308.
   doi:10.1016/j.foreco.2011.12.013.
7. IPCC (2006). *2006 IPCC Guidelines for National Greenhouse Gas Inventories, Vol. 4 (AFOLU), Ch. 4 (Forest
   Land)*, Table 4.3 (carbon fraction 0.47), Table 4.13 (wood densities).
   https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/4_Volume4/V4_04_Ch4_Forest_Land.pdf
8. Fomekong Tane AB, et al. (2024). Site-Species Allometry Equation for Theobroma cacao L. Biomass Estimation
   in Agroforestry Systems of Cameroon. *International Journal of Ecology* 2024: 7096854.
   doi:10.1155/2024/7096854.
9. Saj S, Jagoret P, Todem Ngnogue H, et al. (2013). Carbon storage and density dynamics of associated trees
   in three contrasting Theobroma cacao agroforests of Central Cameroon. *Agroforestry Systems* 87(6):
   1309–1320. doi:10.1007/s10457-013-9639-4.
10. Díaz-Cháux JT, Velasquez-Valencia A, Casanoves F (2025). Influence of landscape structure on carbon
    storage in agroforestry systems with cacao and silvopastoral systems in the Colombian Amazon. *PLoS ONE*
    20(6): e0325477. doi:10.1371/journal.pone.0325477.
11. Niether W, et al. (2016). Cacao cultivation under diverse shade tree cover allows high carbon storage and
    sequestration without yield losses. *PLoS ONE* 11(2): e0149949. doi:10.1371/journal.pone.0149949.
12. Saj S, Jagoret P, Todem Ngnogue H, et al. (2022). Above-ground biomass models for dominant tree species in
    cacao agroforestry systems in Talamanca, Costa Rica. *Agroforestry Systems* 96: 655–670.
    doi:10.1007/s10457-022-00741-y.
13. Fayolle A, Doucet J-L, Gillet J-F, et al. (2013). Tree allometry in Central Africa: testing the validity
    of pantropical multi-species allometric equations for estimating biomass and carbon stocks. *Forest
    Ecology and Management* 305: 29–37. doi:10.1016/j.foreco.2013.05.036.
14. Baskerville GL (1972). Use of logarithmic regression in the estimation of plant biomass. *Canadian
    Journal of Forest Research* 2(1): 49–53. doi:10.1139/x72-009.
15. Winrock International (2016). *Allometric Equation Evaluation Guidance Document* (Walker SM et al.).
    https://winrock.org/wp-content/uploads/2018/08/Winrock-AllometricEquationGuidance-2016.pdf
