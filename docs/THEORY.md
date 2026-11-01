# Theoretical Foundation of Non-Linear Effects (NLE)

## 1. Introduction and Phenomenological Definition

In asymmetric catalysis, the enantiomeric excess of the reaction product ($ee_{\text{prod}}$)
is frequently assumed to follow a linear relationship with the enantiomeric excess of the
chiral catalyst ($ee_{\text{cat}}$):

$$
ee_{\text{prod}} = ee_0 \cdot ee_{\text{cat}}
$$

where $ee_0$ represents the intrinsic enantioselectivity of the pure chiral catalyst ($ee_{\text{cat}} = 1.0$).
When $ee_0 = 1.0$, the ideal linear relationship simplifies to $ee_{\text{prod}} = ee_{\text{cat}}$.

However, seminal investigations by **Henri B. Kagan** (1986), followed by extensive work by
**Ryoji Noyori** (1989) and **Donna G. Blackmond** (2000), revealed that multi-component equilibria,
catalyst self-assembly, and oligomerization induce pronounced deviations from linearity:

- **Positive Non-Linear Effect ((+)NLE, Asymmetric Amplification):**
  $ee_{\text{prod}} > ee_0 \cdot ee_{\text{cat}}$. A catalyst of modest optical purity
  yields product with extraordinarily high enantiomeric purity.
- **Negative Non-Linear Effect ((-)NLE, Asymmetric Depletion):**
  $ee_{\text{prod}} < ee_0 \cdot ee_{\text{cat}}$. The product enantiopurity is severely degraded
  compared to the ideal linear expectation.

---

## 2. Kagan ML2 Reservoir Model

The most ubiquitous and foundational mechanistic framework for non-linear effects is the
**Kagan $\text{ML}_2$ Reservoir Model**, where a chiral monomer catalyst $M$ ($M_R$ or $M_S$)
equilibrates reversibly with dimeric aggregates.

### 2.1 Chemical Equilibria

The system comprises three dimeric species:

$$
2 M_R \xrightleftharpoons{K_{\text{homo}}} (R,R)
$$

$$
2 M_S \xrightleftharpoons{K_{\text{homo}}} (S,S)
$$

$$
M_R + M_S \xrightleftharpoons{K_{\text{hetero}}} (R,S)
$$

The corresponding association constants are defined as:

$$
K_{\text{homo}} = \frac{[(R,R)]}{[M_R]^2} = \frac{[(S,S)]}{[M_S]^2}
$$

$$
K_{\text{hetero}} = \frac{[(R,S)]}{[M_R][M_S]}
$$

### 2.2 Thermodynamic Relation & Statistical Symmetry

The free energy difference between the heterochiral and homochiral dimers is:

$$
\Delta G = G_{(R,S)} - G_{(R,R)}
$$

Because homochiral dimers $(R,R)$ and $(S,S)$ have a rotational $C_2$ symmetry axis
($\sigma_{\text{homo}} = 2$), whereas the heterochiral dimer $(R,S)$ possesses an inversion
center or plane of symmetry without a $C_2$ axis ($\sigma_{\text{hetero}} = 1$), a fundamental
statistical symmetry factor of 2 emerges:

$$
\frac{K_{\text{hetero}}}{K_{\text{homo}}} = 2 \exp\left(-\frac{\Delta G}{R T}\right)
$$

- If $\Delta G < 0$: $K_{\text{hetero}} > 2 K_{\text{homo}}$. The heterochiral dimer is
  thermodynamically favored over the homochiral dimer.
- If $\Delta G > 0$: $K_{\text{hetero}} < 2 K_{\text{homo}}$. The homochiral dimer is favored.
- If $\Delta G = -R T \ln(2) \approx -0.41\text{ kcal/mol}$ (at 298.15 K): $K_{\text{hetero}} = 2 K_{\text{homo}}$,
  yielding the purely statistical, non-interacting linear limit.

---

## 3. Mass Balance and Numerical Solver

Let $C_0$ be the total formal concentration of chiral catalyst ($C_0 = C_R + C_S = 1.0$).
For any specified catalyst enantiomeric excess $ee_{\text{cat}} \in [-1, 1]$:

$$
C_R = \frac{1 + ee_{\text{cat}}}{2}, \quad C_S = \frac{1 - ee_{\text{cat}}}{2}
$$

Conservation of mass requires:

$$
f_1(M_R, M_S) = M_R + 2 K_{\text{homo}} M_R^2 + K_{\text{hetero}} M_R M_S - C_R = 0
$$

$$
f_2(M_R, M_S) = M_S + 2 K_{\text{homo}} M_S^2 + K_{\text{hetero}} M_R M_S - C_S = 0
$$

### 3.1 Analytical Jacobian Matrix

To achieve quadratic convergence and avoid numerical instability across extreme association
constants ($K_{\text{homo}} \ge 10^7$), `nle-predictor` evaluates the exact analytical Jacobian:

$$
J = \begin{pmatrix}
\frac{\partial f_1}{\partial M_R} & \frac{\partial f_1}{\partial M_S} \\
\frac{\partial f_2}{\partial M_R} & \frac{\partial f_2}{\partial M_S}
\end{pmatrix}
= \begin{pmatrix}
1 + 4 K_{\text{homo}} M_R + K_{\text{hetero}} M_S & K_{\text{hetero}} M_R \\
K_{\text{hetero}} M_S & 1 + 4 K_{\text{homo}} M_S + K_{\text{hetero}} M_R
\end{pmatrix}
$$

### 3.2 Product Enantiopurity

In the classical reservoir regime, dimers are catalytically inert (or serve as inactive resting states).
The active monomer pool determines the product enantiopurity:

$$
ee_{\text{prod}} = ee_0 \cdot \frac{M_R - M_S}{M_R + M_S}
$$

---

## 4. Kagan Active Dimer Model

When dimeric complexes are themselves catalytically competent with relative turnover rates
$g = k_{\text{dimer}} / k_{\text{monomer}}$ and $g' = k_{\text{meso}} / k_{\text{monomer}}$,
the instantaneous enantiomeric excess is given by:

$$
ee_{\text{prod}} = \frac{(M_R - M_S) \left[ ee_0 + 2 g \cdot ee_D \cdot K_{\text{homo}} (M_R + M_S) \right]}{(M_R + M_S) + 2 g K_{\text{homo}} (M_R^2 + M_S^2) + 2 g' K_{\text{hetero}} M_R M_S}
$$

where $ee_D$ is the intrinsic enantioselectivity of the homochiral dimer.

---

## 5. Mathematical Curve Diagnostics

`nle-predictor` provides quantitative differential and integral metrics:

1. **Initial Slope at Racemic State ($S_0$):**
   $$S_0 = \left. \frac{d(ee_{\text{prod}})}{d(ee_{\text{cat}})} \right|_{ee_{\text{cat}} = 0}$$
   - $S_0 > 1.0$: Initial asymmetric amplification.
   - $S_0 < 1.0$: Initial asymmetric depletion.

2. **Maximum Amplification Factor ($\alpha_{\text{max}}$):**
   $$\alpha_{\text{max}} = \max_{ee_{\text{cat}} > 0} \left( \frac{ee_{\text{prod}}}{ee_{\text{cat}}} \right)$$

3. **Area Under Asymmetry (AUA):**
   $$\text{AUA} = \int_0^1 (ee_{\text{prod}}(x) - x) \, dx$$
   Provides a scalar index of net amplifying power ($\text{AUA} > 0$) or depletion ($\text{AUA} < 0$).

---

## 6. Standard-State Concentration Corrections

Quantum chemistry calculations typically compute ideal gas Gibbs free energies at standard
pressure ($P^\circ = 1\text{ atm}$, $T = 298.15\text{ K}$). Solution-phase kinetics require standard
molar concentration ($C^\circ = 1\text{ M} = 1\text{ mol/L}$).

For dimerization ($2 M \to D$, $\Delta n = -1$):

$$
\Delta G^\circ_{\text{soln}} = \Delta G^\circ_{\text{gas}} + \Delta G_{\text{corr}}
$$

$$
\Delta G_{\text{corr}} = -\Delta n \cdot R T \ln\left( \frac{R_{\text{atm}} T}{1.0} \right) = +R T \ln(24.4654) \approx +1.894\text{ kcal/mol}
$$

---

## References

1. Kagan, H. B.; Dang, T. P. *J. Am. Chem. Soc.* **1972**, *94*, 6429–6433.
2. Puchot, C.; Samuel, O.; Dunach, E.; Zhao, S.; Agami, C.; Kagan, H. B. *J. Am. Chem. Soc.* **1986**, *108*, 2353–2357.
3. Guillaneux, D.; Zhao, S.-H.; Samuel, O.; Rainford, D.; Kagan, H. B. *J. Am. Chem. Soc.* **1994**, *116*, 9430–9439.
4. Kitamura, M.; Okada, S.; Suga, S.; Noyori, R. *J. Am. Chem. Soc.* **1989**, *111*, 4028–4036.
5. Blackmond, D. G. *Acc. Chem. Res.* **2000**, *33*, 402–411.
6. Walsh, P. J.; Lurain, J. R. *Chem. Rev.* **2003**, *103*, 3297–3344.
7. Denmark, S. E.; Beutner, G. L. *Angew. Chem. Int. Ed.* **2008**, *47*, 1560–1638.
