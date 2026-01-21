
## 1. Overview
This case study investigates how sensitive the Time‑Varying Jackknife Model Averaging (TVJMA) estimator is to the choice of bandwidth in the kernel smoothing step.  
The reference TVJMA paper uses a global rule‑of‑thumb bandwidth:

$h  = c_h  T^{-0.2}$ where $c_h^{opt}=2.34$

but also notes that forecast performance can vary depending on the bandwidth.  
Our goal is to test whether this rule‑of‑thumb is robust across different structural change scenarios and sample sizes.

---

## 2. Research Question
**How does the choice of bandwidth affect TVJMA performance under smooth vs. abrupt structural change and different sample sizes?**

Sub‑questions:
- Does a **smooth** DGP prefer a **larger** bandwidth?
- Does an **abrupt** DGP require a **smaller** bandwidth?
- How does the performance change when \(T\) increases from 50 to 200?

---

## 3. Bandwidth Grid

With compact support and the scaling used in the paper, observations enter the local fit when
```math
\left|\frac{s - t}{T\,h}\right| \le 1
\;\;\Longleftrightarrow\;\;
|s - t| \le T h,
```
so the **effective local window size** is approximately **$$T h$$**:
- **Ultra‑local**: very few points, high variance
- **Moderately local**: balanced bias–variance
- **Global/time‑constant**: all observations enter every local fit


### T= 50 (from ~1 point to global)
**Ultra-local**:  $h = \\{1/50\\} = \\{0.02\\} \Rightarrow window \approx \\{1\\}$ i.e. for $t$, $s=\\{t+1, t, t-1\\}$

**Small local**: $h = \\{10/50=0.2, 0.4\\}  \Rightarrow window \approx \\{10, 20\\}$ 

**Moderate local**:  $h = \\{ 0.60, 0.80\\}  \Rightarrow window \approx \\{30, 40\\}$

**Global benchmarks**:  $h = \\{1.00\\} \Rightarrow window \approx \\{50\\}$ i.e. for $t$, $s=\\textit{all observations}$

### T= 200 (from ~1 point to global)
for the same h ranges! 

## 4. Simulation Design

### DGPs
We consider two distinct structural change patterns:

#### **DGP 1 — Smooth change** (DGP 1 in the original paper)
Parameters vary smoothly over time (e.g., macroeconomic drift).

#### **DGP 3 — Abrupt regime switches** (Variation of DGP 3 in the original paper)
A clear, abrupt break in regression parameters (e.g., financial regime shifts).
```


def DGP3(T, J, R2, alpha=1.5):
    c = np.sqrt(R2 / (1 - R2))

    t = np.arange(1, T + 1)
    tau = t / T

    F_tau = np.where(tau < 0.3, 0.1,
                     np.where(tau < 0.8, 1.0, -0.5))

    X = np.random.randn(T, J)
    X[:, 0] = 1.0

    j = np.arange(1, J+1)
    theta = c * np.sqrt(2 * alpha) * j**(-(alpha + 0.5))

    mu = F_tau * (X @ theta)
    eps = np.random.randn(T)

    Y = mu + eps
    return Y, X, mu

```
### Sample Sizes
- \(T = 50\)  
- \(T = 200\)

### Performance Metric
- Mean Squared Error (MSE) of estimated conditional mean  
- Optionally normalize against an oracle model, following the TVJMA paper

---

## 5. Expected Outcomes
- **Smooth DGP:**  
  Larger bandwidths should perform better (more smoothing reduces noise).

- **Abrupt DGP:**  
  Smaller bandwidths should outperform (avoid oversmoothing across breaks).

- **Effect of increasing (T):**
  - Smooth changes → performance improves as $Th$ grows  
  - Abrupt changes → risk of oversmoothing increases unless $h$ shrinks

---

## 6. Figures

### (1) One Simulation Example per DGP
For each DGP, include a single representative simulation run to illustrate the qualitative behavior of the underlying system
clarifing the structural features of each DGP—whether the evolution is smooth or abrupt

### (2) Performance Curves (fig 1. in the original paper)
* x‑axis: population R2, varied over a grid
* y‑axis: RMSE
* separate curves: one per bandwidth in the experimental grid
* separate panels: one per sample size (T)

### (3) Distributional Analysis via Boxplots of MSE (if time allows)

In addition to average performance, examine the distribution of MSE across Monte Carlo replications using boxplots. 
For each DGP and each sample size, produce boxplots comparing all bandwidths side by side.
This distributional perspective reveals features that average RMSE alone cannot capture, such as:
* variability across replications
* presence of outliers
* stability versus volatility of each bandwidth
* situations where a bandwidth has a low mean error but high dispersion (or vice versa)


