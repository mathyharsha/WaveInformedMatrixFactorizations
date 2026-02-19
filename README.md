# Wave-Informed Matrix Factorization (WIMF)

A MATLAB implementation of Wave-Informed Matrix Factorization, a physics-informed dictionary learning algorithm for 3D tensor data. The method embeds wave-equation regularization into the matrix factorization objective, making it well-suited for applications such as structural health monitoring, seismic analysis, and spatiotemporal signal decomposition.

---

## Overview

WIMF decomposes a spatiotemporal tensor $\mathcal{Y} \in \mathbb{R}^{N_y \times N_x \times N_t}$ into a dictionary of wave-equation-constrained atoms and their corresponding coefficients. The key idea is to work in a spectral domain defined by the eigenvectors of discrete Laplacian operators, where the wave equation takes a particularly simple diagonal form. Dictionary atoms are encouraged to satisfy the wave equation by penalizing deviations from the dispersion relation $(\Delta_x + \Delta_y) u = \frac{1}{c^2} \partial_{tt} u$ for some wave speed $c$.

**The objective function** (minimized over dictionary $D$, coefficients $x$, and wave speeds $C$) is:

$$\min_{D, x, C} \; \frac{1}{2} \| \hat{y} - \mathcal{M}(Dx) \|^2 + \frac{\lambda}{2} \left( \sum_j \gamma \| (A_1 + A_2 - c_j^{-2} A_3) d_j \|^2 + \|D\|_F^2 + \|x\|^2 \right)$$

where $\mathcal{M}(\cdot)$ is a spectral masking operator, $A_1, A_2, A_3$ are Kronecker-structured matrices encoding the spatial and temporal Laplacian eigenvalues, and $\hat{y}$ is the spectrally transformed data.

---

## File Structure

| File | Description |
|---|---|
| `waveInformedMatFac.m` | Main function implementing the greedy atom extraction loop and joint gradient-descent refinement |
| `spectral_wave_transform.m` | Forward spectral transform: maps a spatial tensor to the joint spectral domain via the Laplacian eigenbases $(V_y, V_x, V_t)$ |
| `inv_spectral_wave_transform.m` | Inverse spectral transform: maps a spectral-domain vector back to the spatial tensor domain |
| `spectralMask.m` | Applies an elementwise spatial mask $T$ to each dictionary column in the spectral domain |
| `forPolar_actual.m` | Evaluates the polarization functional $z^\top B^{-1} z$ for a given wave speed $c$, used during greedy atom extraction |

---

## Algorithm

The algorithm proceeds in two stages per iteration:

**1. Greedy Atom Extraction**

For each iteration, the residual in the spectral domain is computed and a new wave speed $c$ is found by maximizing the polarization functional $z^\top B(c)^{-1} z$ via simulated annealing (`simulannealbnd`). The optimal $c$ determines a new dictionary atom $d$ via the wave-regularized filter, and a coefficient $\tau$ is appended.

**2. Joint Refinement (optional)**

If `gradient_descent` is enabled, all dictionary columns, wave speeds, and coefficients are jointly refined using gradient descent with backtracking line search to minimize the full objective.

---

## Usage

```matlab
[Dc, x, C] = waveInformedMatFac(data, Trans, V_x, V_y, V_t, Name, Value, ...)
```

### Inputs

| Argument | Type | Description |
|---|---|---|
| `data` | $N_y \times N_x \times N_t$ array | Input spatiotemporal tensor |
| `Trans` | $N_y \times N_x \times N_t$ array | Binary or soft spatial mask applied in the spatial domain |
| `V_x`, `V_y`, `V_t` | matrices | Eigenvector matrices of the 1D discrete Laplacians along $x$, $y$, $t$ (computed internally if not pre-supplied) |

### Optional Name-Value Arguments

| Name | Default | Description |
|---|---|---|
| `'count'` | `5` | Maximum number of dictionary atoms to extract |
| `'threshold'` | `0.1` | Stopping tolerance on the polarization condition ($\text{polar} \leq 1 + \text{threshold}$) |
| `'tolerance'` | `0.1` | Convergence tolerance for the gradient descent refinement step |
| `'gradient_descent'` | `false` | Whether to run joint gradient-descent refinement after each greedy step |

### Outputs

| Output | Description |
|---|---|
| `Dc` | Dictionary matrix; each column is a spectral-domain atom |
| `x` | Coefficient vector |
| `C` | Diagonal matrix of wave speeds $c_j$ for each atom |

### Example

```matlab
% Load or generate a spatiotemporal tensor
data = rand(32, 32, 50);   % e.g. 32x32 spatial grid, 50 time steps
Trans = ones(size(data));  % no masking

% Run WIMF with greedy extraction only
[Dc, x, C] = waveInformedMatFac(data, Trans, [], [], [], ...
    'count', 8, 'threshold', 0.05);

% Run with joint gradient-descent refinement
[Dc, x, C] = waveInformedMatFac(data, Trans, [], [], [], ...
    'count', 8, 'threshold', 0.05, 'gradient_descent', true, 'tolerance', 0.01);
```

---

## Dependencies

- MATLAB R2019b or later
- **Global Optimization Toolbox** — required for `simulannealbnd` (wave-speed search)
- **Signal Processing Toolbox** — uses `del2` for constructing discrete Laplacians

---

## Notes

- The Laplacian eigenbases $(V_x, V_y, V_t)$ are computed internally from `del2`-based finite-difference operators scaled by $4$. Any pre-computed eigenbases passed as arguments are overwritten inside the function.
- The wave-speed search range is derived from the ratio of spatial-to-temporal Laplacian eigenvalues: $c \in [1/\sqrt{b},\; 1/\sqrt{a}]$ where $a, b$ are the min/max of $(\lambda_x + \lambda_y)/\lambda_t$.
- Regularization parameters $\gamma$ and $\lambda$ are set automatically based on the data size and spectral norm, but can be modified directly in `waveInformedMatFac.m`.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{wimf,
  author = {Harsha Vardhan Tetali, Kang Yang, Joel B. Harley},
  title  = {Wave Physics-Informed Tensor Completion of Sparse Wavefield Videos},
  year   = {},
  url    = {}
}
```

---

## License

[MIT License](LICENSE)
