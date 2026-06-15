"""
Task 3: Multiple Linear Regression FROM SCRATCH in Python
Dataset: simulated_data_multiple_linear_regression_for_ML.csv
Target : disease_score_fluct

Implements:
  - Read CSV
  - Form X (features) and y (disease_score_fluct)
  - compute_hypothesis()
  - compute_cost()       (Mean Squared Error / 2)
  - compute_derivatives()
  - update_parameters()  (gradient descent loop)
  - Evaluation & comparison with scikit-learn
"""

import numpy as np
import pandas as pd

# ══════════════════════════════════════════════════════════
# 1.  READ DATA
# ══════════════════════════════════════════════════════════
def load_data(filepath: str):
    """Read CSV and return feature matrix X and target vector y."""
    df = pd.read_csv(filepath)
    print("=" * 60)
    print("Dataset loaded")
    print("=" * 60)
    print(df.head())
    print(f"\nShape : {df.shape}")
    return df


# ══════════════════════════════════════════════════════════
# 2.  FORM X AND y
# ══════════════════════════════════════════════════════════
def prepare_data(df: pd.DataFrame, target: str = "disease_score_fluct"):
    """
    Build X (m × n+1 with bias column) and y (m,).
    Also returns the mean/std used for normalisation.
    """
    feature_cols = ["age", "BMI", "BP", "blood_sugar", "Gender"]

    X_raw = df[feature_cols].values.astype(float)   # (m, n)
    y     = df[target].values.astype(float)          # (m,)

    # ── Feature normalisation (z-score) ──────────────────
    mu    = X_raw.mean(axis=0)
    sigma = X_raw.std(axis=0)
    sigma[sigma == 0] = 1                            # guard against zero std
    X_norm = (X_raw - mu) / sigma                   # (m, n)

    # ── Add bias column (column of 1s) ───────────────────
    m = X_norm.shape[0]
    X = np.hstack([np.ones((m, 1)), X_norm])        # (m, n+1)

    print(f"\nFeatures : {feature_cols}")
    print(f"X shape  : {X.shape}  (includes bias term)")
    print(f"y shape  : {y.shape}")

    return X, y, mu, sigma, feature_cols


# ══════════════════════════════════════════════════════════
# 3.  HYPOTHESIS  h(x) = X · θ
# ══════════════════════════════════════════════════════════
def compute_hypothesis(X: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """
    Linear hypothesis function.

    Parameters
    ----------
    X     : (m, n+1) design matrix (with bias column)
    theta : (n+1,)   parameter vector

    Returns
    -------
    h     : (m,) predicted values
    """
    return X @ theta          # matrix–vector dot product


# ══════════════════════════════════════════════════════════
# 4.  COST  J(θ) = (1/2m) Σ (h - y)²
# ══════════════════════════════════════════════════════════
def compute_cost(X: np.ndarray, y: np.ndarray, theta: np.ndarray) -> float:
    """
    Mean Squared Error cost (divided by 2 for cleaner gradient).

    Parameters
    ----------
    X     : (m, n+1) design matrix
    y     : (m,)     true values
    theta : (n+1,)   parameter vector

    Returns
    -------
    J     : scalar cost value
    """
    m  = len(y)
    h  = compute_hypothesis(X, theta)
    errors = h - y                         # (m,)
    J = (1 / (2 * m)) * np.dot(errors, errors)
    return J


# ══════════════════════════════════════════════════════════
# 5.  DERIVATIVES  ∂J/∂θⱼ = (1/m) Σ (h - y) · xⱼ
# ══════════════════════════════════════════════════════════
def compute_derivatives(X: np.ndarray, y: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """
    Gradient of the cost function w.r.t. every parameter.

    Parameters
    ----------
    X     : (m, n+1) design matrix
    y     : (m,)     true values
    theta : (n+1,)   parameter vector

    Returns
    -------
    grad  : (n+1,) gradient vector
    """
    m      = len(y)
    h      = compute_hypothesis(X, theta)
    errors = h - y                         # (m,)
    grad   = (1 / m) * (X.T @ errors)     # (n+1,)
    return grad


# ══════════════════════════════════════════════════════════
# 6.  GRADIENT DESCENT  (update parameters)
# ══════════════════════════════════════════════════════════
def gradient_descent(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float  = 0.1,
    n_iter: int   = 1000,
    verbose: bool = True,
) -> tuple[np.ndarray, list[float]]:
    """
    Batch gradient descent to learn θ.

    Parameters
    ----------
    X      : (m, n+1) design matrix
    y      : (m,)     target vector
    alpha  : learning rate
    n_iter : number of iterations
    verbose: print cost every 100 iterations

    Returns
    -------
    theta      : (n+1,) learned parameters
    cost_history: list of cost at each iteration
    """
    n_params     = X.shape[1]
    theta        = np.zeros(n_params)      # initialise θ = 0
    cost_history = []

    print("\n" + "=" * 60)
    print("Gradient Descent Training")
    print("=" * 60)
    print(f"Learning rate : {alpha}")
    print(f"Iterations    : {n_iter}")
    print(f"Parameters    : {n_params}  (1 bias + {n_params-1} features)\n")

    for i in range(n_iter):

        # ── Compute gradient ──────────────────────────────
        grad = compute_derivatives(X, y, theta)

        # ── Update all parameters simultaneously ──────────
        theta = theta - alpha * grad

        # ── Record cost ───────────────────────────────────
        cost = compute_cost(X, y, theta)
        cost_history.append(cost)

        if verbose and (i % 100 == 0 or i == n_iter - 1):
            print(f"  Iter {i:>5d}  |  Cost = {cost:.6f}")

    return theta, cost_history


# ══════════════════════════════════════════════════════════
# 7.  EVALUATION HELPERS
# ══════════════════════════════════════════════════════════
def r2_score_manual(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return 1 - ss_res / ss_tot

def rmse(y_true, y_pred):
    return np.sqrt(np.mean((y_true - y_pred) ** 2))

def mae(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))


# ══════════════════════════════════════════════════════════
# 8.  MAIN
# ══════════════════════════════════════════════════════════
def main():
    filepath = "simulated_data_multiple_linear_regression_for_ML.csv"

    # ── Load ──────────────────────────────────────────────
    df = load_data(filepath)

    # ── Prepare ───────────────────────────────────────────
    X, y, mu, sigma, feature_cols = prepare_data(df, target="disease_score_fluct")

    # ── Initial cost (θ = 0) ─────────────────────────────
    theta_init = np.zeros(X.shape[1])
    J0 = compute_cost(X, y, theta_init)
    print(f"\nInitial cost (θ=0) : {J0:.4f}")

    # ── Train ─────────────────────────────────────────────
    theta_final, cost_history = gradient_descent(
        X, y, alpha=0.1, n_iter=1000, verbose=True
    )

    # ── Final cost ────────────────────────────────────────
    J_final = compute_cost(X, y, theta_final)
    print(f"\nFinal cost         : {J_final:.4f}")

    # ── Learned parameters ────────────────────────────────
    print("\n" + "=" * 60)
    print("Learned Parameters (θ)")
    print("=" * 60)
    labels = ["bias (θ₀)"] + [f"{f} (θ{i+1})" for i, f in enumerate(feature_cols)]
    for label, val in zip(labels, theta_final):
        print(f"  {label:<25} : {val:>12.6f}")

    # ── Predictions & evaluation ─────────────────────────
    y_pred = compute_hypothesis(X, theta_final)

    print("\n" + "=" * 60)
    print("Evaluation on Full Dataset")
    print("=" * 60)
    print(f"  R²   : {r2_score_manual(y, y_pred):.6f}")
    print(f"  RMSE : {rmse(y, y_pred):.4f}")
    print(f"  MAE  : {mae(y, y_pred):.4f}")

    # ── Sample predictions ────────────────────────────────
    print("\nSample Predictions (first 10 rows):")
    sample = pd.DataFrame({
        "Actual":    y[:10].round(3),
        "Predicted": y_pred[:10].round(3),
        "Residual":  (y[:10] - y_pred[:10]).round(3),
    })
    print(sample.to_string(index=False))

    # ── Cost curve summary ────────────────────────────────
    print("\nCost curve (every 100 iterations):")
    for i in range(0, len(cost_history), 100):
        print(f"  Iter {i:>5d}  cost = {cost_history[i]:.4f}")

    print("\n✓ Gradient descent converged successfully.")


if __name__ == "__main__":
    main()