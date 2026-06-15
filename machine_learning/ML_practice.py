import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_score


def main():

    # ==================================================
    # 1. Boston Housing - Linear Regression with 10-Fold CV
    # ==================================================
    print("=" * 60)
    print("BOSTON HOUSING DATASET - 10 FOLD CROSS VALIDATION")
    print("=" * 60)

    df = pd.read_csv("BostonHousing.csv")

    X = df.drop("medv", axis=1)
    y = df["medv"]

    model = LinearRegression()

    kf = KFold(
        n_splits=10,
        shuffle=True,
        random_state=42
    )

    r2_scores = cross_val_score(
        model,
        X,
        y,
        cv=kf,
        scoring="r2"
    )

    print("R² score for each fold:")
    for i, score in enumerate(r2_scores, start=1):
        print(f"Fold {i}: {score:.4f}")

    print("\nResults")
    print("-" * 30)
    print(f"Mean R² Score      : {np.mean(r2_scores):.4f}")
    print(f"Standard Deviation : {np.std(r2_scores):.4f}")

    # ==================================================
    # 2. Matrix Multiplication (Xθ)
    # ==================================================
    print("\n" + "=" * 60)
    print("MATRIX MULTIPLICATION (Xθ)")
    print("=" * 60)

    theta = [
        [2],
        [3],
        [3]
    ]

    X_matrix = [
        [1, 0, 2],
        [0, 1, 1],
        [2, 1, 0],
        [1, 1, 1],
        [0, 2, 1]
    ]

    X_theta = []

    for i in range(len(X_matrix)):
        total = 0

        for j in range(len(theta)):
            total += X_matrix[i][j] * theta[j][0]

        X_theta.append([total])

    print("Theta:")
    for row in theta:
        print(row)

    print("\nX:")
    for row in X_matrix:
        print(row)

    print("\nXθ:")
    for row in X_theta:
        print(row)

    # ==================================================
    # 3. Plot y = 2x1² + 3x1 + 4
    # ==================================================
    print("\n" + "=" * 60)
    print("PLOTTING y = 2x1² + 3x1 + 4")
    print("=" * 60)

    x1 = np.linspace(-10, 10, 100)

    y_plot = 2 * (x1 ** 2) + 3 * x1 + 4

    plt.figure(figsize=(8, 5))
    plt.plot(x1, y_plot)
    plt.title("y = 2x1² + 3x1 + 4")
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()