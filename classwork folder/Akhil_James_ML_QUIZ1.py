
import numpy as np
from sklearn.datasets import load_diabetes
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_score




def main():
    # ==================================================
    # 1. Diabetes Data Set - Linear Regression with 10-Fold CV
    # ==================================================
    print("=" * 60)
    print("Diabetes Data Set - Linear Regression with 10-Fold CV")
    print("=" * 60)

    daib= load_diabetes(as_frame=True)
    diabetes_X = daib.data
    y = daib.target
    X=diabetes_X.iloc[:,[2]]








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
    # 2. dot product
    # ==================================================
    print("\n" + "=" * 60)
    print("dot product of X and Y")
    print("=" * 60)

    x = [2, 1, 2]
    y = [1, 2, 2]

    result = 0
    print("\n")

    for i in range(len(x)):
        result += x[i] * y[i]
    print("Vector x:")
    for i in x:
        print(i)

    print("Vector y:")
    for i in y:
        print(i)
    print("Dot Product =", result)



    # ==================================================
    # 2b.Ploting
    # ==================================================

    print("plot = 3*(x1 ** 2)-2*x1+1.5")
    print("=" * 60)

    x1 = np.linspace(-10, 10, 100)

    y_plot = 3*(x1 ** 2)-2*x1+1.5

    plt.figure(figsize=(8, 5))
    plt.plot(x1, y_plot)
    plt.title("y = 3x1² - 2x1 + 1.5")
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()