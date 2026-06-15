import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error


# =====================================================
# Gradient Descent Functions (From Scratch)
# =====================================================

def hypothesis(X, theta):
    return np.dot(X, theta)


def compute_cost(X, y, theta):

        predictions = hypothesis(X, theta)
        cost = np.sum((predictions - y) ** 2) / 2
        return cost


def compute_derivative(X, y, theta):
    m = len(y)
    predictions = hypothesis(X, theta)
    gradient = np.dot(X.T, (predictions - y)) / m
    return gradient


def main():

    # =====================================================
    # Read CSV
    # =====================================================
    df = pd.read_csv(
        "simulated_data_multiple_linear_regression_for_ML.csv"
    )

    print("Dataset Shape:", df.shape)

    # =====================================================
    # PART 1
    # Scikit-Learn Linear Regression
    # Target = disease_score
    # =====================================================

    print("\n========== SCIKIT-LEARN : disease_score ==========")

    X1 = df.drop(
        columns=["disease_score", "disease_score_fluct"]
    )

    y1 = df["disease_score"]

    X1_train, X1_test, y1_train, y1_test = train_test_split(
        X1,
        y1,
        test_size=0.2,
        random_state=42
    )

    model1 = LinearRegression()

    model1.fit(X1_train, y1_train)

    y1_pred = model1.predict(X1_test)

    print("R² Score:",
          r2_score(y1_test, y1_pred))

    print("MSE:",
          mean_squared_error(y1_test, y1_pred))

    # =====================================================
    # PART 2
    # Scikit-Learn Linear Regression
    # Target = disease_score_fluct
    # =====================================================

    print("\n========== SCIKIT-LEARN : disease_score_fluct ==========")

    X2 = df.drop(
        columns=["disease_score", "disease_score_fluct"]
    )

    y2 = df["disease_score_fluct"]

    X2_train, X2_test, y2_train, y2_test = train_test_split(
        X2,
        y2,
        test_size=0.2,
        random_state=42
    )

    model2 = LinearRegression()

    model2.fit(X2_train, y2_train)

    y2_pred = model2.predict(X2_test)

    print("R² Score:",
          r2_score(y2_test, y2_pred))

    print("MSE:",
          mean_squared_error(y2_test, y2_pred))

    # =====================================================
    # PART 3
    # Gradient Descent From Scratch
    # Target = disease_score_fluct
    # =====================================================

    print("\n========== GRADIENT DESCENT FROM SCRATCH ==========")

    X = df.drop(
        columns=["disease_score_fluct"]
    ).values

    y = df["disease_score_fluct"].values

    # Feature Scaling
    X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

    # Add Bias Column
    X = np.c_[np.ones(len(X)), X]

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.20,
        random_state=42
    )

    print("Training Samples:", len(X_train))
    print("Testing Samples :", len(X_test))

    theta = np.zeros(X_train.shape[1])

    learning_rate = 0.01
    iterations = 1000

    print(
        "\nInitial Cost:",
        compute_cost(X_train, y_train, theta)
    )

    for i in range(iterations):

        gradient = compute_derivative(
            X_train,
            y_train,
            theta
        )

        theta = theta - learning_rate * gradient

        if i % 100 == 0:

            cost = compute_cost(
                X_train,
                y_train,
                theta
            )

            print(
                f"Iteration {i:4d} Cost = {cost:.4f}"
            )

    print("\nFinal Theta Values:")
    print(theta)

    y_pred = hypothesis(X_test, theta)

    print("\nGradient Descent Performance")

    print(
        "R² Score:",
        r2_score(y_test, y_pred)
    )

    print(
        "MSE:",
        mean_squared_error(y_test, y_pred)
    )

    print("\nFirst 10 Predictions")
    print(y_pred[:10])

    print("\nFirst 10 Actual Values")
    print(y_test[:10])


if __name__ == "__main__":
    main()







































