import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score


# =====================================================
# FROM SCRATCH FUNCTIONS
# =====================================================

def hypothesis(X, theta):
    return np.dot(X, theta)


def compute_cost(X, y, theta):
    m = len(y)
    predictions = hypothesis(X, theta)
    cost = np.sum((predictions - y) ** 2) / (2 * m)
    return cost


def compute_derivative(X, y, theta):
    m = len(y)
    predictions = hypothesis(X, theta)
    errors = predictions - y
    gradient = np.dot(X.T, errors) / m
    return gradient


# =====================================================
# MAIN
# =====================================================

def main():

    # Load dataset
    df = pd.read_csv("simulated_data_multiple_linear_regression_for_ML.csv")

    print("Dataset Preview:")
    print(df.head())

    # =================================================
    # PART 1 : SCIKIT-LEARN FOR disease_score
    # =================================================

    print("\n==============================")
    print("Predicting disease_score")
    print("==============================")

    X1 = df.drop("disease_score", axis=1)
    y1 = df["disease_score"]

    X1_train, X1_test, y1_train, y1_test = train_test_split(
        X1, y1, test_size=0.2, random_state=42
    )

    model1 = LinearRegression()
    model1.fit(X1_train, y1_train)

    y1_pred = model1.predict(X1_test)

    print("MSE:", mean_squared_error(y1_test, y1_pred))
    print("R² Score:", r2_score(y1_test, y1_pred))

    print("\nIntercept:")
    print(model1.intercept_)

    print("\nCoefficients:")
    for feature, coef in zip(X1.columns, model1.coef_):
        print(f"{feature}: {coef}")

    print("\nFirst 5 Predictions:")
    print(model1.predict(X1.head()))

    # =================================================
    # PART 2 : SCIKIT-LEARN FOR disease_score_fluct
    # =================================================

    print("\n==============================")
    print("Predicting disease_score_fluct")
    print("==============================")

    X2 = df.drop("disease_score_fluct", axis=1)
    y2 = df["disease_score_fluct"]

    X2_train, X2_test, y2_train, y2_test = train_test_split(
        X2, y2, test_size=0.2, random_state=42
    )

    model2 = LinearRegression()
    model2.fit(X2_train, y2_train)

    y2_pred = model2.predict(X2_test)

    print("MSE:", mean_squared_error(y2_test, y2_pred))
    print("R² Score:", r2_score(y2_test, y2_pred))

    print("\nIntercept:")
    print(model2.intercept_)

    print("\nCoefficients:")
    for feature, coef in zip(X2.columns, model2.coef_):
        print(f"{feature}: {coef}")

    print("\nFirst 5 Predictions:")
    print(model2.predict(X2.head()))
    # =================================================
    # PART 3 : FROM SCRATCH GRADIENT DESCENT
    # =================================================

    print("\n==============================")
    print("Gradient Descent From Scratch")
    print("==============================")

    X = df.drop("disease_score_fluct", axis=1).values


    y = df["disease_score_fluct"].values

    # Feature Scaling
    X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

    # Add bias column
    X = np.c_[np.ones(X.shape[0]), X]
    v=X.shape[1]
    print(v)

    # Initialize theta
    theta = np.zeros(X.shape[1])

    alpha = 0.01
    iterations = 1000

    for i in range(iterations):

        gradient = compute_derivative(X, y, theta)

        # Update parameters
        theta = theta - alpha * gradient

        if i % 100 == 0:
            cost = compute_cost(X, y, theta)
            print(f"Iteration {i}: Cost = {cost:.4f}")

    print("\nFinal Theta:")
    print(theta)

    predictions = hypothesis(X, theta)

    print("\nFirst 5 Predictions:")
    print(predictions[:5])

    print("\nFirst 5 Actual Values:")
    print(y[:5])

    print("\nFinal Cost:")
    print(compute_cost(X, y, theta))

if __name__ == "__main__":
    main()




































# import pandas as pd
# import numpy as np
#
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import mean_squared_error, r2_score
#
#
# # =====================================================
# # FROM SCRATCH FUNCTIONS
# # =====================================================
#
# def hypothesis(X, theta):
#     return np.dot(X, theta)
#
#
# def compute_cost(X, y, theta):
#     predictions = hypothesis(X, theta)
#     cost = np.sum((predictions - y) ** 2) / 2
#     return cost
#
#
# def compute_derivative(X, y, theta):
#     predictions = hypothesis(X, theta)
#     errors = predictions - y
#     gradient = np.dot(X.T, errors)
#     return gradient
#
#
# # =====================================================
# # MAIN
# # =====================================================
#
# def main():
#
#     # Load dataset
#     df = pd.read_csv("simulated_data_multiple_linear_regression_for_ML.csv")
#
#     print("Dataset Preview:")
#     print(df.head())
#
#     # =================================================
#     # PART 1 : SCIKIT-LEARN FOR disease_score
#     # =================================================
#
#     print("\n==============================")
#     print("Predicting disease_score")
#     print("==============================")
#
#     X1 = df.drop("disease_score", axis=1)
#     y1 = df["disease_score"]
#
#     X1_train, X1_test, y1_train, y1_test = train_test_split(
#         X1, y1, test_size=0.2, random_state=42
#     )
#
#     model1 = LinearRegression()
#     model1.fit(X1_train, y1_train)
#
#     y1_pred = model1.predict(X1_test)
#
#     print("MSE:", mean_squared_error(y1_test, y1_pred))
#     print("R² Score:", r2_score(y1_test, y1_pred))
#
#     print("\nIntercept:")
#     print(model1.intercept_)
#
#     print("\nCoefficients:")
#     for feature, coef in zip(X1.columns, model1.coef_):
#         print(f"{feature}: {coef}")
#
#     print("\nFirst 5 Predictions:")
#     print(model1.predict(X1.head()))
#
#     # =================================================
#     # PART 2 : SCIKIT-LEARN FOR disease_score_fluct
#     # =================================================
#
#     print("\n==============================")
#     print("Predicting disease_score_fluct")
#     print("==============================")
#
#     X2 = df.drop("disease_score_fluct", axis=1)
#     y2 = df["disease_score_fluct"]
#
#     X2_train, X2_test, y2_train, y2_test = train_test_split(
#         X2, y2, test_size=0.2, random_state=42
#     )
#
#     model2 = LinearRegression()
#     model2.fit(X2_train, y2_train)
#
#     y2_pred = model2.predict(X2_test)
#
#     print("MSE:", mean_squared_error(y2_test, y2_pred))
#     print("R² Score:", r2_score(y2_test, y2_pred))
#
#     print("\nIntercept:")
#     print(model2.intercept_)
#
#     print("\nCoefficients:")
#     for feature, coef in zip(X2.columns, model2.coef_):
#         print(f"{feature}: {coef}")
#
#     print("\nFirst 5 Predictions:")
#     print(model2.predict(X2.head()))
#
#     # =================================================
#     # PART 3 : FROM SCRATCH GRADIENT DESCENT
#     # =================================================
#
#     print("\n==============================")
#     print("Gradient Descent From Scratch")
#     print("==============================")
#
#     X = df.drop(
#         ["disease_score", "disease_score_fluct"],
#         axis=1
#     ).values
#
#     y = df["disease_score_fluct"].values
#
#     # Add bias column
#     X = np.c_[np.ones(X.shape[0]), X]
#
#     # Initialize theta
#     theta = np.zeros(X.shape[1])
#
#     alpha = 0.0001
#     iterations = 1000
#
#     for i in range(iterations):
#
#         gradient = compute_derivative(X, y, theta)
#
#         # Update parameters
#         theta = theta - alpha * gradient
#
#         if i % 100 == 0:
#             cost = compute_cost(X, y, theta)
#             print(f"Iteration {i}: Cost = {cost}")
#
#     print("\nFinal Theta:")
#     print(theta)
#
#     predictions = hypothesis(X, theta)
#
#     print("\nFirst 5 Predictions:")
#     print(predictions[:5])
#
#     print("\nFirst 5 Actual Values:")
#     print(y[:5])
#
#
# if __name__ == "__main__":
#     main()