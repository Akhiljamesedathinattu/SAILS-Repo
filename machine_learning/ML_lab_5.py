import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from sklearn.linear_model import (
    LogisticRegression,
    RidgeClassifier
)

from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    classification_report
)
import warnings

warnings.filterwarnings("ignore")


##################################################
# SIGMOID FUNCTION
##################################################

def sigmoid(x):
    return 1 / (1 + np.exp(-x))


##################################################
# DERIVATIVE OF SIGMOID
##################################################

def sigmoid_derivative(x):

    s = sigmoid(x)

    return s * (1 - s)


##################################################
# L1 NORM
##################################################

def l1_norm(vector):

    total = 0

    for value in vector:
        total += abs(value)

    return total


##################################################
# L2 NORM
##################################################

def l2_norm(vector):

    total = 0

    for value in vector:
        total += value ** 2

    return np.sqrt(total)


##################################################
# MAIN
##################################################

def main():

    ##################################################
    # PART 1 : SIGMOID VISUALIZATION
    ##################################################

    x = np.linspace(-10, 10, 100)

    y_sigmoid = sigmoid(x)

    plt.figure(figsize=(8,5))
    plt.plot(x, y_sigmoid)
    plt.title("Sigmoid Function")
    plt.xlabel("x")
    plt.ylabel("sigmoid(x)")
    plt.grid()
    plt.show()


    ##################################################
    # PART 2 : SIGMOID DERIVATIVE VISUALIZATION
    ##################################################

    y_derivative = sigmoid_derivative(x)

    plt.figure(figsize=(8,5))
    plt.plot(x, y_derivative)
    plt.title("Derivative of Sigmoid")
    plt.xlabel("x")
    plt.ylabel("sigmoid'(x)")
    plt.grid()
    plt.show()


    ##################################################
    # PART 3 : L1 AND L2 NORM
    ##################################################

    vector = np.array([1, -2, 3, -4, 5])

    print("\nVector")
    print(vector)

    print("\nL1 Norm")
    print(l1_norm(vector))

    print("\nL2 Norm")
    print(l2_norm(vector))


    ##################################################
    # PART 4 : LOAD WISCONSIN DATASET
    ##################################################

    data = load_breast_cancer()

    X = data.data
    y = data.target

    print("\nDataset Shape")
    print(X.shape)

    print("\nTarget Classes")
    print(np.unique(y))


    ##################################################
    # FEATURE SCALING
    ##################################################

    scaler = StandardScaler()

    X = scaler.fit_transform(X)


    ##################################################
    # TRAIN TEST SPLIT
    ##################################################

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.20,
        random_state=999
    )

    print("\nTraining Samples :", len(X_train))
    print("Testing Samples  :", len(X_test))


    ##################################################
    # PART 5 : LOGISTIC REGRESSION
    ##################################################

    print("\n================================")
    print("LOGISTIC REGRESSION")
    print("================================")

    logistic_model = LogisticRegression(
        max_iter=5000
    )

    logistic_model.fit(
        X_train,
        y_train
    )

    y_pred_logistic = logistic_model.predict(
        X_test
    )

    logistic_accuracy = accuracy_score(
        y_test,
        y_pred_logistic
    )

    print(
        f"Accuracy = {logistic_accuracy:.4f}"
    )

    print("\nConfusion Matrix")
    print(
        confusion_matrix(
            y_test,
            y_pred_logistic
        )
    )


    ##################################################
    # PART 6 : RIDGE CLASSIFIER
    ##################################################

    print("\n================================")
    print("RIDGE CLASSIFIER")
    print("================================")

    ridge_model = RidgeClassifier()

    ridge_model.fit(
        X_train,
        y_train
    )

    y_pred_ridge = ridge_model.predict(
        X_test
    )

    ridge_accuracy = accuracy_score(
        y_test,
        y_pred_ridge
    )

    print(
        f"Accuracy = {ridge_accuracy:.4f}"
    )


    ##################################################
    # PART 7 : LASSO CLASSIFIER
    ##################################################

    print("\n================================")
    print("LASSO CLASSIFIER")
    print("================================")

    lasso_model = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        max_iter=5000
    )

    lasso_model.fit(
        X_train,
        y_train
    )

    y_pred_lasso = lasso_model.predict(
        X_test
    )

    lasso_accuracy = accuracy_score(
        y_test,
        y_pred_lasso
    )

    print(
        f"Accuracy = {lasso_accuracy:.4f}"
    )


    ##################################################
    # PERFORMANCE COMPARISON
    ##################################################

    print("\n================================")
    print("MODEL COMPARISON")
    print("================================")

    print(
        f"Logistic Regression : {logistic_accuracy:.4f}"
    )

    print(
        f"Ridge Classifier    : {ridge_accuracy:.4f}"
    )

    print(
        f"Lasso Classifier    : {lasso_accuracy:.4f}"
    )

    print("\nClassification Report")
    print(
        classification_report(
            y_test,
            y_pred_logistic
        )
    )


if __name__ == "__main__":
    main()