#Logistic Regression on Sonar Dataset (Without Normalization)
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score


def logistic_regression():

    # Read dataset
    df = pd.read_csv("sonar.csv", header=None)

    X = df.iloc[:, :-1].values
    y = df.iloc[:, -1].values

    # Number of folds
    k = 10

    # Create and shuffle indices
    indices = np.arange(len(X))
    np.random.seed(42)
    np.random.shuffle(indices)

    fold_size = len(X) // k

    scores = []

    print("K-Fold Cross Validation\n")

    for i in range(k):

        start = i * fold_size

        if i == k - 1:
            end = len(X)
        else:
            end = start + fold_size

        # Testing indices
        test_idx = indices[start:end]

        # Training indices
        train_idx = np.concatenate((indices[:start], indices[end:]))

        # Split the data
        X_train = X[train_idx]
        X_test = X[test_idx]

        y_train = y[train_idx]
        y_test = y[test_idx]

        # Train Logistic Regression
        model = LogisticRegression(max_iter=1000)

        model.fit(X_train, y_train)

        pred = model.predict(X_test)

        acc = accuracy_score(y_test, pred)

        scores.append(acc)

        print("Fold", i + 1)
        print("Training samples:", len(train_idx))
        print("Testing samples :", len(test_idx))
        print("Accuracy =", acc)
        print()

    print("Average Accuracy =", sum(scores) / len(scores))


def mymain():
    logistic_regression()


if __name__ == "__main__":
    mymain()