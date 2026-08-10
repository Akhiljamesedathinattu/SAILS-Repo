#Data Standardization (From Scratch)


import numpy as np

def standardize():
    X = np.array([[1, 20],
                  [2, 30],
                  [3, 40],
                  [4, 50]], dtype=float)

    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)

    X_std = (X - mean) / std

    print("Original Array:")
    print(X)

    print("\nStandardized Array:")
    print(X_std)


def mymain():
    standardize()


if __name__ == "__main__":
    mymain()