#Data Normalization (From Scratch)


import numpy as np

def normalize():
    X = np.array([[1, 20],
                  [2, 30],
                  [3, 40],
                  [4, 50]], dtype=float)

    min_val = np.min(X, axis=0)
    max_val = np.max(X, axis=0)

    X_norm = (X - min_val) / (max_val - min_val)

    print("Original Array:")
    print(X)

    print("\nNormalized Array:")
    print(X_norm)


def mymain():
    normalize()


if __name__ == "__main__":
    mymain()