#K-Fold Cross Validation (From Scratch)
import numpy as np

def k_fold_split():
    X = np.array([[1, 20],
                  [2, 30],
                  [3, 40],
                  [4, 50],
                  [5, 60],
                  [6, 70],
                  [7, 80],
                  [8, 90],
                  [9, 100],
                  [10, 110]])

    y = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1])

    k = 5

    indices = np.arange(len(X))
    np.random.shuffle(indices)

    fold_size = len(X) // k
    print(y)
    print(indices)



    print("K-Fold Splits:\n")

    for i in range(k):
        start = i * fold_size

        if i == k - 1:
            end = len(X)
        else:
            end = start + fold_size

        test_idx = indices[start:end]
        train_idx = np.concatenate((indices[:start], indices[end:]))

        print("Fold", i + 1)
        print("Training Indices:", train_idx)
        print("Testing Indices :", test_idx)
        print()


def mymain():
    k_fold_split()


if __name__ == "__main__":
    mymain()