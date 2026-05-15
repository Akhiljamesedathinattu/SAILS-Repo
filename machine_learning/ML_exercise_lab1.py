import matplotlib as mpl
import matplotlib.pyplot as plt
def transpose():
    A = [
        [1, 2, 3],
        [4, 5, 6]
    ]


    # Find transpose of A
    AT = []

    for j in range(len(A[0])):   # columns of A
        row = []
        for i in range(len(A)):  # rows of A
            row.append(A[i][j])
        AT.append(row)
    print(AT)


    # Compute ATA = AT × A
    ATA = []

    for i in range(len(AT)):
        row = []
        for j in range(len(A[0])):
            sum = 0
            for k in range(len(A)):
                sum += AT[i][k] * A[k][j]
            row.append(sum)
        ATA.append(row)
    print("Matrix A:")
    for row in A:
        print(row)

    print("\nTranspose of A:")
    for row in AT:
        print(row)

    print("\nTranspose of A*A:")
    for row in ATA:
        print(row)
def mymain():
    transpose()
    # Print result



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()