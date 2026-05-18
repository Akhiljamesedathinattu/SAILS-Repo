import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
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


    #  ATA = AT × A
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


    x1 = np.linspace(-100, 100, 100)
    y = 2 * x1 + 3

    plt.plot(x1, y)
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.title("y = 2x1 + 3")
    plt.grid(True)
    plt.show()
    # Print result


    x1 = np.linspace(-10, 10, 100)
    y = 2 * x1 ** 2 + 3 * x1 + 4

    plt.plot(x1, y)
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.title("y = 2x1^2 + 3x1 + 4")
    plt.grid(True)
    plt.show()

    mu = 0
    sigma = 15

    x = np.linspace(-100, 100, 100)

    pdf = (1 / (sigma * np.sqrt(2 * np.pi))) * \
          np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

    plt.plot(x, pdf)
    plt.xlabel("x")
    plt.ylabel("PDF")
    plt.title("Gaussian PDF")
    plt.grid(True)
    plt.show()
    

    x1 = np.linspace(-10, 10, 100)
    y = x1 ** 2

    plt.plot(x1, y)
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.title("y = x1^2")
    plt.grid(True)
    plt.show()

    # Derivatives
    points = [-5, -3, 0, 3, 5]

    for x in points:
        derivative = 2 * x
        print(f"x = {x}, derivative = {derivative}")




if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()