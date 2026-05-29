import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
matplotlib.use('TkAgg')
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

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

#Implement y = 2x1 + 3 and plot x1, y [start=-100, stop=100, num=100]


    x1 = np.linspace(-100, 100, 100)
    y= 2 * x1 + 3


    plt.plot(x1, y)
    plt.title("y = 2x1 + 3")
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()
#Implement y = 2x12 + 3x1 + 4 and plot x1, y in the range [start=--10, stop=10, num=100]
    x1 = np.linspace(-10, 10, 100)
    y = 2 * x1 ** 2 + 3 * x1 + 4

    plt.plot(x1, y)
    plt.title("y = 2x1² + 3x1 + 4")
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()
#Implement Gaussian PDF - mean = 0, sigma = 15 in the range[start=-100, stop=100, num=100]
    mu = 0
    sigma = 15

    x = np.linspace(-100, 100, 100)

    pdf = (1 / (sigma * (2 * 3.14159) ** 0.5)) *(2.71828 ** (-((x - mu) ** 2) / (2 * sigma ** 2)))

    print("PDF =", pdf)
    plt.figure(figsize=(6, 4))
    plt.plot(x, pdf)
    plt.title("Gaussian PDF")
    plt.xlabel("x")
    plt.ylabel("PDF")
    plt.grid(True)
    plt.show()
#Implement y = x1^2, plot x1, y in the range [start=--10, stop=10, num=100]. Compute the value of derivatives at these points, x1 = -5, -3, 0, 3, 5.  What is the value of x1 at which the function value (y) is zero. What do you infer from this?
    x1 = np.linspace(-10, 10, 100)
    y = x1 ** 2

    plt.figure(figsize=(6, 4))
    plt.plot(x1, y)
    plt.title("y = x1²")
    plt.xlabel("x1")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()

    points = [-5, -3, 0, 3, 5]

    print("\nDerivatives:")
    for p in points:
        derivative = 2 * p
        print(f"At x = {p}, derivative = {derivative}")
    print("\nFunction y = x² is zero at x = 0")
#Implement y = y=2x1+3x2^2+3x3^3+4, where x1, x2 and x3 are three independent variables. Compute the gradient of y at a few points and print the values.


def gradient():
    x1 = int(input("Enter x1: "))
    x2 = int(input("Enter x2: "))
    x3 = int(input("Enter x3: "))


    # Function value
    y = 2 * x1 + 3 * (x2 ** 2) + 3 * (x3 ** 3) + 4

    #Gradient calculations
    dy_dx1 = 2
    dy_dx2 = 6 * x2
    dy_dx3 = 9 * x3

    print(f"\nFor point ({x1}, {x2}, {x3})")
    print(f"y = {y}")

    print("\nGradient of y:")
    print([dy_dx1, dy_dx2, dy_dx3])


    #7.Linearmodel: Xθ
    theta = [2, 3, 3]

    X1 = [
        [1, 0, 2],
        [0, 1, 1],
        [2, 1, 0],
        [1, 1, 1],
        [0, 2, 1]
    ]

    Xθ = []

    for i in range(len(X1)):

        sum = 0

        for j in range(len(theta)):
            sum += X1[i][j] * theta[j]

        Xθ.append(sum)

    print("theta:")
    for row in theta:
        print(row)

    print("\nX1:")
    for row in X1:
        print(row)

    print("\nXθ:")
    for row in Xθ:
        print("|",row,"|")


    # theta = np.array([2,3,3])
    #
    # X1 = np.array([[1, 0, 2],[0, 1, 1], [2, 1, 0],[1, 1, 1],[0, 2, 1]])
    #
    # result = X1 @ theta
    #
    # print("\nXθ =")
    #
    # for row in result:
    #     print("⥏",row,"⥑")




def mymain():
    transpose()
    gradient()
    # [X, y] = load_data()
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=10)
    # print('----------TRAINING----------')
    # print("N=%d" % (len(X)))
    # # training  LINEAR REGRESSION
    # model = LinearRegression()
    # # train the model
    # model.fit(X_train, y_train)
    # # prediction on a test set
    # y_pred = model.predict(X_test)
    # # perfomance measure   r2 score  this will give us a number if it close to 1 the model is good
    # r2 = r2_score(y_test, y_pred)
    # print("r2 score is %0.2f (closer to 1 is good)" % r2)
    # print("done")

    # Print result



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()