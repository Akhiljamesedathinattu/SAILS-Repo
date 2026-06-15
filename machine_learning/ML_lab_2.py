import numpy as np


def covariance():
    X = [
        [1, 0, 2],
        [0, 1, 1],
        [2, 1, 0],
        [1, 1, 1],
        [0, 2, 1]
    ]

    n = len(X)
    m = len(X[0])

    # Find column means
    means = []

    for j in range(m):# this is for column
        total = 0
        for i in range(n):# this is for row
            total= total+ X[i][j]
            mu=total/n
        means.append(mu)
    print(means)
    #i am finding the difference of  the mean with each row
    mu_c = []

    for i in range(n):
        r= []
        for j in range(m):
            value = round(X[i][j] - means[j], 2)
            r.append(value)
        mu_c.append(r)

    print("\nAfter subtracting")
    for row in mu_c:
        print("|", *row, "|")


 # next is the transpose of mu_c
    mu_c_T = []

    for j in range(len( mu_c[0])):  # columns of A
        row=[]
        for i in range(len( mu_c)):  # rows of A
            row.append( mu_c[i][j])
        mu_c_T.append(row)
    print("\nthe transpose")
    for row in mu_c_T:
        print("|", *row, "|")
    # then the multiplication of mu_c_T*mu_c
    COV=[]


    for i in range(len(mu_c_T)):
        row=[]

        for j in range(len(mu_c[0])):
            total = 0
            for k in range(len(mu_c)):
                total += round(mu_c_T[i][k] * mu_c[k][j],2)
            row.append(total/(n-1))

        COV.append(row)

    print("\nTranspose of matric * matrix:")
    for row in COV:
        print("|", *row, "|")


    import numpy as np

    X = np.array([
        [1, 0, 2],
        [0, 1, 1],
        [2, 1, 0],
        [1, 1, 1],
        [0, 2, 1]
    ])

    # Covariance matrix
    cov_matrix = np.cov(X, rowvar=False,ddof=1)
    print("Covariance matrix verification with numpy")
    print(cov_matrix)





#Compute the dot product of two vectors, x and y given below
# x = [2  1  2]T and y = [1  2  2]T . What is the meaning of the dot product of two vectors? Illustrate that with your own example.
def dot_product():

        x = [2, 1, 2]
        y = [1, 2, 2]

        result = 0
        print("\n")

        for i in range(len(x)):
            result += x[i] * y[i]
        print("Vector x:")
        for i in x:
            print(i)

        print("Vector y:")
        for i in y:
            print(i)
        print("Dot Product =", result)




def mymain():
    covariance()
    dot_product()


if __name__ == "__main__":
    mymain()

