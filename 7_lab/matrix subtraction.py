def matrix_subtraction():
    m1 =[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    m2 =[[10, 20, 30], [40, 50, 60], [70, 80, 90]]
    m=[]

    s = len(m1)
    print("Matrix A:")
    for row in m1:
        print(row)
    print("\n matrix B:")
    for row in m2:
        print(row)





    for item1, item2 in zip(m1, m2):
        item = []
        for i in range(s):
            v = item1[i] - item2[i]
            item.append(v)
        m.append(item)

    print("\nsubstraction A-B:")
    for row in m:
        print(row)


def mymain():
    matrix_subtraction()



if __name__ == "__main__":
    mymain()