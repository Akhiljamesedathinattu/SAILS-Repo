def matrix_subtraction():
    m1 =[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    m2 =[[10, 20, 30], [40, 50, 60], [70, 80, 90]]
    m=[]
    item=[]





    for item1, item2 in zip(m1, m2):
        s = len(item1)
        for i in range(s):
            v = item1[i] - item2[i]
            item.append(v)

    m.append(item)
    print(m)
def mymain():
    matrix_subtraction()



if __name__ == "__main__":
    mymain()