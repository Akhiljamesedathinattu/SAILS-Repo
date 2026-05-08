def duplicatelist():
    n =[1, 2, 3, 4, 2, 7, 8, 8, 3, 9, 1]
    l1=[]
    for i in n:
        if n.count(i)>1:
            l1.append(n[i-1])
    print(set(l1))





def mymain():
    duplicatelist()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
