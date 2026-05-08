def ktimes():
    n = [10,5,30,40,150,60,60,10,10]
    k=1
    l1=[]
    for i in n:
        if n.count(i)>=k:
            l1.append(i)
    print(set(l1))




def mymain():
    ktimes()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
