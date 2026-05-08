def sumavg():
    n=list(input("Enter a list: "))
    v=len(n)
    s1=sum(n)
    s2=sum(n)/v
    print(s1)
    print(s2)


def mymain():
    sumavg()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()

