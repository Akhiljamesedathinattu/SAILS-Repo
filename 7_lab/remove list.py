def duplicatelist():
    # n1=[2,5,2,4,2,6,9,5,3,4,6]
    # r=int(input())
    # print(type(r))
    # n=[x for x in n1 if x!=r]
    # print(n)
    n1 = n1=[2,5,2,4,2,6,9,5,3,4,6]
    print(type(n1[1]))
    r = int(input())
    n = [x for x in n1 if int(x) != r]
    print(n1)
    print(n)


def mymain():
    duplicatelist()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
