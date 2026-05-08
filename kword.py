def duplicatelist():

    n1 = list(input("enter the words:").split())
    print(len(n1))
    r = str(input())
    n = [x for x in n1 if str(x[0])==r]
    #print(n1)

    print(" total number of words with", r ," is ",len(n))
    print(n)


def mymain():
    duplicatelist()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
