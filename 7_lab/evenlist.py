def evenlist():
    n=[5,5,8,5,7,9,5,6,9,7,7,5,7,5,7,5,5,5,7,5]
    for i in range(0,len(n)):
            if n[i]%2==0:
                print(n[i],"even")
            else:
                print(n[i],"odd")

    


def mymain():
    evenlist()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
