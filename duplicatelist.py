def duplicatelist():
    n = [10,20,30,40,30,20,10,10]
    l1=[]
    for i in range(len(n)):
        if n[i]==n[i+1]:
            l1.append(n[i])
            print(set(l1))
            break




def mymain():
    duplicatelist()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
