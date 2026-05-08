def prime():
    n=int(input("enter the number"))
    if (n<=1):
        print("not prime")
        return
    for i in range (2,n):
        if(n%i==0):
            print("not prime")
            return
    print("prime")
        
def mymain():
     prime()
    
if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()

