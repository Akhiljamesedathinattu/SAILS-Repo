def individual():
    num=int(input("enter the number"))
    while(num>0):
        digit=num%2
        print(str(digit))
        print(digit[::-1])
        num=num//2
    
def mymain():
     individual()
     
     
if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
