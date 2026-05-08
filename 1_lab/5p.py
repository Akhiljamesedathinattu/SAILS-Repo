def fibinocci():
    n=int(input())
    a=0
    b=1
    z=0
    print("fib series:")
    while(z<=n):
        c=a+b
        a=b
        b=c
        z=z+1
        print(a)
        
        
def mymain():
    fibinocci()

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
    
    
    
