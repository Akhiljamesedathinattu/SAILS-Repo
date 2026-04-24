def motif():
    m=input("enter the sequence:")
    n=input("enter the motif :")
    l=len(m)
    l1=len(n)
    for i in range(0,l-l1+1):
        if(m[i:i+l1]==n):
            print(i)
        
  
        
        
def mymain():
    motif()


if __name__ == "__main__":
    mymain()


