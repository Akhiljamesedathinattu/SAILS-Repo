def compression():
    s=str(input("enter the sentance  :"))
    n=len(s)
    v=''
    c=0
    g=1
    for i in range(0,n-1):
        if s[i]==s[i+1]:
            c=c+2
            sttrg=s[i]+str(c)
            print(sttrg)

        else:
            print(s[i],g)


                
            
    
    
    
        
        























    
def mymain():
    compression()
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
