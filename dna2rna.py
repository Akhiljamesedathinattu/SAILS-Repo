def DNA2RNA():
    s1=input("enter the DNA sequence: ")

    s1=s1.upper()
    n=len(s1)
    s2=s1.replace('T','U')
    print(s2)
            
def motif():
    m=input("enter the sequence:")
    n=input("enter the motif :")
    l=len(m)
    l1=len(n)
    p=""
    for i in m:
        if(m[i:n]==n):
            p.append(i)
        print(p)
        
        
        
    
        
        
    
    
    
    
    
    
    
    
    
    
    def mymain():
    #DNA2RNA()
    motif()
    
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
