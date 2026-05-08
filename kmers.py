def kmer():
    s1=str(input("enter sequence for kmers: "))
    v=int(input(' enter the kmer needed;  '))
    n=len(s1)
    
    for i in range(n-v+1):
        s2=s1[i:i+v]
        print(s2)
            
        
        
    
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def mymain():
    kmer()
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
