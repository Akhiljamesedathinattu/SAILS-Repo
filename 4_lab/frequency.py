def frequency():
    string=str(input("enter your string :  "))
    val=0
   
    for i in string:
        s=string.count(i)
        print(i,"is in",s,"times")
        if(s>=val):
            val=s
            new=i
    print("maximum value =",val,"and the character is :",new)
            
                
            

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def mymain():
    frequency()
    
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
