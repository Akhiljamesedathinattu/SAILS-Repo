def anagram():
    s=str(input("enter the word 1 :"))
    s1=str(input("enter the word 2 :"))
    v=len(s)
    v1=len(s1)
    if v!=v1:
        print("not anagram")
    
    for i in s1:
        if s.count(i)!=s1.count(i):
            print(" not anagram")
            break
        else:
            print(" anagram")
            break
    
        
        























    
def mymain():
    anagram()
if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
        
