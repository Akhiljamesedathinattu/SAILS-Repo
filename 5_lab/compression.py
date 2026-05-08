def compression():
    s ="aabbccddddddeffff"
    c = 1
    e = ''
    for i in range(len(s)-1):
        if s[i] == s[i + 1]:
            c=c+1
        else:
            e = e + s[i] + str(c)
            c = 1

    e = e + s[i] + str(c)
    print(e)




    
        
        























    
def mymain():
    compression()
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
    
