

def Helloworld():
    print ("Hello world")
def celcius_to_fahrehiet():
    f=15487
    c=((9/5)*(f-32))
    return c
def swap():
    a=5
    b=8
    c=a
    d=b
    print("a and b")
    print(a,b)
    print("b and a")
    print(d,c)
    
def fahrehiet_to_celcius():
    c=894
    f=((c-32)*(5/9))
    return f
    
def mymain():
    Helloworld()
    c=celcius_to_fahrehiet()
    f=fahrehiet_to_celcius()
    print(c)
    print(f)
    #swap()
    a=5
    b=8
    z=((a+3),(b-3))
    print(a,b,z)
    
    

    
if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()

