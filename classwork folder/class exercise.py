def classexercise():

    # dictplayers={
    # "rohit":70,
    # "sachin":52,
    # "sanju":152,
    # }
    # dictplayers["ishan"]=90
    # dictplayers["sreyas"]=125
    # dictplayers["rinku"]=130
    # print(dictplayers)
    #
    #
    # 'rohit'in dictplayers
    # k=dictplayers.values()
    # print(type(k))
    # print(list(dictplayers.keys()))
    # print(list(dictplayers.values()))
    # print(list(dictplayers.items()))
    # empty_set=set()
    # even_numbers={0,2,4,6,8,10}
    # print(even_numbers)
    # new_set=set("letters")
    # print(new_set)
    # print(len(new_set))
    # lset=set(["viart","rohit","sachin","sanju"])
    # print(lset)
    # tset=set(("rohit","sachin","sanju"))
    # print(tset)



   # a={1,2}
   # b={2,3}
   # print(a & b)
   # print(a.intersection(b))
   #  #union
   # print(a|b)
   # print(a-b)
   # print(a.difference(b))
   # print(a^b)
   # print(a-b)
   # print(a.symmetric_difference(b))
   # print(a<=b)
   # print(a.issubset(b))

    print("Class exercise")
def comprehension():
    print("Class exercise")
    # number_list=[numbers for numbers in range(1,6)]
    # print(number_list)
    # number_list=[numbers-1 for numbers in range(1,6)]
    # print(number_list)
    # a_list=[numbers**2 for numbers in range(1,6) ]
    # print(a_list)
    # rows=range(1,4)
    # cols=range(1,3)
    # for row in rows:
    #     for col in cols:
    #         print(row,col)
    # cells=[(rows,cols) for row in rows for col in cols]
    # rows=range(1,4)
    # cols=range(1,3)
    # cells=[(rows,cols) for row in rows for col in cols]
    # for cell in cells:
    #     print(cell)
    # for row,col in cells:
    #     print(row,col)
    # word='letter'
    # letter_counts={letter:word.count(letter) for letter in word}
    # print(letter_counts)
def min(a,b):
    c=a-b
    return(c)
def even(x):
    if x%2==0:
        return(x)
    else:
        return(x)
def odd(x):
    if x%2==1:
        return(x)
    else:
        return(x)
def akhilpos(maths,science,english):
    maths=maths.lower()
    science=science.lower()
    english=english.lower()
def akhilkeyword(maths,science,english):

    text = '''
        python is a good general-purpose,high level labguage its design make it readable '''
    print(text)
    print(len(text))
    fout=open('python.txt',wt)
































def mymain():
    classexercise()
    comprehension()

    x=4
    y=3
    o=min(x,y)
    print(o)
    c=even(x)
    z=odd(x)
    print(c,"even")
    print(z,"odd")
    c="rathish"
    b="balu"
    e="sussan"

    print(akhilpos(c,b,e))
    print(akhilkeyword(science=e,maths=b,english=c))
    text='''
    python is a good general-purpose,high level labguage '''
    print(text)
    print(len(text))







if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
