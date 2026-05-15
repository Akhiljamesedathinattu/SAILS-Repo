def most_unique():
    test_dict = {"Gfg" : [5, 7, 7, 7, 7],"Best" : [9, 9, 6, 5, 5], "is" : [6, 7, 7, 8,9,3,4,7], }
    #print(test_dict)
    c=0
    b=''
    for key, value in test_dict.items():
        if len(set(value))>c:
            print(len(set(value)))
            c=len(set(value))
            b=key
    print(b)







def mymain():
    most_unique()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()