def most_unique():
    test_dict = {"Gfg" : [5, 7, 7, 7, 7],"Best" : [6,7,7,9,10], "is" : [9,9,5,10,11,5], }
    #print(test_dict)

    c= 0
    result =[]

    for key, value in test_dict.items():
        count = len(set(value))

        if count > c:
            c =count
            result = [key]
        elif c ==count:
                result.append(key)

    print("Key:", result)












def mymain():
    most_unique()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()