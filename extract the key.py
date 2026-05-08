def matrixi():
    test_dict = {
        "Apple": [1, 2, 2, 3, 4],
        "Ball": [5, 5, 6, 7, 7],
        "Cat": [8, 9, 9, 10, 10],
        "Dog": [11, 12, 12, 13],
        "Egg": [14, 14, 15, 16],
        "Fish": [17, 18, 18, 19],
        "Green": [2, 4, 4, 6, 8],
        "Blue": [1, 3, 5, 5, 7],
        "Red": [9, 9, 10, 11],
        "Sun": [20, 21, 21, 22],
        "Moon": [23, 24, 24, 25],
        "Star": [26, 27, 27, 28],
        "Book": [3, 3, 4, 5, 6],
        "Pen": [7, 8, 8, 9],
        "Bag": [10, 10, 11, 12],
        "King": [1, 1, 2, 3],
        "Queen": [4, 5, 5, 6],
        "Prince": [7, 8, 8, 9],
        "Car": [10, 20, 20, 30],
        "Bike": [40, 40, 50, 60]
    }
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
    matrixi()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()