def tupleops():
    empty_tuple = ()
    print(empty_tuple)
    persons = 'Rohit', 'Shubham', 'Virat', 'Ishan'
    print(persons)
    a, b, c, d = persons
    print(a, b, c, d)
    a = 5
    b = 6
    a, b = b, a
    print(a, b)
    mylist = ['a', 5, 'c']
    mytuple = tuple(mylist)
    print(mytuple)
def dictops():
    empty_dict = {}
    scores = {"Rahul": 86, "Shubham": 128, "Virat": 110, "Dhoni": 86}
    print(scores)
    print(scores["Shubham"])
    lol = [['a', 'b'], ['c', 'd'], ['e', 'f']]
    print(dict(lol))
    lot = [('a', 'b'), ('c', 'd'), ('e', 'f')]
    print(dict(lot))
    tol = (['a', 'b'], ['c', 'd'], ['e', 'f'])
    print(dict(tol))
    los = ['ab', 'cd', 'ef']
    print(dict(los))
    tos = ('ab', 'cd', 'ef')
    print(dict(tos))
    dictplayers = {
        "Rohit": 78,
        "Gill": 80,
        "Kohli": 45,
    }
    dictplayers["ishan"] = 90
    dictplayers["sanju"] = 125
    dictplayers["dube"] = 130
    print(dictplayers)
    dictnew = {
        "bumrah": 30,
        "siraj": 20
    }
    dictplayers.update(dictnew)
    print(dictplayers)
    del dictplayers["siraj"]
    temp = dictplayers.copy()
    dictplayers.clear()
    dictplayers = temp.copy()
    print(dictplayers)


    dictplayers.keys()
    print(list(dictplayers.keys()))
def setops():
    empty_set = set()
    even_numbers = {0, 2, 4, 6, 8}
    new_set = set("letters")
    print(new_set)
    print(even_numbers)
    playerscores = {
        "Rohit": {53, 87, 96},
        "Virat": {123, 21, 78},
        "Shubham": {23, 44, 85}
    }
    for name, scores in playerscores.items():
        if 123 in scores:
            print(name)
        for v in scores:
            if v >= 50:
                print(name)
                break
    houses = {(44.79, -93.14, 285): 'My House', (38.89, -77.03, 13): 'The White House'}
    print(houses[(44.79, -93.14, 285)])

def swapList():
    newList = [12, 35, 9, 56, 24]
    size = len(newList)

    # Swapping
    temp = newList[0]
    newList[0] = newList[size - 1]
    newList[size - 1] = temp

    print(newList)

    return newList
def matrix_addition(m1, m2):
    result = []
    m1_size = len(m1)

    for i in range(m1_size):
        item1 = m1[i]
        item2 = m2[i]
        row_size = len(item1)
        row = []
        for j in range(row_size):
            s = item1[j] + item2[j]
            row.append(s)

        result.append(row)

    return result
def matrix_sum(m1, m2):
    m, item = [], []

    for item1, item2 in zip(m1, m2):
        s = len(item1)
        item = []
        for i in range(s):
            item = []
            v = item1[i] + item2[i]
            item.append(v)

        m.append(item)

    return m



def mymain():
    # tupleops()
    #dictops()
    # setops()
    # swapList()
    m1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    m2 = [[10, 20, 30], [40, 50, 60], [70, 80, 90]]
    m= matrix_sum(m1, m2)
    m = matrix_addition(m1, m2)
    print(m)
    m = matrix_addition(m1, m2)
    print(m)
































if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()