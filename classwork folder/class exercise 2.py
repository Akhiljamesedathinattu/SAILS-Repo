def classexercise2():
# dictplayers={
#     "rohit":78,
#     "gill":80,
#     "kohli":45,
#     "shyam":56,
# }
#
# "rohit"in dictplayers
# dictplayers["gill"]
# dictplayers["shyam"]
# dictplayers.get("gill")
# dictplayers.get("shyam",'not present')
#
# k=dictplayers.keys()
# print(list(dictplayers.keys()))
# print(list(dictplayers.values()))
# print(list(dictplayers.items()))
#
# empty_set=set()
# even_numbers={0,2,4,6,8}
# new_set=set("letters")
# print(new_set)
# print(set())
# lset=set(["virat","rohit","shubham"])
# print(lset)

# playerscores={
#     "rohit":{53,87,96},
#     "virat":{123,21,78},
#     "shubham":{23,44,85},
# }
# for name,scores in playerscores.items():
#      if 123 in scores:
#          print(name)

# for name,scores in playerscores.items():
#     for v in scores:
#         if v >=50:
#             print(name)
#             # print(type(scores))
#             # print(scores)
#             # print(type(name))
#             # print(name)
#

#     # days=['monday','tuesday','wednesday']
#     # fruits=['orange','banana','peach']
#     # drinks=['coffee','tea','milk']
#     # desserts=['ice cream','pie','tiramisu','pudding']
#     # for day,fruit,drink,dessert in zip(days,fruits,drinks,desserts):
#     #     print(day,":drink",drink,"-eat",fruit,"-enjoy",dessert)
#     print("helloworld")
#
# def comprehensions():
# #     number_list=[number for number in range(1,6)]
# #     print(number_list)
# #     number_list=[number-1  for number in range (1,6)]
# #     print(number_list)
# #     a_list=[number for number in range (1,6) if number % 2 == 1]
# #     print(a_list)
#
#     rows=range(1,4)
#     cols=range(1,3)
#     for row in rows:
#         for col in cols:
#             print(row,col)
#     cells=[(rows,col)for row in rows for col in cols]
#     for cell in cells:
#         print(cell)
#     word='letters'
#     letter_counts={letter:word.count(letter) for letter in word}
#     print(letter_counts)
#
#
#
# #def players(bat, bowl, keep):
# def get_price(amount, tax=0.07, discount=0):
#     total = amount + (amount * tax) - discount
#     return(total)
# def print_kwargs(**kwargs):
#     print('Keyword arguments:', kwargs)
#     print(kwargs.keys())


def text():
    text='''python is a good general-purpose,high level language its design make it readable '''
    print(text)
    print(len(text))
    fout = open('akhil.txt', 'wt')
    bytes_written = fout.write(text)
    print(bytes_written)
    fout.close()
    fin=open('akhil.txt','rt')
    ftext = fin.read()
    fin.close()
    print(len(ftext))









def mymain():
    text()
    classexercise2()
    #comprehensions()
    #positional arguments
    # b=("sanju")
    # bo=("sarah")
    # ke=("san")
    #players(b, bo, ke)
    # keyword arguments
    #print (players(bat='Rohit', bowl='bumrah', keep='ishan'))
    #default argument
    # print(greet("Alice"))
    # print(greet("Bob", "Good morning"))
    # print(get_price(int(input("Enter your price: "))))
    # print(get_price(100, 0.10))
    # print(get_price(100, 0.05, 10))
    # print_kwargs(batsman='Rohit', bowler='bumrah', keeper='ishan')




if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()