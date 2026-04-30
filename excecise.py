def listops():
    birthday = '1/6/1952'
    blist = birthday.split('/')
    print(blist)
#
def ano1():
    players = ['Sachin', 'Dhoni', 'Virat']
    print(players[0], players[1], players[2])
    print( players[-1],  players[-2],  players[-3])
def birds():
    small_birds = ['hummingbird', 'finch']
    extinct_birds = ['dodo', 'passenger pigeon', 'Norwegian Blue']
    carol_birds = [3, 'French hens', 2, 'turtledoves']
    all_birds = [small_birds, extinct_birds, 'macaw', carol_birds]
    print(all_birds)
    print( all_birds[0])
    print(all_birds[1][0])
def individual():
    num = int(input("enter the number: "))
    a = ''

    while num > 0:
        digit = num % 10
        num = num // 10
        a = a + str(digit)

    v = a[::-1]

    for i in v:
        print(i)


def mymain():
    individual()
    
 def cou():
    s="we are the people"
    
    v=" are you happy"
    f=s.count(v)
    print(v)   

def compress():
s = "aabcecccaaa"
result = ""
count = 1

for i in range(len(s)):
    if i < len(s)-1:
        if s[i] == s[i+1]:
            count += 1
        else:
            result = result + s[i] + str(count)
            count = 1
    else:
        result = result + s[i] + str(count)

print(result)






        
def mymain():
   #listops()
   #ano1()
   #birds()
   #individual()
   cou()
    
    

if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
