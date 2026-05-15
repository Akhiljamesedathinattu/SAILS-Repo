def text():
    # fin = open("akhil.txt", "rt")
    # ftext = fin.read()
    # fin.close()
    # print(ftext)




    #read a line at time
    ftext=''
    fin = open("akhil.txt", "rt")
    while True:
        line = fin.readline()
        if not line:
            break
        ftext+=line
    fin.close()
    print(ftext)
    print(len(ftext))






def mymain():
    text()
if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()