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


if __name__ == "__main__":
    mymain()


