def matrix():
    numbers= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    squares = {num: num**2 for num in numbers if num % 2 == 0}

    print(squares)

def mymain():
    matrix()


if __name__ == "__main__":

    mymain()