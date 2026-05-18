def matrix():
    numbers= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    squares = {num: num**2 for num in numbers if num % 2 == 0}

    print(squares)

def mymain():
    matrix()

    sentence = "Hello, how are you?"

    reverse_words = {word: word[::-1] for word in sentence.split()}

    print(reverse_words)


if __name__ == "__main__":

    mymain()