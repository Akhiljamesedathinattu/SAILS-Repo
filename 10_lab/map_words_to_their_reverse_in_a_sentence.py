def reverse():
    sentence = "Hello, how are you?"

    reverse_words = {word: word[::-1] for word in sentence.split()}

    print(reverse_words)

def mymain():
    reverse()


if __name__ == "__main__":

    mymain()
