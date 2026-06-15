
def comprehension():

        fruits = ['mango', 'kiwi', 'strawberry', 'guava', 'pineapple', 'mandarin orange']

        vowels = "aeiou"

        fruits_with_only_two_vowels = [
            fruit
            for fruit in fruits
            if sum(1 for char in fruit if char in vowels) == 2
        ]

        print(fruits_with_only_two_vowels)

    def mymain():
        comprehension()

    if __name__ == "__main__":
        mymain()


def mymain():
    comprehension()


if __name__ == "__main__":
    mymain()