def comprehension():
    fruits = ['mango', 'kiwi', 'strawberry', 'guava', 'pineapple', 'mandarin orange']

    capitalized_fruits = [fruit.capitalize() for fruit in fruits]

    print(capitalized_fruits)
    
def mymain():
    comprehension()


if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()