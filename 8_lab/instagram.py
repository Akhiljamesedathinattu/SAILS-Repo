def insta():
    # Input followers
    followers = [
        "alex", "ben", "charlie", "david", "emma", "frank", "grace", "henry",
        "isabella", "jack", "karen", "liam", "mia", "noah", "olivia", "paul",
        "quinn", "ryan", "sophia", "thomas", "uma", "victor", "william",
        "xavier", "yara", "zane", "aarav", "anaya", "aditya", "bhavya",
        "chirag", "deepak", "esha", "farhan", "gautam", "harini", "isha",
        "jatin", "kavya", "lokesh", "manoj", "nandini", "omkar", "priya",
        "rakesh", "sanjana", "tarun", "uday", "varun", "yesha", "zoya"
    ]
    followers=set(followers)

    # Input following
    following = [
        "ben", "charlie", "emma", "frank", "grace", "henry", "jack", "liam",
        "mia", "noah", "olivia", "paul", "ryan", "sophia", "thomas", "victor",
        "william", "xavier", "aarav", "anaya", "aditya", "deepak", "esha",
        "farhan", "gautam", "harini", "isha", "kavya", "lokesh", "manoj",
        "omkar", "priya", "rakesh", "sanjana", "tarun", "varun", "zoya",
        "neha", "rahul", "sneha", "vikas"
    ]
    following=set(following)

    # Mutuals
    mutuals = followers & following

    # Fans
    fans = followers - following

    print("\nMutuals:", mutuals)
    print("Fans:", fans)

def mymain():
    insta()



if __name__ == "__main__":
    # print("this is the beginning of my program")
    mymain()
