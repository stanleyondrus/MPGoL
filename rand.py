import random

def read_nums(pathname):
    with open(pathname, "r") as file:
        nums = []
        for line in file:
            line = line.replace("\n", "").split()
            num = float(line[0])
            nums.append(num)
        return nums

def add_error(num, err):
    parity = random.choice([-1, 1])
    scale = random.random()
    err_n = err * num * scale * parity
    return num + err_n

if __name__ == "__main__":
    nums = read_nums("nums.txt")
    thresh = 0.05
    new_nums = [add_error(num, thresh) for num in nums]
    print(nums)
    print(new_nums)