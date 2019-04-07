import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    nums = []
    with open("uni_out.ibin", "rb") as file:
        data = file.read(4)
        while data:
            num = int.from_bytes(data, "little")
            nums.append(num)
            data = file.read(4)
    print(len(nums))
    print(nums)

    a = np.array(nums).reshape((128, 128))
    plt.imshow(a, cmap='hot', interpolation='nearest')
    plt.show()