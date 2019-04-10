import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    nums = []
    with open("output_mastiff/heatmap_out.ibin", "rb") as file:
        data = file.read(4)
        while data:
            num = int.from_bytes(data, "little")
            nums.append(num)
            data = file.read(4)
    print(len(nums))

    a = np.array(nums).reshape((1024, 1024))
    plt.imshow(a, cmap='hot', interpolation='nearest')
    plt.show()

    universe = []
    with open("output_mastiff/uni_out.usbin", "rb") as file:
        data = file.read(2)
        while data:
            cell = int.from_bytes(data, "little")
            universe.append(cell)
            data = file.read(2)
    b = np.array(universe).reshape((4096, 4096))
    print(b)
    plt.imshow(b, cmap='gray')
    plt.show()
