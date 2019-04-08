import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    nums = []
    with open("output_q/heatmap_out.ibin", "rb") as file:
        data = file.read(4)
        while data:
            num = int.from_bytes(data, "big")
            nums.append(num)
            data = file.read(4)
    print(len(nums))
    print(nums)

    # universe = []
    # with open("uni_out.usbin", "rb") as file:
    #     data = file.read(2)
    #     while data:
    #         cell = int.from_bytes(data, "little")
    #         universe.append(cell)
    #         data = file.read(2)
    # b = np.array(universe).reshape((4096, 4096))
    # plt.imshow(b, cmap='grey')
    # plt.show()

    a = np.array(nums).reshape((1024, 1024))
    plt.imshow(a, cmap='hot', interpolation='nearest')
    plt.show()