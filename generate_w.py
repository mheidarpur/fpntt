import time
import math
import os
import numpy as np
from nttcore import find_modulus
from nttcore import find_primitive_root
from nttcore import reciprocal


def gen_w(n, mod, data_num_of_bits, ntt_modulus, prim_root, logging):
    """
    Check if a  folder exists, otherwise create it
    """

    def check_and_create_folder(folder_path):
        # Check if the folder exists
        if not os.path.exists(folder_path):
            # If it doesn't exist, create it
            try:
                os.makedirs(folder_path)
                logging.info(f"Folder created successfully: {folder_path}")
            except OSError as e:
                logging.error(f"Error creating folder: {e}")
        else:
            logging.info(f"Folder already exists: {folder_path}")

    check_and_create_folder("./w/fwd/bin")
    check_and_create_folder("./w/fwd/dec")
    check_and_create_folder("./w/inv/bin")
    check_and_create_folder("./w/inv/dec")

    w = prim_root
    w_array_length = int(n / 2)
    number_of_ntt_stages = int(math.log(n, 2))

    """
        generating forward ntt w
    """
    w_array = np.array([0] * (w_array_length))
    for i in range(0, w_array_length):
        w_array[i] = pow(w, i, mod)

    file_w = open("debug/w.txt", "w")
    for i in range(0, w_array_length):
        l = str(w_array[i]) + "\n"
        file_w.writelines(l)

    num_addr_bits = int(math.log(w_array_length, 2))
    # creating  array of w for each stage
    for j in range(0, number_of_ntt_stages):
        file_w = open("./w/fwd/dec/w_stage" + str(j) + ".txt", "w")
        file_w_bin = open("./w/fwd/bin/w_stage" + str(j) + ".txt", "w")
        for i in range(0, int(pow(2, j))):
            bits = num_addr_bits
            ry = reverse(i, bits)
            w_reordered = w_array[int(ry)]

            L = str(w_reordered) + "\n"
            file_w.writelines(L)

            bin_w = bin(w_reordered)[2:].zfill(data_num_of_bits)
            BL = str(bin_w) + "\n"
            file_w_bin.writelines(BL)
        file_w.close()
        file_w_bin.close()

    """
        inverse ntt w
    """
    w_inv = reciprocal(prim_root, ntt_modulus)
    logging.info(f"w inverse = {w_inv}")

    """
    Create an array of all w inverse values
    """
    w_inv_array = np.array([0] * (w_array_length))
    for i in range(0, w_array_length):
        w_inv_array[i] = pow(w_inv, i, mod)

    file_w_inv = open("debug/w_inv.txt", "w")
    for i in range(0, w_array_length):
        l = str(w_inv_array[i]) + "\n"
        file_w_inv.writelines(l)

    """
    creating  array of w inverse for each stage
    """
    for j in range(0, number_of_ntt_stages):
        file_w_inv = open("./w/inv/dec/inv_w_stage" + str(j) + ".txt", "w")
        file_w_inv_bin = open("./w/inv/bin/w_inv_stage" + str(j) + ".txt", "w")
        k = int(n / (pow(2, j + 1)))
        order_of_read = 0
        for i in range(0, int(pow(2, j))):
            bits = num_addr_bits
            w_inv_reordered = w_inv_array[int(order_of_read)]
            order_of_read = order_of_read + k
            L = str(w_inv_reordered) + "\n"
            file_w_inv.writelines(L)

            bin_w_inv = bin(w_inv_reordered)[2:].zfill(data_num_of_bits)
            IL = str(bin_w_inv) + "\n"
            file_w_inv_bin.writelines(IL)
        file_w_inv.close()
        file_w_inv_bin.close()


def reverse(x, bits):
    y = 0
    for i in range(bits):
        y = (y << 1) | (x & 1)
        x >>= 1
    return y
