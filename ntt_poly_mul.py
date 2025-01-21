import time
import numpy as np
import math
from ntt import ntt_func as fntt
from nttcore import find_modulus
from nttcore import find_primitive_root
from nttcore import reciprocal
from generate_w import gen_w
import logging
from prettytable import PrettyTable
import os

"""
n: polynomial number of coefficients
polynomial_q: polynomial coefficients modulus
phi: phi parameter used to weight polynomial coefficients
"""
n = 256
polynomial_q = 1049089
phi = 2016

"""
setting up print parameters.
ntt_x_print_stage: for x polynomial, it prints specified NTT stage number coefficients and results
ntt_y_print_stage: for y polynomial, it prints specified NTT stage number coefficients and results
inv_ntt_x_print_stage: it prints specified inverse NTT stage number coefficients and results
"""
ntt_x_print_stage = 7
ntt_y_print_stage = 5
inv_ntt_x_print_stage = 7

"""'
setting up logging 
"""
formatter = logging.Formatter("  %(funcName)s -  %(levelname)s: %(message)s")
logger = logging.getLogger("my_logger")
logging.getLogger("").setLevel(logging.DEBUG)
fh = logging.FileHandler("ntt_poly_mul.log")
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logging.getLogger("").addHandler(fh)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logging.getLogger("").addHandler(ch)

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


def fixlen(inp, size):
    inp = inp.replace("[", "")
    inp = inp.replace("]", "")
    length = len(inp)
    if length < size:
        string_val = " " * (size - length)
        new_string = inp + string_val
        return new_string


"""
Find NTT modulus and primitive root
"""
ntt_modulus = find_modulus(n, polynomial_q)
prim_root = find_primitive_root(n, ntt_modulus - 1, ntt_modulus)
prim_root_inv = reciprocal(prim_root, ntt_modulus)
"""
Finding the word length required for RTL implementation of NTT modulus
Finding the number of bits required for NTT weights
"""
ntt_word_length = math.ceil(math.log2(ntt_modulus))
ad_num_of_bits = math.ceil(math.log2(n))

"""
Generating two random test polynomials to be multiplied
"""
x = np.random.randint(polynomial_q, size=(n))
y = np.random.randint(polynomial_q, size=(n))

"""
Check if mem_files folder exists, otherwise create it
"""
check_and_create_folder("./mem_files")

"""
Convert the generated X polynomial coefficients to binary and write them to a .mem file to be loaded into FPGA
"""
filexbin = open("mem_files/xram.mem", "w")
for i in range(0, n):
    binc = bin(x[i])[2:].zfill(ntt_word_length)
    L = str(binc) + "\n"
    filexbin.writelines(L)

"""
Convert the generated Y polynomial coefficients to binary and write them to a .mem file to be loaded into FPGA
"""
fileybin = open("mem_files/yram.mem", "w")
for i in range(0, n):
    binc = bin(y[i])[2:].zfill(ntt_word_length)
    L = str(binc) + "\n"
    fileybin.writelines(L)

"""
write x and y polynomial coefficients into a text file for debug
"""
check_and_create_folder("./debug")
file_x = open("debug/x.txt", "w")
file_y = open("debug/y.txt", "w")
for i in range(0, n):
    L = str(x[i]) + "\n"
    file_x.writelines(L)
    L = str(y[i]) + "\n"
    file_y.writelines(L)

"""
reference result for polynomial multiplication
Here two polynomials are multiplied, reminder of each coefficient modulus polynomial q is calculated and then
the resulted polynomial is divided by reduction polynomial (X^(n-1) + 1)
Result is used as a reference to verify NTT polynomial multiplication
"""
reduction_polynomial = [1] + ((n - 1) * [0]) + [1]
multi_res = np.polymul(np.flip(x, 0), np.flip(y, 0)) % polynomial_q
refntt_mul_res = (
    (np.polydiv(multi_res, reduction_polynomial)[1]) % polynomial_q
).astype(int)
refntt_mul_res = np.flip(refntt_mul_res, 0)
filepmr = open("debug/numpy_poly_ntt_poly_mul_result.txt", "w")
for i in range(0, len(refntt_mul_res)):
    l = str(refntt_mul_res[i]) + "\n"
    filepmr.writelines(l)

"""
Generating an array of phi values
"""
phi_array = np.array([0] * n)
for i in range(0, n):
    phi_array[i] = pow(phi, i, ntt_modulus)

"""
converting array of phy values to binary and writing to a mem file for loading into FPGA
"""
filex_phi_mem = open("mem_files/phi.mem", "w")
for i in range(0, len(phi_array)):
    binc = bin(phi_array[i])[2:].zfill(ntt_word_length)
    L = str(binc) + "\n"
    filex_phi_mem.writelines(L)

"""
writing phi array values to text file for debug
"""
filephi = open("debug/phi.txt", "w")
for i in range(0, n):
    l = str(phi_array[i]) + "\n"
    filephi.writelines(l)

"""
multiplying phi array and x and y coefficients to create a weighted array
"""
weighted_x = np.array([0] * n)
weighted_y = np.array([0] * n)
for i in range(0, n):
    weighted_x[i] = ((x[i] % ntt_modulus) * phi_array[i]) % ntt_modulus
    weighted_y[i] = ((y[i] % ntt_modulus) * phi_array[i]) % ntt_modulus

"""
writing weighted x and y array values to text file for debug
"""

file_ntt_x_phi = open("debug/weighted_x.txt", "w")
file_ntt_y_phi = open("debug/weighted_y.txt", "w")
for i in range(0, n):
    l = str(weighted_x[i]) + "\n"
    file_ntt_x_phi.writelines(l)
    l = str(weighted_y[i]) + "\n"
    file_ntt_y_phi.writelines(l)

"""
Generating twiddling factors
"""
gen_w(n, ntt_modulus, ntt_word_length, ntt_modulus, prim_root, logging)

"""
Taking NTT of weighted X and weighted Y
"""

NX = fntt(weighted_x, "ntt_x", ntt_modulus, prim_root, "fw", logging)
NY = fntt(weighted_y, "ntt_y", ntt_modulus, prim_root, "fw", logging)

"""
Writing NTT of weighted X and weighted Y to file for debug
"""
file_ntt_x = open("debug/nnt_of_weighted_x.txt", "w")
file_ntt_y = open("debug/nnt_of_weighted_y.txt", "w")
for i in range(0, n):
    l = str(NX[i]) + "\n"
    file_ntt_x.writelines(l)
    l = str(NY[i]) + "\n"
    file_ntt_y.writelines(l)

"""
Element wise multiplication NTT(weighted x)* NTT(weighted x)
"""
nx_ny_mul = np.array([0] * n)
for i in range(0, n):
    nx_ny_mul[i] = (NX[i].item() * NY[i].item()) % ntt_modulus

"""
Writing Element wise multiplication result to file for debug
"""
file_nx_ny_mul = open("debug/point_mul_ntt_x_y_phi.txt", "a")
for i in range(0, len(nx_ny_mul)):
    line = str(nx_ny_mul[i])
    line = line.strip("[")
    line = line.strip("]")
    line = line.strip(" ")
    file_nx_ny_mul.writelines(line + "\n")

"""
Taking inverse NTT from the element wise multiplication
"""
mul_inv = fntt(nx_ny_mul, "ntt", ntt_modulus, prim_root_inv, "inv", logging)

filemp = open("debug/inv_nnt_of_point_mul.txt", "a")
for i in range(0, n):
    l = str(mul_inv[i]) + "\n"
    filemp.writelines(l)

"""
Calculating the value of phi inverse. This is used to create inverse phi array and remove the weight initially a added 
to the polynomial to avoid zero padding
"""
phi_inverse = 0
for i in range(0, ntt_modulus):
    if (phi * i) % ntt_modulus == 1:
        phi_inverse = i
        break
if phi_inverse == 0:
    logging.error("Phi inverse not found")
    exit()

"""
Creating phi inverse array
"""
phi_inverse_array = [0] * n
for i in range(0, n):
    phi_inverse_array[i] = pow(phi_inverse, i, ntt_modulus)

"""
Writing Element wise multiplication result to file for debug
"""
filephinv = open("debug/phi_inv.txt", "w")
for i in range(0, n):
    l = str(phi_inverse_array[i]) + "\n"
    filephinv.writelines(l)

"""
calculating the value of scale
"""
scaler = reciprocal(n, ntt_modulus)

"""
Multiplying scale into phi inverse to avoid extra calculation required to scale the NTT multiplication result
"""
file_scaled_phi_inv = open("debug/scaled_phi_inv.txt", "w")
sca_phi_inv = np.array([0] * n)
for i in range(0, n):
    sca_phi_inv[i] = (phi_inverse_array[i] * scaler) % ntt_modulus
for i in range(0, n):
    l = str(sca_phi_inv[i]) + "\n"
    file_scaled_phi_inv.writelines(l)

"""
converted scaled phi inverse values to binary and write to mem file to load into FPGA
"""
file_sc_inv_phi_mem = open("mem_files/scaled_inv_phi.mem", "w")
for i in range(0, len(phi_array)):
    binc = bin(sca_phi_inv[i])[2:].zfill(ntt_word_length)
    L = str(binc) + "\n"
    file_sc_inv_phi_mem.writelines(L)

"""
Writing scaled phi inverse values to file for debug
"""
filephix = open("debug/debug_inv_ntt_point_mul_and_scaled_inv_phi.txt", "a")
for i in range(0, n):
    l = (
        "r "
        + fixlen(str(mul_inv[i]), 10)
        + " phi "
        + fixlen(str(sca_phi_inv[i]), 10)
        + "\n"
    )
    filephix.writelines(l)

"""
multiplying the scaled phi inverse values and inverse NTT results to calculate the final NTT multiplication result
"""
ntt_poly_mul_res = [0] * n
for i in range(0, len(mul_inv)):
    ntt_poly_mul_res[i] = (mul_inv[i] * sca_phi_inv[i]) % ntt_modulus

"""
Writing the final NTT polynomial multiplication result to file for debug
"""
filenttmul = open("debug/python_ntt_ntt_poly_mul_res.txt", "a")
for i in range(0, len(mul_inv)):
    line = str(ntt_poly_mul_res[i])
    line = line.strip("[")
    line = line.strip("]")
    line = line.strip(" ")
    filenttmul.writelines(line + "\n")

"""
Print NTT parameters
"""
columns = ["n", "Poly. q", "NTT modulus", "Prim. Root", "Phi", "Phi inv.", "Scale"]
myTable = PrettyTable()
# Add Columns
myTable.add_column(columns[0], [n])
myTable.add_column(columns[1], [polynomial_q])
myTable.add_column(columns[2], [ntt_modulus])
myTable.add_column(columns[3], [prim_root])
myTable.add_column(columns[4], [phi])
myTable.add_column(columns[5], [phi_inverse])
myTable.add_column(columns[6], [scaler])
logging.info("\n" + str(myTable))

"""
Comparing NTT polynomial multiplication result with reference results
"""
num_errs = 0
for i in range(0, len(mul_inv)):
    if ntt_poly_mul_res[i] != refntt_mul_res[i]:
        num_errs = num_errs + 1
if num_errs == 0:
    logging.info(f" Number of errors: {num_errs}")
else:
    logging.error(f" Number of errors: {num_errs}")
