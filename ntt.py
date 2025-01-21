import time
from nttcore import *
from prettytable import PrettyTable
import os

# Define debug arrays
debug_array_256 = np.empty((0, 5))
debug_array_128 = np.empty((0, 5))
debug_array_64 = np.empty((0, 5))
debug_array_32 = np.empty((0, 5))
debug_array_16 = np.empty((0, 5))
debug_array_8 = np.empty((0, 5))
debug_array_4 = np.empty((0, 5))
debug_array_2 = np.empty((0, 5))
"""
This file calculates NTT for an input polynomial. It will prints all values of inputs, weights and outputs per stage
which was used for debuging. The output of function NTT of the input polynomial
x: input polynomial
name" either x or y, this is used for prinitng the debug files
polynomial_q: polynomial q
logging: logging functions used to print debug tables
"""


def ntt_func(x, name, ntt_modeulus, prim_root, ntt_mode, logging):
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

    logging.info("----------------------")
    logging.info(f"Forward NTT python run...")

    def ntt256(x, r, polynomial_q):
        global debug_array_256

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt128(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt128(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array_256
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_256 = np.vstack((debug_array_256, new_row))
        return quotient.reshape(N, 1)

    def ntt128(x, r, polynomial_q):
        global debug_array_128

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt64(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt64(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_128 = np.vstack((debug_array_128, new_row))
        return quotient.reshape(N, 1)

    def ntt64(x, r, polynomial_q):
        global debug_array_64
        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt32(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt32(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_64 = np.vstack((debug_array_64, new_row))

        return quotient.reshape(N, 1)

    def ntt32(x, r, polynomial_q):
        global debug_array_32

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt16(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt16(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_32 = np.vstack((debug_array_32, new_row))

        return quotient.reshape(N, 1)

    def ntt16(x, r, polynomial_q):
        global debug_array_16

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt8(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt8(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]
        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_16 = np.vstack((debug_array_16, new_row))
        return quotient.reshape(N, 1)

    def ntt8(x, r, polynomial_q):
        global debug_array_8

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt4(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt4(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_8 = np.vstack((debug_array_8, new_row))
        return quotient.reshape(N, 1)

    def ntt4(x, r, polynomial_q):
        global debug_array_4

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]

        # Split input into even and odd indices
        xe = x[::2]
        xo = x[1::2]

        # Find modulus and primitive root for smaller NTT
        length = len(xe)
        inter_ntt_mode = find_modulus(length, polynomial_q)
        inter_ntt_root = find_primitive_root(length, inter_ntt_mode - 1, inter_ntt_mode)
        if ntt_mode == "inv":
            inter_ntt_root = reciprocal(inter_ntt_root, inter_ntt_mode)

        # Recursively compute NTT for even and odd parts
        X_even = ntt2(x[::2], inter_ntt_root, inter_ntt_mode)
        X_odd = ntt2(x[1::2], inter_ntt_root, inter_ntt_mode)

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_4 = np.vstack((debug_array_4, new_row))
        return quotient.reshape(N, 1)

    def ntt2(x, r, polynomial_q):
        global debug_array_2

        # Convert input to numpy array of integers
        x = np.asarray(x, dtype=int)
        N = x.shape[0]
        n = np.arange(N)

        # Recursively compute NTT for even and odd parts
        X_even = x[0]
        X_odd = x[1]

        # Calculate twiddle factors for next step
        TS = [pow(int(r), k, int(polynomial_q)) for k in range(N // 2)]
        T = [
            pow(int(r), k, int(polynomial_q)) * X_odd[k] % polynomial_q
            for k in range(N // 2)
        ]

        # Calculate even and odd parts
        sme1 = [(X_even[k] + T[k]) % polynomial_q for k in range(N // 2)]
        sme2 = [(X_even[k] - T[k]) % polynomial_q for k in range(N // 2)]

        # Concatenate results
        quotient = np.concatenate([sme1, sme2])

        # Debugging: Store intermediate values in debug_array
        T = [pow(int(r), k, int(polynomial_q)) % polynomial_q for k in range(N // 2)]
        for i in range(0, N // 2):
            new_row = np.array(
                [
                    [
                        str(X_odd[i]).replace("[", "").replace("]", ""),
                        str(X_even[i]).replace("[", "").replace("]", ""),
                        str(T[i]).replace("[", "").replace("]", ""),
                        str(sme1[i]).replace("[", "").replace("]", ""),
                        str(sme2[i]).replace("[", "").replace("]", ""),
                    ]
                ]
            )
            debug_array_2 = np.vstack((debug_array_2, new_row))
        return quotient.reshape(N, 1)

    # Get the length of the input vector
    n = len(x)

    # Reshape the input vector to a column vector
    x = x.reshape((n, 1))

    # Set up the folder path for debug files
    folder_path = "debug/fw_" + name
    check_and_create_folder(folder_path)
    file_table_256 = open(folder_path + "/" + name + "_stage_7.txt", "w")
    file_table_128 = open(folder_path + "/" + name + "_stage_6.txt", "w")
    file_table_64 = open(folder_path + "/" + name + "_stage_5.txt", "w")
    file_table_32 = open(folder_path + "/" + name + "_stage_4.txt", "w")
    file_table_16 = open(folder_path + "/" + name + "_stage_3.txt", "w")
    file_table_8 = open(folder_path + "/" + name + "_stage_2.txt", "w")
    file_table_4 = open(folder_path + "/" + name + "_stage_1.txt", "w")
    file_table_2 = open(folder_path + "/" + name + "_stage_0.txt", "w")

    # Choose the appropriate NTT function based on the input size
    if n == 256:
        NTT_of_X = ntt256(x, prim_root, ntt_modeulus)
    elif n == 128:
        NTT_of_X = ntt128(x, prim_root, ntt_modeulus)
    elif n == 64:
        NTT_of_X = ntt64(x, prim_root, ntt_modeulus)
    elif n == 32:
        NTT_of_X = ntt32(x, prim_root, ntt_modeulus)
    elif n == 16:
        NTT_of_X = ntt16(x, prim_root, ntt_modeulus)
    elif n == 8:
        NTT_of_X = ntt8(x, prim_root, ntt_modeulus)
    elif n == 4:
        NTT_of_X = ntt4(x, prim_root, ntt_modeulus)
    elif n == 2:
        NTT_of_X = ntt2(x, prim_root, ntt_modeulus)
    else:
        logging.error("NTT Error, wrong input size")
        exit()

    # Define column names for debug tables
    debug_columns = ["X_odd", "X_even", "w", "out_1", "out_2"]

    # Function to write debug information to a file
    def write_debug_table_to_file(debug_Table, debug_array, targ_file):
        debug_Table.add_column(debug_columns[0], debug_array[:, 0])
        debug_Table.add_column(debug_columns[1], debug_array[:, 1])
        debug_Table.add_column(debug_columns[2], debug_array[:, 2])
        debug_Table.add_column(debug_columns[3], debug_array[:, 3])
        debug_Table.add_column(debug_columns[4], debug_array[:, 4])
        targ_file.write(debug_Table.get_string())

    # Write debug tables to files for each NTT stage
    debug_Table_256 = PrettyTable()
    write_debug_table_to_file(debug_Table_256, debug_array_256, file_table_256)

    debug_Table_128 = PrettyTable()
    write_debug_table_to_file(debug_Table_128, debug_array_128, file_table_128)

    debug_Table_64 = PrettyTable()
    write_debug_table_to_file(debug_Table_64, debug_array_64, file_table_64)

    debug_Table_32 = PrettyTable()
    write_debug_table_to_file(debug_Table_32, debug_array_32, file_table_32)

    debug_Table_16 = PrettyTable()
    write_debug_table_to_file(debug_Table_16, debug_array_16, file_table_16)

    debug_Table_8 = PrettyTable()
    write_debug_table_to_file(debug_Table_8, debug_array_8, file_table_8)

    debug_Table_4 = PrettyTable()
    write_debug_table_to_file(debug_Table_4, debug_array_4, file_table_4)

    debug_Table_2 = PrettyTable()
    write_debug_table_to_file(debug_Table_2, debug_array_2, file_table_2)

    logging.info(f"Done Forward NTT python")
    logging.info("----------------------")
    return NTT_of_X
