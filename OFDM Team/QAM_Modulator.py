import matplotlib.pyplot as plt
#16QAM Generator
def mapping(matrix_len4_arrays): #can upgrade to 64QAM, takes in a matrix with arrays all of length 4
    mapping = {
        (0, 0, 0, 0): -3-3j,
        (0, 0, 0, 1): -3-1j,
        (0, 0, 1, 0): -3+3j,
        (0, 0, 1, 1): -3+1j,
        (0, 1, 0, 0): -1-3j,
        (0, 1, 0, 1): -1-1j,
        (0, 1, 1, 0): -1+3j,
        (0, 1, 1, 1): -1+1j,
        (1, 0, 0, 0): 3-3j,
        (1, 0, 0, 1): 3-1j,
        (1, 0, 1, 0): 3+3j,
        (1, 0, 1, 1): 3+1j,
        (1, 1, 0, 0): 1-3j,
        (1, 1, 0, 1): 1-1j,
        (1, 1, 1, 0): 1+3j,
        (1, 1, 1, 1): 1+1j
        }
    return list(mapping[tuple(i)] for i in matrix_len4_arrays)

def demapping(complex_value_array):
    demapping = {
        -3-3j: (0, 0, 0, 0), 
        -3-1j: (0, 0, 0, 1), 
        -3+3j: (0, 0, 1, 0), 
        -3+1j: (0, 0, 1, 1), 
        -1-3j: (0, 1, 0, 0), 
        -1-1j: (0, 1, 0, 1), 
        -1+3j: (0, 1, 1, 0), 
        -1+1j: (0, 1, 1, 1), 
        3-3j: (1, 0, 0, 0), 
        3-1j: (1, 0, 0, 1), 
        3+3j: (1, 0, 1, 0), 
        3+1j: (1, 0, 1, 1), 
        1-3j: (1, 1, 0, 0), 
        1-1j: (1, 1, 0, 1), 
        1+3j: (1, 1, 1, 0), 
        1+1j: (1, 1, 1, 1)
        }
    return list(demapping[i] for i in complex_value_array)

def matrix_of_len4_arrays(seq):
    matrix_len4_arrays = []
    for i in range(len(seq)):
        if(i % 4 == 0):
            array_len4 = []
            matrix_len4_arrays.append(array_len4)
        array_len4.append(seq[i])
    return matrix_len4_arrays

seq = [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1 ,1, 0, 1, 0, 0, 1, 0, 1]
