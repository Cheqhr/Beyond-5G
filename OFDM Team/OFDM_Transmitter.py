import FFT_IFFT_DFT_IDFT as fourier
import random as rnd

def weighted_bernoulli(p):
    return 1 if rnd.random() < p else 0

def serial_to_parallel(seq):
    matrix_len4_arrays = []
    for i in range(len(seq)):
        if(i % 4 == 0):
            array_len4 = []
            matrix_len4_arrays.append(array_len4)
        array_len4.append(seq[i])
    return matrix_len4_arrays

def combing_two_arrays(arr1, arr2): #adds contents of array1 to the front of array2
    len1 = len(arr1)
    len2 = len(arr2)
    length = len1 + len2
    new_arr = []
    for i in range(length):
        if(i > len1 - 1):
            new_arr.append(arr2[i - len1])
        else:
            new_arr.append(arr1[i])
    return new_arr


#16QAM Modulator
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

def cyclic_prefix_addition(seq, cyclic_prefix_len):
    cp = seq[-cyclic_prefix_len:]
    return combing_two_arrays(cp, seq)



seq = []
n = 5 #Number of QAM Symbols to be Made
p = 0.5 #Probability of a 1 from the weighted bernoulli
for i in range(4*n):
    seq.append(weighted_bernoulli(p))
cyclic_prefix_len = len(seq) // 10 #This is an assumption
matrix = serial_to_parallel(seq)
QAM_Modulated_Seq = mapping(matrix)
time_domain_QAM_seq = fourier.IDFT(QAM_Modulated_Seq)
transmit_signal = cyclic_prefix_addition(time_domain_QAM_seq, cyclic_prefix_len)
