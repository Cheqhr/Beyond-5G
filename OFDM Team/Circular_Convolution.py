import numpy as np

x_vec = [1, 2, 2, 1]
h_vec = [1, 2, 3, 1]

def cycle(index, array):
    return array[-index:] + array[:-index]

def circular_convolution(x_vec, h_vec):
    N = len(x_vec)
    array = np.zeros((N, N))
    for i in range(N):
        if(i == 0):
            array[:, 0] = np.transpose(x_vec)
        col_i = cycle(i, x_vec)
        array[:, i] = np.transpose(col_i)
    
    circular_convolution_res = np.matmul(array, h_vec)

    return circular_convolution_res