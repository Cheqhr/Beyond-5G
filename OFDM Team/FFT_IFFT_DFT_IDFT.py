import cmath as cm
import numpy as np
 
def DFT(seq): #Works for all lists
    N = len(seq)
    fourier_coefficient = 0
    fourier_coefficients = []
    for i in range(N):
        for j in range(len(seq)):
            exponential = cm.exp(((cm.pi*2) / N) * i * j * -1j)
            fourier_coefficient += exponential * seq[j]
        fourier_coefficients.append(fourier_coefficient)
        fourier_coefficient = 0
    return fourier_coefficients

def IDFT(seq):
    N = len(seq)
    fourier_coefficient = 0
    fourier_coefficients = []
    for i in range(N):
        for j in range(len(seq)):
            exponential = cm.exp(((cm.pi*2) / N) * i * j * 1j)
            fourier_coefficient += exponential * seq[j]
        fourier_coefficients.append(fourier_coefficient / N)
        fourier_coefficient = 0
    return fourier_coefficients

def FFT(seq): #Only works for lists when len(N) % 2 == 0. 
    N = len(seq)
    if N == 1:
        return [seq[0]]

    w_n = list(complex(np.cos(((-2*cm.pi) / N) * i), np.sin(((-2*cm.pi) / N) * i)) for i in range(N))   
     
    even_elem = seq[0::2]
    odd_elem = seq[1::2]
        
    y_even = FFT(even_elem)
    y_odd = FFT(odd_elem)
    
    y = [0]*N
    
    for i in range(N//2):
        y[i] =  y_even[i] + w_n[i] * y_odd[i]
        y[i + N//2] = y_even[i] - w_n[i] * y_odd[i]
    
    return y

def IFFT(seq):
    N = len(seq)
    if N == 1:
        return [seq[0]]

    w_n = list(complex(np.cos(((2*cm.pi) / N) * i), np.sin(((2*cm.pi) / N) * i)) for i in range(N))   
     
    even_elem = seq[0::2]
    odd_elem = seq[1::2]
        
    y_even = FFT(even_elem)
    y_odd = FFT(odd_elem)
    
    y = [0]*N
    
    for i in range(N//2):
        y[i] =  (y_even[i] + w_n[i] * y_odd[i]) / N
        y[i + N//2] = (y_even[i] - w_n[i] * y_odd[i]) / N
    
    return y

def fourier_transform_cleanup(seq):
    for i in range(len(seq)):
        real_num = seq[i].real
        if(abs(real_num) < 1e-13):
            real_num = 0
            
        complex_num = seq[i].imag
        if(abs(complex_num) < 1e-13):
            complex_num = 0
        
        seq[i] = real_num + complex_num*1j
        
    return seq