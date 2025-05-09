import FFT_IFFT_DFT_IDFT as fourier
import random as rnd
import numpy as np
import matplotlib.pyplot as plt

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

def combing_two_arrays(arr1, arr2):
    return arr1 + arr2


#16QAM Modulator
def mapping(matrix_len4_arrays): #Upgrade to 64QAM, takes in a matrix with arrays all of length 4
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

def cyclic_prefix_removal(seq, cyclic_prefix_len):
    return seq[cyclic_prefix_len:cyclic_prefix_len + len(seq)]

def rounding_list(seq):
    for i in range(len(seq)):
        real_num = seq[i].real
        complex_num = seq[i].imag
        seq[i] = round(real_num) + round(complex_num) * 1j
    return seq

def plotter_fn(seq, seq2, scatter_plot): #decomposes 2 arrays into their complex and real parts and plots them
#This code is unnecessary for the transmitter, but is kept in since it is a useful visualization and debugging tool
    complex= []
    real = []
    temp = []
    complex1 = []
    real1 = []
    temp1 = []
    if(scatter_plot == 'scatter'):
        for i in range(len(seq)):
            complex.append(seq[i].imag)
            real.append(seq[i].real)
        for i in range(len(seq2)):
            complex1.append(seq2[i].imag)
            real1.append(seq2[i].real)
        for i in range(len(real)):
            temp.append(i)
        for i in range(len(real1)):
            temp1.append(i)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True)
        f.supylabel("Amplitude (V)")
        f.supxlabel("n")
        ax1.scatter(temp, complex)
        ax1.set_title('Complex Orginal')
        ax1.grid()
        ax2.scatter(temp, real)
        ax2.set_title('Real Orginal')
        ax2.grid()
        ax3.scatter(temp1, complex1)
        ax3.set_title('Complex Upsampled')
        ax3.grid()
        ax4.scatter(temp1, real1)
        ax4.set_title('Real Upsampled')
        ax4.grid()
        #plt.show()
    elif(scatter_plot == "plot"):
        for i in range(len(seq)):
            complex.append(seq[i].imag)
            real.append(seq[i].real)
        for i in range(len(seq2)):
            complex1.append(seq2[i].imag)
            real1.append(seq2[i].real)
        for i in range(len(real)):
            temp.append(i)
        for i in range(len(real1)):
            temp1.append(i)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True)
        f.supylabel("Amplitude (V)")
        f.supxlabel("n")
        ax1.plot(temp, complex)
        ax1.set_title('Pulse Shaped Complex Component')
        ax1.grid()
        ax2.plot(temp, real)
        ax2.set_title('Pulse Shaped Real Component')
        ax2.grid()
        ax3.plot(temp1, complex1)
        ax3.set_title('Complex Upsampled')
        ax3.grid()
        ax4.plot(temp1, real1)
        ax4.set_title('Real Upsampled')
        ax4.grid()
        #plt.show()
    else:
        for i in range(len(seq)):
            complex.append(seq[i].imag)
            real.append(seq[i].real)
        for i in range(len(seq2)):
            complex1.append(seq2[i].imag)
            real1.append(seq2[i].real)
        for i in range(len(real)):
            temp.append(i)
        for i in range(len(real1)):
            temp1.append(i)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True)
        f.supylabel("Amplitude (V)")
        f.supxlabel("n")
        ax1.plot(temp, complex)
        ax1.set_title('Pulse Shaped Complex Component')
        ax1.grid()
        ax2.plot(temp, real)
        ax2.set_title('Pulse Shaped Real Component')
        ax2.grid()
        ax3.scatter(temp1, complex1)
        ax3.set_title('Complex Upsampled')
        ax3.grid()
        ax4.scatter(temp1, real1)
        ax4.set_title('Real Upsampled')
        ax4.grid()
        #plt.show()

def upsampling(seq, interp): 
    upsampled_seq = []
    zeros = interp * [0]
    for i in range(len(seq)):
        upsampled_seq.append(seq[i])
        upsampled_seq = upsampled_seq + zeros
    print(upsampled_seq)
    return upsampled_seq

def raised_root_cosine_filter(n, T_s, beta):
    if(n == 0):
        return (1 / T_s) * (1 + beta * (4 / np.pi - 1))
    if(n == (T_s / (4 * beta))):
        return (beta / (T_s * (2 ** 0.5))) * ((1 + 2 / np.pi) * np.sin(np.pi / (4 * beta)) + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta)))
    numerator = np.sin((np.pi * n * (1 - beta)) / T_s) + ((4 * beta * n) / T_s) * np.cos((np.pi * n * (1 + beta)) / T_s)
    denominator = ((np.pi * n) / T_s) * (1 - ((4 * beta * n) / T_s) ** 2)
    return (numerator / denominator) * 1 / T_s

def pulsed_shaped_filter(seq, T_s, beta):
    x = np.linspace(-50, 50, 100)
    raised_cosine_output = []
    for i in x:
        raised_cosine_output.append(raised_root_cosine_filter(i, T_s, beta))
    return 3*np.convolve(seq, raised_cosine_output)

def modulator(seq, carrier_freq):
    new_seq = []
    time = np.linspace(0, 100, len(seq)) #ask on this
    for i in range(len(seq)):
        new_seq.append(seq[i] * np.cos(2 * np.pi * carrier_freq * time[i]))
    return new_seq

seq = []
n = 4 #Number of QAM Symbols to be Made
p = 0.5 #Probability of a 1 from the weighted bernoulli
upsample_factor = 20 #Interpolation value for transmitted OFDM symbol
T_s = 18 #sampling period
beta = 1 #roll off constant, get as close to zero without causing issues
for i in range(4*n):
    seq.append(weighted_bernoulli(p))
cyclic_prefix_len = len(seq) // 10 #This is an assumption, needs to be longer than impulse response of the channel
binary_symbol_matrix = serial_to_parallel(seq)
QAM_Modulated_Seq = mapping(binary_symbol_matrix)
time_domain_QAM_seq = fourier.IDFT(QAM_Modulated_Seq)
OFDM_Symbol = cyclic_prefix_addition(time_domain_QAM_seq, cyclic_prefix_len)
print(OFDM_Symbol)
OFDM_Symbol_Upsampled = upsampling(OFDM_Symbol, upsample_factor)
pulsed_signal = pulsed_shaped_filter(OFDM_Symbol_Upsampled, T_s, beta)
transmitted_signal = modulator(pulsed_signal, 1000)
plotter_fn(OFDM_Symbol, OFDM_Symbol_Upsampled, "scatter")
plotter_fn(pulsed_signal, OFDM_Symbol_Upsampled, "")
plotter_fn(pulsed_signal, transmitted_signal, "plot")
#Channel needs to be built

#Receiver Architechture
received_signal = cyclic_prefix_removal(OFDM_Symbol, cyclic_prefix_len) #Make sure receiver can take in upsampled and pulsed OFDM sequence 
fourier_transform = fourier.DFT(received_signal)
fourier_transform_rounded = rounding_list(fourier_transform) #Needs to be removed later
QAM_Demod_Seq = demapping(fourier_transform_rounded)
plt.show()
