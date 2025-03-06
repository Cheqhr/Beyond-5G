#2 LFSRs with 2**N - 1 states that shift right every clock cycle
#I will use Fibonacci Implementation
#taps 6, 5
import time as ti

def clock(shift): #1Hz Clock, Hz = 1/(2 * delay_time)
    delay_time = 0
    if(shift != True):
        ti.sleep(delay_time)
        return True
    else:
        ti.sleep(delay_time)
        return False

def bit_wanted(n, bit_num, binary_len): #only works from 0 to binary_len - 1
    if(binary_len < bit_num | bit_num <= 0):
        return
    return (n >> bit_num) & 1

def bitShiftRight(lfsr, bit, binary_len):
    return ((lfsr >> 1) | (bit << binary_len - 1))

def taps6and5_lfsr(clock, lfsr, binary_len):
    if(clock):
        bit = (lfsr ^ (lfsr >> 1)) & 1 #implements taps 6 and 5
        lfsr = bitShiftRight(lfsr, bit, binary_len)
    return lfsr

def generate_gold_codes(binary_len):
    start_state1 = 1 << binary_len - 1 | 1 #make sure start sequences have low cross correlation
    start_state2 = start_state1 | 1 << 3
    lfsr1 = start_state1 #could seed the lfsr's to generate different sequences each time, problem is cross correlation 
    lfsr2 = start_state2
    shift = False
    gold_codes = []
    while True:
        shift = clock(shift)
        lfsr1 = taps6and5_lfsr(shift, lfsr1, binary_len)
        lfsr2 = taps6and5_lfsr(shift, lfsr2, binary_len)
        if(shift):
            gold_codes.append(lfsr1 ^ lfsr2)
        if(lfsr1 == start_state1 or lfsr2 == start_state2):
            break
    return gold_codes

def gold_code_time_shifted(gold_code, binary_len):
    gold_codes_shifted = []
    for i in range(binary_len - 1):
        bit = bit_wanted(gold_code, 0, binary_len)
        gold_code_shifted = bitShiftRight(gold_code, bit, binary_len)
        gold_code = gold_code_shifted
        gold_codes_shifted.append(gold_code)
    return gold_codes_shifted


def gold_codes_autocorrelation(binary_len): #binary_len is the length of your binary number. If you have '0001' the binary_len is 4
    gold_codes = generate_gold_codes(binary_len)
    auto_correlation_matrix = []
    for i in range(len(gold_codes)):
        auto_correlation_per_gold_code = []
        gold_codes_time_shift = gold_code_time_shifted(gold_codes[i], binary_len)
        for j in range(len(gold_codes_time_shift)):
            agreements = 0
            for k in range(binary_len):
                if(bit_wanted(gold_codes[i], k, binary_len) == bit_wanted(gold_codes_time_shift[j], k, binary_len)):
                    agreements += 1
            disagreements = binary_len - agreements
            auto_correlation_per_gold_code.append(agreements - disagreements)
        auto_correlation_matrix.append(auto_correlation_per_gold_code)
    return auto_correlation_matrix

def gold_codes_crosscorrelation(binary_len):
    gold_codes = generate_gold_codes(binary_len)
    cross_correlation_matrix = []
    for i in range(len(gold_codes)):
        cross_correlation_per_gold_code = []
        for j in range(len(gold_codes)):
            agreements = 0
            for k in range(binary_len):
                if(bit_wanted(gold_codes[i], k, binary_len) == bit_wanted(gold_codes[j], k, binary_len)):
                    agreements += 1
            disagreements = binary_len - agreements
            if(i != j):
                cross_correlation_per_gold_code.append(agreements - disagreements)
        cross_correlation_matrix.append(cross_correlation_per_gold_code)
    return cross_correlation_matrix