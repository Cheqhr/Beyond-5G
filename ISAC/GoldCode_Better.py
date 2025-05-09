#for lfsr length 5
import time as ti
import random as rnd
import numpy as np

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

def taps5and2_lfsr(clock, lfsr, binary_len):
    if(clock):
        bit = (lfsr ^ (lfsr >> 3)) & 1 #implements taps 6 and 5
        lfsr = bitShiftRight(lfsr, bit, binary_len)
    return lfsr

def taps5and4and3and2_lfsr(clock, lfsr, binary_len):
    if(clock):
        bit = (lfsr ^ (lfsr >> 1) ^ (lfsr >> 2) ^ (lfsr >> 3)) & 1 #implements taps 6 and 5
        lfsr = bitShiftRight(lfsr, bit, binary_len)
    return lfsr

def generate_gold_codes(binary_len, start_state1, start_state2):
    lfsr1 = start_state1
    lfsr2 = start_state2
    shift = False
    gold_code = []
    m_seq1 = []
    m_seq2 = []
    while True:
        shift = clock(shift)
        lfsr1 = taps5and2_lfsr(shift, lfsr1, binary_len)
        lfsr2 = taps5and4and3and2_lfsr(shift, lfsr2, binary_len)
        if(shift):
            first_bit = bit_wanted(lfsr1, binary_len - 1, binary_len)
            second_bit = bit_wanted(lfsr2, binary_len - 1, binary_len)
            m_seq1.append(first_bit)
            m_seq2.append(second_bit)
            gold_code.append(first_bit ^ second_bit)
        if(lfsr1 == start_state1 or lfsr2 == start_state2):
            break
    return gold_code, m_seq1, m_seq2

def gold_code_lfsr_seed_time_shifted(lfsr_seed, binary_len):
    lfsr_seeds_shifted = []
    lfsr_seeds_shifted.append(lfsr_seed)
    for i in range(binary_len - 1):
        bit = bit_wanted(lfsr_seed, 0, binary_len)
        lfsr_code_shifted = bitShiftRight(lfsr_seed, bit, binary_len)
        lfsr_seed = lfsr_code_shifted
        lfsr_seeds_shifted.append(lfsr_seed)
    return lfsr_seeds_shifted

def cross_correlation(seq1, seq2, shift):
    cross_correlation = 0
    for i in range(len(seq1)):
        cross_correlation += (-1) ** (seq1[(i + shift) % len(seq1)] + seq2[i])
    return cross_correlation

seed_lfsr1 = 0b01001
seed_lfsr2 = 0b10010
seeds_lfsr1 = gold_code_lfsr_seed_time_shifted(seed_lfsr1, 5)
seeds_lfsr2 = gold_code_lfsr_seed_time_shifted(seed_lfsr2, 5)
gold_codes = []
m_seq1s = []
m_seq2s = []
for i in range(len(seeds_lfsr1)):
    gold_code, m_seq1, m_seq2 = generate_gold_codes(5, seeds_lfsr1[i], seeds_lfsr2[i])
    gold_codes.append(gold_code)
    m_seq1s.append(m_seq1)
    m_seq2s.append(m_seq2)
    
cross_correlations_time_shift = []
for i in range(len(gold_codes[0]) + 1):
    cross_correlations_time_shift.append(cross_correlation(gold_codes[0], gold_codes[0], i))

print(cross_correlations_time_shift)