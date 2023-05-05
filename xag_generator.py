import numpy as np
import random

class xag:
    def __init__(self) -> None:
        pass
    

def tt_to_string(tt):
    
    tt_string = ''
    for i in range(len(tt)):
        tt_string += '1' if tt[i] else '0'
    
    return tt_string

def generate_xag(num_inputs, num_gates, num_outputs):
    
    tts = []
    tt_strings = set()
    referenced_outputs = set()
    
    # prepare the base truth table
    tts.append(np.zeros(2**num_inputs, dtype=bool))
    for idx in range(num_inputs):        
        tt = np.zeros(2**num_inputs, dtype=bool)
        for i in range(2**num_inputs):
            tt[i] = (i>>idx) & 0b1;
        tts.append(tt)
        tt_strings.add(tt_to_string(tt))
    
    for _ in range(num_gates):
        tt1, tt2 = random.choices(tts, k=2)
        
        referenced_outputs.add(tt_to_string(tt1))
        referenced_outputs.add(tt_to_string(tt2))
        
        op = random.choice(['and', 'xor'])
        
        polar1 = random.choice([True, False])
        polar2 = random.choice([True, False])
        
        tt1 = np.logical_xor(tt1, polar1)
        tt2 = np.logical_xor(tt2, polar2)
        
        if op == 'and':
            tt = np.logical_and(tt1, tt2)
        elif op == 'xor':
            tt = np.logical_xor(tt1, tt2)
        
        if tt_to_string(tt) in tt_strings:
            continue
        tts.append(tt)
        tt_strings.add(tt_to_string(tt))
        
    candidates = []
    for cand in tts:
        cand_string = tt_to_string(cand)
        if cand_string not in referenced_outputs:
            continue
        candidates.append(cand)
        
    outputs = random.choices(candidates, k=num_outputs)
    
    for tt in outputs:
        print(tt_to_string(tt) + "\n")
    
generate_xag(8, 2000, 4)