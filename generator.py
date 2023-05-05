import numpy as np
import random

class xag:
    def __init__(self) -> None:
        self.tts = []
        self.tt_string_to_node = {}
        self.referenced_outputs = set()
    
    def insert_tt(self, tt, node_tuple=None):
        tt_string = tt_to_string(tt)
        if tt_string in self.tt_string_to_node:
            return False
        self.tts.append(tt)
        self.tt_string_to_node[tt_string] = (tt, node_tuple)
        return True

    def create_inputs(self, num_inputs):
        # prepare the base truth table
        self.insert_tt(np.zeros(2**num_inputs, dtype=bool))
        for idx in range(num_inputs):        
            tt = np.zeros(2**num_inputs, dtype=bool)
            for i in range(2**num_inputs):
                tt[i] = (i>>idx) & 0b1;
            self.insert_tt(tt)

    def creat_random_node(self):
        tt1, tt2 = random.choices(self.tts, k=2)
        
        self.referenced_outputs.add(tt_to_string(tt1))
        self.referenced_outputs.add(tt_to_string(tt2))
        
        op = random.choice(['and', 'xor'])
        
        polar1 = random.choice([True, False])
        polar2 = random.choice([True, False])
        
        _tt1 = np.logical_xor(tt1, polar1)
        _tt2 = np.logical_xor(tt2, polar2)
        
        if op == 'and':
            tt = np.logical_and(_tt1, _tt2)
        elif op == 'xor':
            tt = np.logical_xor(_tt1, _tt2)
        
        tts1 = tt_to_string(tt1)
        tts2 = tt_to_string(tt2)

        return self.insert_tt(
            tt,
            (tts1, tts2, op, polar1, polar2))
    
    def create_random_xag(self, num_gates, num_outputs):
        for _ in range(num_gates):
            while not self.creat_random_node():
                pass
    
    def create_problem(self, num_outputs):
        candidates = []
        for cand in self.tts:
            cand_string = tt_to_string(cand)
            if cand_string not in self.referenced_outputs:
                continue
            candidates.append(cand)
            
        outputs = random.choices(candidates, k=num_outputs)
        
        problem = ''
        for tt in outputs:
            problem += tt_to_string(tt) + '\n'

        solution = 0
        visited = set()

        def dfs(tt_string: str):
            nonlocal solution
            nonlocal visited
            assert tt_string in self.tt_string_to_node
            _, node = self.tt_string_to_node[tt_string]
            if node is None:
                return
            left, right, op, _, _ = node
            if tt_string in visited:
                return
            visited.add(tt_string)
            dfs(left)
            dfs(right)
            if op == 'and':
                solution += 1
            elif op == 'xor':
                solution += 1
            else:
                raise Exception('Unknown operator')
            
        for tt in outputs:
            dfs(tt_to_string(tt))
        
        return problem, solution

def tt_to_string(tt):
    
    tt_string = ''
    for i in range(len(tt)):
        tt_string += '1' if tt[i] else '0'
    
    return tt_string

def generate_xag(num_inputs, num_gates, num_outputs):
    
    xag_obj = xag()
    xag_obj.create_inputs(num_inputs)
    xag_obj.create_random_xag(num_gates, num_outputs)
    return xag_obj.create_problem(num_outputs)

if __name__ == "__main__": 
  generate_xag(8, 2000, 4)