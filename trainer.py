#!/usr/bin/env python
# -*- encoding=utf8 -*-

'''
Author: Hanyu Wang
Created time: 2023-05-05 16:39:37
Last Modified by: Hanyu Wang
Last Modified time: 2023-05-05 19:11:40
'''

import os
import subprocess
from generator import generate_xag
from analyzer import analyze_results, print_results

def run_exp(skip_make: bool = False):
    
    solutions = {}

    for i in range(100):
        filename = f"./experiments/training_benchmarks/ex{i:02d}.truth"
        problem, solution = generate_xag(6, 1000, 1)
        with open(filename, 'w') as f:
            f.write(problem)
        solutions[f"ex{i:02d}"] = solution

    while len(os.listdir('./experiments/training_benchmarks')) != 100:
        print(f"{len(os.listdir('./experiments/training_benchmarks')):02d} / 100.", end='\r')
        subprocess.run('sleep 1', shell=True)
        
    subprocess.run('rm -f ./experiments/agent_xag.json', shell=True)

    if not skip_make:
        subprocess.run('rm -r build && mkdir build && cd build && cmake -DMOCKTURTLE_EXPERIMENTS=ON -DCMAKE_BUILD_TYPE=RELEASE ..  && make agent_xag && ./experiments/agent_xag', shell=True)
    else:
        subprocess.run('cd build && ./experiments/agent_xag', shell=True)
    
    # we check the number of files in the folder is 100
    while len(os.listdir('./experiments/training_results')) != 100:
        print(f"{len(os.listdir('./experiments/training_results')):02d} / 100.", end='\r')
        subprocess.run('sleep 1', shell=True)

    results = analyze_results('./experiments/agent_xag.json')

    useful_results = {}

    # get the difference
    for i in results:
        results[i] = solutions[i] - results[i]
        if results[i] < 0:
            print(f"INFO: {i} is a counter example.")
            useful_results[i] = results[i]

    print_results(results)

    # remove the benhmarks
    for i in range(100):
        filename = f"./experiments/training_benchmarks/ex{i:02d}.truth"
        subprocess.run(f"rm -f {filename}", shell=True)
    
    # remove the results
    for i in range(100):
        filename = f"./experiments/training_results/ex{i:02d}.aig"
        subprocess.run(f"rm -f {filename}", shell=True)
    
    # remove the json
    subprocess.run('rm -f ./experiments/agent_xag.json', shell=True)

if __name__ == "__main__":
    run_exp(skip_make=True)
