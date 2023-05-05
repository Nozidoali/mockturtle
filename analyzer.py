import json
import pandas as pd


def analyze_results(filename: str):

    results = {}

    with open(filename) as f:
        data = json.load(f)[0] # the first entry is the results (only one entry)
        
        for entry in data['entries']:
            bmark = entry['benchmark']
            score = entry['#gates']
            
            results[bmark] = score

    return results        
            
def print_results(results):
    row = 0
    for i in results:
        score = results[i]
        print(f"{i}: {score:5d}\t", end='')
        row += 1
        if row % 10 == 0:
            print()

