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


def print_large_ids(results):
    ids = []
    for i in results:
        if results[i] > 1000:

            id = int(i.replace("ex", ""))
            ids.append(str(id))
    
    print("{" + ",".join(ids) + "}")

if __name__ == "__main__":
    results = analyze_results("./experiments/contest_xag.json")
    # print_results(results)
    print_large_ids(results)