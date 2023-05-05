import json
import pandas as pd

def analyze_results(filename: str):
    
    df = pd.DataFrame(columns=['benchmark', 'score'])
    
    with open(filename) as f:
        data = json.load(f)[0] # the first entry is the results (only one entry)
        
        for entry in data['entries']:
            bmark = entry['benchmark']
            score = entry['#gates']
            
            df.loc[len(df)] = [bmark, score]
    
    df.to_csv('results.csv', index=False)
            
def print_results():
    
    df = pd.read_csv('results.csv')
    
    for i, row in df.iterrows():
        print(f"{row['benchmark']}: {row['score']:5d}\t", end='')
        if i % 10 == 9:
            print()

analyze_results("experiments/contest_xag.json")
print_results()