# Evaluating four different seeding algorithms:
1. All overlapping seeds,
2. Minimizer seeds,
3. Spaced seeds,
4. Genome-on-Diet seeds.

## How to use it?
python kc-py1.py [SequenceFile] [Kmer Size] [Minimizer Window Size] [1] [ComparisonMode (0: 1-to-1, 1: 1-to-many] [Pattern Sequence] [Debugging Mode?]

## Examples
Using a pattern sequence of '110':
```
python kc-py1.py 1000.5 8 6 1 0 110 0 > out2.csv
```
Using a pattern sequence of '10':
```
python kc-py1.py 1000.5 8 6 1 0 10 0 > out2.csv
```
Using a pattern sequence of '101001':
```
python kc-py1.py 1000.5 8 6 1 0 101001 0 > out2.csv
```
Using a pattern sequence of '100':
```
python kc-py1.py 1000.5 8 6 1 0 100 0 > out2.csv
```

## Visualization
You can copy the CSV results you got from our program and paste them into the first five columns in any sheets of EvalResults-June2023.xlsx and you will get similar plots to what we have in the Genome-on-Diet paper (Figure 4).

## Acknowledgment
- We build this tool on top of https://github.com/lh3/kmer-cnt/blob/master/kc-py1.py, hence the name of the file.
- We use the following sequence pairs: https://github.com/lorrainea/Seedability/blob/main/data/synthetic/1000.5, however our tool uses only the first sequence from each pair and builds its own mutated copy (with a mutation rate randomly chosen between 0% (no mutation) to 50% (highly mutated) of the sequence length for each input sequence pair).


## Feedback
We welcome any input or feedback on this evaluation.
