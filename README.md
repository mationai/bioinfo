# Cogniterra/Stepik Bioinformatics/Genomics with Python

Solutions (spoilers!) for same course in 3 platforms:
- [Cogniterra](https://cogniterra.org/) (was Stepik) Genome Sequencing
- [Coursera Bioinformatics](https://www.coursera.org/specializations/bioinformatics)


## Functional Packages and Usage

Functional Programming and String Algorithms is like Peanut Butter meets Bread, so many are used here.

### pipe lib usage
```py
from pipe import permutations

for n in 'AC' | permutations(2):
   print(n)
```


## Setup / Python version

Use python 3.12

## Running

```sh
python src/cogniterra/1-1.seq     # run all exercises in 1-1.seq
python src/cogniterra/1-1.seq 1   # run exercise 1 in 1-1.seq
python src/cogniterra/1-1.seq 2:4 # run exercise 2-4 in 1-1.seq
```