PconsFold
=========

A pipeline for protein folding using predicted contacts from PconsC and a  Rosetta folding protocol.

![PconsFold pipeline](https://github.com/ElofssonLab/pcons-fold/blob/release-1.0/pipeline_horiz.png)


Instructions:
-------------

Pipeline overview:

1. Input: fasta file containing one sequence
2. Prepare input for PconsC
3. Contact prediction with PconsC
4. Prepare input for Rosetta folding
5. Rosetta folding
6. Extract and relax structures with lowest Rosetta energy
7. Output: the predicted contact map (also as a plot) and the lowest energy structure.


Dependencies:

- Rosetta v3.5 or weekly
- Jackhmmer from HMMER v3.0 or higher
- HHblits from HHsuite v2.0.16 or higher
- PSICOV v1.11 or higher
- either MATLAB v8.1 or higher
- or MATLAB Compiler Runtime (MCR) v8.1 or higher 

MATLAB is needed to run plmDCA. However, if MATLAB is not available you can also use a compiled version of plmDCA. For the compiled version to run you need to provide a path to MCR.


How to run it:

Make sure all paths are correct in localconfig.py. Remark: except for plmDCA, dependencies are not included in this repository.

```
./run_pipeline.py  [-c n_cores] [-n n_decoys] [--norelax] [--nohoms] hhblits_database jackhmmer_database sequence_file
```

```
python pconsc/predictAll_1.0.py [-c cores] hhblits_database jackhmmer_database sequence_file
```


Where "pconsc_output" is the contact map as given by PconsC "(sequence_file).pconsc.out"


``` 
python folding/rosetta/prepare_input.py sequence_file pconsc_output factor

python folding/rosetta/fold.py [-c cores] sequence_file rosetta_constraints

python folding/rosetta/extract.py [-c cores] number_of_extracted_structures relax_flag
```

Where "factor" is a float denoting the fraction (with respect to the length of the sequence) of top ranked contacts to use during folding.
This script generates the file `(pconsc_output)-(factor).constraints` which is then used by Rosetta in the next step. 
"rosetta_constraints" is the output from `prepare_input.py` and "relax_flag" is either 1 or 0 deciding whether the output structures should be relaxed or not respectively. We recommend to set it to 1 for better stereochemistry.
