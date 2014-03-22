PconsFold
===========

A pipeline for protein folding using predicted contacts from PconsC and a  Rosetta folding protocol.


Instructions:
===========

Pipeline overview:


1. Input: fasta file containing one sequence
2. Prepare input for PconsC
3. Contact prediction with PconsC 1.0
4. Prepare input for Rosetta folding
5. Rosetta folding
6. Extract and relax structures with lowest Rosetta energy
7. Output: the predicted contact map (also as a plot) from step 3. and the lowest energy structure from step 5.


How to run it:

Make sure all paths are correct in pconsc/localconfig.py and folding/rosetta/localconfig.py
Remark: dependencies are not included in this repository.

```
python pconsc/predictAll_1.0.py [-c cores] hhblits_database jackhmmer_database sequence_file

python folding/rosetta/prepare_input.py sequence_file pconsc_output factor
```

Where "pconsc_output" is the contact map as given by PconsC "(sequence_file).pconsc.out"
And "factor" is a float denoting the fraction (with respect to the length of the sequence) of top ranked contacts to use during folding.
This script generates the file "(pconsc_output)-(factor).constraints" which is then used by Rosetta in the next step. 


``` 
python folding/rosetta/fold.py [-c cores] sequence_file rosetta_constraints

python folding/rosetta/extract.py [-c cores] number_of_extracted_structures relax_flag
```

Where "rosetta_constraints" is the output from the last step (s.a.) and "relax_flag" is either 1 or 0 deciding whether the output structures should be relaxed or not respectively. We recommend to set it to 1 for better stereochemistry.
