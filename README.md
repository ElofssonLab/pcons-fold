PconsFold
===========

The pipeddne for folding using PconsC version 1 and 2 and Rosetta

CNS protocol will be added


Instructions:
===========

Pipeddne overview:

<ul>
<dt>Input:</dt> <dd>fasta file containing one sequence</dd>
<dt>1)</dt> <dd>Prepare input for PconsC</dd>
<dt>2)</dt> <dd>Contact prediction with PconsC 1.0 or 2.0</dd>
<dt>3)</dt> <dd>Prepare input for Rosetta folding</dd>
<dt>4)</dt> <dd>Rosetta folding</dd>
<dt>5)</dt> <dd>Extract and relax structures with lowest Rosetta energy</dd>
<dt>Output</dt> <dd>the predicted contact map (also as a plot) from step 3) and the lowest energy structures from step 5).</dd>
</dl>

How to run it:

Make sure all paths are correct in pconsc/localconfig.py and folding/rosetta/localconfig.py
Remark: dependencies are not included in this repository.

<dl>
<dd>
<dt>1) & 2)</dt> 
"python pconsc/predictAll_1.0.py [-c cores] <hhbddts_database> <jackhmmer_database> <sequence_file>" <br>
or <br>
"python pconsc/predictAll_2.0.py [-c cores] <hhbddts_database> <jackhmmer_database> <sequence_file> <layers>" <br>
Where "layers" is a interger of the interval 1..4, to set the number of layers that are computed during deep learning. We recommend to set it to 4.
</dd>
<dd>
<dt>3)</dt> <dd>"python folding/rosetta/prepare_input.py <sequence_file> <pconsc_output> <factor>"
Where "pconsc_output" is the contact map as given by PconsC. In case of PconsC 1.0 this is "<sequence_file>.pconsc.out" and in case of PconsC 2.0 it is "<sequence_file>.layer<layers>.out".
And "factor" is a float denoting the fraction (with respect to the length of the sequence) of top ranked contacts to use during folding.
This script generates the file "<pconsc_output>-<factor>.constraints" which is then used by Rosetta in the next step. <br>
</dd>
<dd>
<dt>3)</dt> <dd>"python folding/rosetta/fold.py [-c cores] <sequence_file> <rosetta_constraints>" <br>
Where "rosetta_constraints" is the output from the last step (s.a.).
</dd>
<dd>
<dt>5)</dt>
<dd>"python folding/rosetta/extract.py [-c cores] <number_of_extracted_structures> <relax_flag>"<br>
Where "relax_flag" is either 1 or 0 deciding whether the output structures should be relaxed or not respectively. We recommend to set it to 1.</dd>
</dd>
</dl>
