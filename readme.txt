Eulerian Path Genome Assembler
===============================

A Python project that assembles a reference genome from a given spectrum. Assuming that the spectrum does not contain sequencing errors, a de bruijn graph can be built from the spectrum and traversed using a singular Eulerian path. 

The input would be a spectrum and the output a is a predicted genome sequence. The output contains the headers of the spectrum in the sorted order that they appear in the genome. This is done by aligning the spectrum back to the predicted assmbled sequence.

Deliverables:
-------------

euleriangenomeassembly.py -- code for genome assembly 

predictions.txt -- output of the headers of the spectrum in the sorted order that they appear in the genome

predictions.zip -- zipped csv of predictions.txt


Usage
-----
The program takes in one input, the spectrum fasta without the genome positions 

To run the program, navigate to the project directory and run:

> python3 euleriangenomeassembly.py spectrum.fasta

The program takes the following argument:

* `--spectrum.fasta`: A fasta file of the spectrum without genome positions 

Examples
--------

Here is an example of how to run the program: 

>python3 euleriangenomeassembly.py project3a_10000_spectrum.fasta

Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For assembling a reference genome of 10000 nucleotides and spectrum alignment is for a given spectrum of ~9980 20kmers is:

real	0m0.628s
user	0m0.510s
sys	0m0.067s
