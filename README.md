Motif_mark.py is a script that, when called, creates a graphical output depicting where exons 
in genes are found, and where motifs are found. It outputs a single graphic in .svg form, 
with all genes analyzed in the same file.

#######################################

Motif_mark.py is called from the command line, in bash. Before this is done, a set-up needs 
to be completed. This is done by creating an environment.

#######################################

Creating Environment:

In order to properly run the script, pycairo is needed. To install and use the 
pycairo module, an environment must be created. This can be done in the following way, 
in this example, creating an environment named "pycairo_env".

Only the first 3 lines are needed to make the environment - the last line is used to update 
conda if it is not current, and may not be necessarily needed.

	conda create --name pycairo_env
	conda activate pycairo_env/
	conda install -c conda-forge pycairo
	conda update -n base defaults conda

After the environment is created, it must be activated, in the following way: 

	conda activate [/dir/..]/pycairo_env

It is best practice to activate the environment in the same directory that motif_marker.py 
is being run from.

#######################################

Running motif_mark.py:

Motif_mark.py is run from the command line:

./motif_mark.py -fa [fasta file input] -m [motif list file input] [--help]

The first two arguments, -fa and -m, are required.

The first argument, -fa, is for the input of the FASTA file for processing. This should contain 
a list of genes with headers in FASTA format, with the corresponding gene sequences below them. 
The gene sequences in the FASTA file should hold the convention that nucleotides in exons are 
capitalized, and introns are lowercase.

The second argument, -m, is for the input of the motif list file. This file should contain a list 
of the motifs that are being searched for in the genes in the FASTA file. The motifs can contain 
ambiguous IUPAC symbols, and indeed, this program is written with that in mind. Currently, this 
program only has the capacity for 10 motifs. Each motif in the input should be on a new line.

The --help argument is optional, and will display the following message:

'-fa' to specify fasta file; '-m' to specify motif mark file

-h, --help            show this help message and exit
  -fa FASTA_FILE, --fasta_file FASTA_FILE
                        input file with fasta sequences
  -m MOTIF_MARK_FILE, --motif_mark_file MOTIF_MARK_FILE
                        input file with motif marks


#######################################

Output: 

The output of the file will be all the genes outputs from the FASTA file as one combined 
.svg file. Each gene's output will be of the same format:

The gene symbol from the FASTA header line at the top.

A graphic, to scale, of the length of the gene in light blue, with a narrow width for introns, 
and a larger width for exons corresponding to the positions of the two (in nucleotide base 
pairs) in the gene sequence given.

Each motif (up to 10 total) found in the gene is marked in a different color by a 
narrow, long vertical bar at the start of the position at which it is found in the gene, and 
a shorter, wider vertical bar of the same color that spans the length of the motif in the 
gene. The longer bar marks the start of the motif in order to provide more clarity when motifs 
overlap - it is not exactly 1 nucleotide long, so it is not to scale - it is the smallest 
width that the graphical output of pycairo will allow for.

Under the graphic of the gene, the box labeled "Motif Legend" contains each motif that is 
marked on the current gene's graphic - both the color that it is marked in, and the 
corresponding motif sequence (with ambiguous IUPAC symbols substituted for nucleotides where 
appropriate).
