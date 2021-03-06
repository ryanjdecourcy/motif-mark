#!/usr/bin/env python

# Importing modules
# in bash terminal run "pip install pycairo"
import cairo
import regex
import argparse

# Setting up argparse inputs
parser = argparse.ArgumentParser(description="'-fa' to specify fasta file; '-m' to specify motif mark file")
parser.add_argument("-fa", "--fasta_file", help="input file with fasta sequences", required = True)
parser.add_argument("-m", "--motif_mark_file", help="input file with motif marks", required = True)
args = parser.parse_args()

fasta_input = args.fasta_file
motif_input = args.motif_mark_file

###############################################################
# Creating classes
###################################################################

# FASTA header class - did NOT include font (yet)
class FastaHeader:
    def __init__(self, the_header):
        '''FASTA Header class - contains information about the FASTA header:
        start - an integer defining the location on the output of the start of the header's text
        text - a string containing what the header should be'''

        # Data
        self.start = []
        self.text = the_header

#aheader = FastaHeader('121', 'FASTA_header')

#################################


# Gene class - still don't know if "the_stat" should be an int - and whether
# a location/position, or a nucleotide number in the actual genome
class Gene:
    def __init__(self, the_length):
        '''Gene - contains information about the actual gene to be drawn:
        start - the position on the output that the gene should start being drawn at
        length - how far the gene is drawn from the start to its end
        width - how wide the gene should be when it is drawn
        [same dimensional frame of reference as exon and motif objects]'''

        # Data
        self.start = []
        self.length = int(the_length)

#agene = Gene("the start of the gene", "121", "10")

#################################

# Exon class

class Exon:
    def __init__(self, the_length):
        '''Contains the information needed to draw the exon object:
        start - the position on the output that the exon should start being drawn at
        length - the length that the exon should extend to, along the gene, from its start
        width - how wide the exon should be when it's drawn
        [same dimensional frame of reference as gene and motif objects]
        '''

        # Data
        self.start = []
        self.length = int(the_length)

# anexon = Exon(100, 500, 10)

#################################

# Motif class - this is going to need to be made into a plurality

class Motif:
    def __init__(self, the_length, the_sequence):
        '''Contains the information to draw a single motif object:
        (it is possible there will be multiple motif objects per gene):
        start - the position on the output that the motif should start being drawn at
        length - how far from the start position the motif should extend
        width - how wide the motif should be
        [same dimensional frame of reference as gene and exon objects]'''

        # Data
        self.start = []
        self.length = int(the_length)
        self.sequence = the_sequence

'''amotif = Motif(848, 189, 15)
bmotif = Motif(777, 10, 10)
cmotif = Motif(888, 10, 8)'''


# GeneGroup class
class GeneGroup:
    def __init__(self, group_gene, group_exon, group_header, group_rank):
        '''All user created classes combined into one unique class.
        Added in the following order (all previously defined user classes):
        Gene, exon, header, rank (which position in sequence this gene is out of all the genes)
        & the GeneGroupObject.end_motifs; a list, which can be appended to, 
        containing a list of the motifs from the Motif user defined class'''

        # Data
        self.end_gene = group_gene
        self.end_exon = group_exon
        self.end_motifs = []
        self.end_header = group_header
        self.end_rank = int(group_rank)



####### DONE DEFINING USER CREATED CLASSES
###############################################################
################# BUILDING MOTIF CONVERSION FUNCTION BELOW

def motif_conversion(amotif):

    reg_building = ''

    for x in range(len(amotif)):
        if (len(vIUPAC[amotif[x]])) == 1:
            reg_building += (vIUPAC[amotif[x]][0])
        
        elif (len(vIUPAC[amotif[x]])) > 1:
            reg_building += "["
            for y in range(len(vIUPAC[amotif[x]])):
                reg_building += (vIUPAC[amotif[x]][y])
            reg_building += "]"

    return reg_building

########### END OF MOTIF CONVERSION FUNCTION
###############################################################
######## Start of exon_placer function

def exon_placer(aread):

    exon_regex = "[A-Z]+"

    # Just edited the below line to add overlapped = True flag
    test_exon_indices = [match.start() for match in regex.finditer(exon_regex, aread, overlapped = True)]
    test_exon = regex.findall(exon_regex, aread, overlapped = True)

    anexon_reg_found = regex.search(exon_regex, aread, overlapped = True)

    anexon = anexon_reg_found[0]

    exon_start_indices = [match.start() for match in regex.finditer(exon_regex, aread, overlapped = True)]
    exon_start = exon_start_indices[0]
    exon_end = exon_start + len(anexon)

    # list to output: the exon sequence, start position and end position 
    read_exon_start_end = [anexon, exon_start, exon_end]
    
    return read_exon_start_end


# so exon_placer fn() outputs a list of 3 elements:
# the exon found, the start position and the end position

###### End of exon_placer function
#################################################################


# Dict to translate ambiguous IUPAC nt symbols
# Credit to wikipedia: https://en.wikipedia.org/wiki/Nucleic_acid_notation)
# Although the dictionary was expanded to include lowercase nt symbols


vIUPAC = {
    "A":["A", 'a'            ],
    "C":[    "C", 'c'        ],
    "G":[        "G", 'g'    ],
    "T":[            "T", 't', 'U', 'u'],
    "U":[            "U", 'u', 'T', 't'],
    "W":["A", 'a',        "T", 't'],
    "S":[    "C", 'c', "G", 'g'    ],
    "M":["A", 'a', "C", 'c'        ],
    "K":[        "G",'g', "T", 't'],
    "R":["A", 'a',   "G", 'g'  ],
    "Y":[    "C", 'c',   "T", 't'],
    "B":[    "C", 'c', "G", 'g', "T", 't'],
    "D":["A", 'a',    "G", 'g', "T", 't'],
    "H":["A", 'a', "C", 'c',    "T", 't'],
    "V":["A", 'a', "C", 'c', "G", 'g'   ],
    "N":["A", 'a', "C", 'c', "G", 'g', "T", 't'],
    "Z":[               ],

    "a":["a", 'A'            ],
    "c":[    "c", 'C'        ],
    "g":[        "g", 'G'    ],
    "t":[            "t", 'T', 'U', 'u'],
    "u":[            "u", 'U', 'T', 't'],
    "w":["a", 'A',        "t", 'T'],
    "s":[    "c", 'C',"g", 'G'    ],
    "m":["a", 'A', "c", 'C'        ],
    "k":[        "g", 'G', "t", 'T'],
    "r":["a", 'A',     "g", 'G'   ],
    "y":[    "c", 'C',    "t", 'T'],
    "b":[    "c", 'C', "g", 'G', "t", 'T'],
    "d":["a", 'A',    "g", 'G', "t", 'T'],
    "h":["a", 'A', "c", 'C',    "t", 'T'],
    "v":["a", 'A', "c", 'C', "g", 'G'   ],
    "n":["a", 'A', "c", 'C', "g", 'G', "t", 'T'],
    "z":[               ],
}

########################################################################
# Building list of motifs from input file
motif_list = []

with open(motif_input, "r") as mo:
    for line in mo:
        motif_list.append(line.replace('\n',''))

########################################################################
########## USE MOTIF CONV. FN() TO MAKE DICTIONARY OF REGEX

motif_regex_dict = {}

for m in range(len(motif_list)):
    motif_regex_dict[motif_list[m]] = motif_conversion(motif_list[m])

# motif_regex_dict has keys with original motifs (with ambiguous IUPAC symbols)
# and values with their corresponding reg. expressions 
# built with motif_conversion() fn previously defined

########################################################################
# Creating read and header lists to work from

current_read_list = []
all_reads = []
all_read_headers = []

with open(fasta_input, "r") as fa:
    for line in fa:
        if ">" not in line:
            current_read_list.append(line.replace('\n',''))
        else:
            all_reads.append(''.join(current_read_list))
            current_read_list = []
            all_read_headers.append(line.replace('\n',''))

    all_reads.append(''.join(current_read_list))

# removing first (empty) read
all_reads.pop(0)


# Get list of just the gene symbols of the headers from the FASTA header lines
header_list = [x.split()[0][1:] for x in all_read_headers]



##########################################################
##########################################################
# OUTPUT TO GRAPHIC AND LOOPING THROUGH ALL READS 
# FOR PROCESSING STARTS BELOW HERE
##########################################################

#####################################################3
# Defining color functions
def intron_exon_col():
    ctx.set_source_rgb(0.541176470588235, 0.772549019607843, 1.0)

def red():
    ctx.set_source_rgb((247/255), (35/255), (35/255))

def orange():
    ctx.set_source_rgb((255/255), (140/255), 0)

def yellow():
    ctx.set_source_rgb((255/255), (234/255), 0)

def green():
    ctx.set_source_rgb((77/255), (255/255), 0)

def indigo():
    ctx.set_source_rgb((85/255), (43/255), (117/255))

def violet():
    ctx.set_source_rgb((252/255), (50/255), (209/255))

def dark_green():
    ctx.set_source_rgb((46/255), (102/255), (62/255))

def navy_blue():
    ctx.set_source_rgb((10/255), (92/255), (170/255))

def burgundy():
    ctx.set_source_rgb((128/255), (60/255), (60/255))

def purple():
    ctx.set_source_rgb((158/255), (21/255), (142/255))

def black():
    ctx.set_source_rgb(0, 0, 0)
    ######################################################
    ### Setting up combination function to iterate through color-picking
def col_fn_list(x):
    if x == 0:
        intron_exon_col()
        
    elif x == 1:
        red()
        
    elif x == 2:
        orange()
        
    elif x == 3:
        yellow()
        
    elif x == 4:
        green()
        
    elif x == 5:
        indigo()
        
    elif x == 6:
        violet()
        
    elif x == 7:
        dark_green()
        
    elif x == 8:
        navy_blue()
        
    elif x == 9:
        burgundy()
        
    elif x == 10:
        purple()
        
    elif x == 11:
        black()
    #########################################################
    # col_fn_list() function chooses color based on number - pulls from previously defined color functions
    # 0: for intron/exon
    # 1 - 10 for colors of motifs
    # 11: black



# Setting up cairo surface
num_reads = len(all_reads)
width, height = 1000, (500 * num_reads) + 200

# Setting up output to write to [file prefix of fasta input].svg
surface = cairo.SVGSurface( ("%s.svg" % (fasta_input)), width, height)
ctx = cairo.Context(surface)

# Adding text to annotate motif output

# Setting color, size and style of font
ctx.set_source_rgb(0, 0, 0)
ctx.set_font_size(15)
ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

# Getting dimensions of font output for centering
(x, y, width, height, dx, dy) = ctx.text_extents("Longer vertical bars indicate start of motifs for discerning ambiguous overlapping sites")

# Positioning and outputting text
ctx.move_to( (500 - (width / 2)), ( (500 * num_reads) + 100 - height))
ctx.show_text("Longer vertical bars indicate start of motifs for discerning ambiguous overlapping sites")


# Looping through all the FASTA records, using variable rec (for record)
for rec in range(len(all_reads)):
    aread = all_reads[rec]
    aheader = all_read_headers[rec]

    # STARTING HERE EVERYTHING NEEDS TO BE PACKED INTO USER DEFINED CLASSES

    # Gene Class

    end = int(aheader.split()[1].split(":")[1].split("-")[1])
    start = int(aheader.split()[1].split(":")[1].split("-")[0])
    glen = (end - start)

    agene = Gene(glen)
    agene.start = [0, (150 + (rec * 500))]

    # agene (Gene class) is now initiated
    '''
    ### Drawing intron to overlay exon on

    # Variables for placement
    # Intron dimensions
    intron_tl_x = 0
    intron_tl_y = 150 + (rec * 500)
    intron_width = 1000
    intron_height = 10
    '''

    # FastaHeader.start .txt
    fa_header = FastaHeader(aheader.split()[0][1:])
    print(fa_header.text)

    # Variables for placement
    # header_x_placement is picked based of text_extents, to center text in output
    header_y_placement = 25 + (rec * 500)

    # Getting dimensions of text in order to place with respect to center in .svg graphic output
    (x, y, width, height, dx, dy) = ctx.text_extents("Gene Symbol: " + fa_header.text)

    # Placing and outputting text
    fa_header.start = [(500 - (width / 2)), header_y_placement]

    # OK - so fa_header - FastaHeader class - works
    #   .start: [x, y] placement of header


   

    # Getting positions for exon
    # test_exon_placer_output is a list containing the following: exon sequence, start and end (order-specific)
    test_exon_placer_output = (exon_placer(aread))
    start_pct = (test_exon_placer_output[1] / len(aread)) * 1000
    end_pct = ( (test_exon_placer_output[2] - test_exon_placer_output[1]) / len(aread) ) * 1000

     # EXON OBJECT PORTION - .start, >.length<
    anexon = Exon(end_pct)
    anexon.start = [start_pct, (135 + (rec * 500))]

    # anexon - Exon() class object - works

    '''
    # Exon
    # Variables for placement
    y_exon = 135 + (rec * 500)
    vert_len_exon = 40

    
    # Drawing exon
    col_fn_list(0)
    ctx.rectangle(start_pct, y_exon, end_pct, vert_len_exon)
    ctx.fill()
    '''

    # MOTIF OBJECTS PORTION

    ###############################################################
    # Getting indices from reads - where motifs start
    # Adding to empty (temporary) dictionary

    aread_motif_indices = {}
    
    # Looping through motif_list to check through each motif - adding all to temporary dictionary
    for x in range(len(motif_list)):
        current_motif_indices = [match.start() for match in regex.finditer(motif_regex_dict[motif_list[x]], aread, overlapped = True)]
        aread_motif_indices[motif_list[x]] = current_motif_indices


    agenegroup = GeneGroup(agene, anexon, fa_header, rec)

    # Above outputs "aread_motif_indices" dictionary with lists as indices of original motifs
    # empty lists indicate that the motif was not found in the read at all
    

    working_motif_list = []

    for (k, v) in aread_motif_indices.items():
        if len(v) > 0:
            alist = [k, v]
            working_motif_list.append(alist)

    #print(working_motif_list[-1])
    #print(len(working_motif_list))


    

    for m in range(len(working_motif_list)):
        print(working_motif_list[m])
        mseq = (working_motif_list[m][0])
        mlen = len(working_motif_list[m][0])
        amotif = Motif(mlen, mseq)
        amotif.start = working_motif_list[m][1]

        agenegroup.end_motifs.append(amotif)
    # .start, >.length, .sequence<

    #####################################################################
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #####################################################################
    ## GeneGroup object has now been built

    ######################################################
    ### Drawing intron to overlay exon on

    # Variables for placement
    # Intron dimensions
    intron_tl_x = 0
    intron_tl_y = agenegroup.end_gene.start[1]
    intron_width = 1000
    intron_height = 10


    # Drawing intron
    col_fn_list(0)
    ctx.rectangle(intron_tl_x, intron_tl_y, intron_width, intron_height)
    ctx.fill()

    ##### Drawing Exon

    # Exon
    # Variables for placement
    start_pct = agenegroup.end_exon.start[0]
    y_exon = agenegroup.end_exon.start[1]
    end_pct = agenegroup.end_exon.length
    vert_len_exon = 40

    # Drawing exon
    col_fn_list(0)
    ctx.rectangle(start_pct, y_exon, end_pct, vert_len_exon)
    ctx.fill()

    ###########################
    ## Adding Header

    ######################################################
    ### Adding header of read to top of output

    # Variables for placement
    # header_x_placement is picked based of text_extents, to center text in output
    header_x_placement = fa_header.start[0]
    header_y_placement = fa_header.start[1]

    # Making header text - only need gene name - 
    # splitting on whitespace and leaving out the leading ">" character
    gene_name = fa_header.text

    # Setting color, size, style of font
    col_fn_list(11)
    ctx.set_font_size(20)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)


    # Placing and outputting text
    ctx.move_to(header_x_placement, header_y_placement)
    ctx.show_text("Gene Symbol: " + gene_name)
    
    ######################################

    ### MOTIF PORTION - AND LEGEND
    #################################################################
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    #####################################################################
    # Making motif legend

    # Variables for placement
    legend_text_x = 50
    legend_text_y = 280 + (rec * 500)

    # Color, size, style of font
    col_fn_list(11)
    ctx.set_font_size(15)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    # Positioning and outputting font
    ctx.move_to(legend_text_x, legend_text_y)
    ctx.show_text("Motif Legend")

    # Making motif legend outline
    # based on length of working_motif_list 
    # how many motifs are found in this sequence

    # Setting dimensions of outline
    left = 50
    right = 175
    top = 290 + (rec * 500)
    bottom = top + (len(agenegroup.end_motifs) * 25 )


    # Setting color, line width, choosing path and drawing to .svg output
    col_fn_list(11)
    ctx.set_line_width(1)
    ctx.move_to(left, top)
    ctx.line_to(right, top)
    ctx.line_to(right, bottom)
    ctx.line_to(left, bottom)
    ctx.line_to(left, top)
    ctx.close_path()
    ctx.stroke()

    ###########################

    # Setting variable x for further use positioning individual genes in output graphic
    x = 0 + (rec * 500)

    # Setting loop - variable "m" for: 1 motif, and associated information
    for m in range(len(agenegroup.end_motifs)):
        # Setting color; moving through list of colors in pre-defined swatch
        col_fn_list(m + 1)
        motif_len = agenegroup.end_motifs[m].length
        curr_motif_width = (motif_len / len(aread)) * 1000

        # Inner loop to move through each of the occurrences of finding a motif
        for p in range(len(agenegroup.end_motifs[m].start)):
            
            # Setting position for motif start and scaling to size of graphic output
            raw_posn = agenegroup.end_motifs[m].start[p]
            curr_motif_posn = (raw_posn / len(aread) ) * 1000

            # Setting dimensions for motif
            y_top_main = 130 + (rec * 500)
            y_main_length = 50
            y_top_starting = 100 + (rec * 500)
            y_starting_length = 110


            # Drawing motifs, composed of two parts:
            # Main portion - the entire length of the motif is represented by this rectangle
            ctx.rectangle(curr_motif_posn, y_top_main, curr_motif_width, y_main_length)
            # "Starting" portion - the exact start of the motif is marked by this thin "line" (also a rectangle)
            ctx.rectangle(curr_motif_posn, y_top_starting, 1, y_starting_length)

        # Making swatches of color (squares) in the key legend, 
        # to annotate motif colors with sequences from input

        # Variables for placement
        legend_col_x = 60

        # Drawing squares to output
        ctx.rectangle(legend_col_x, (300 + x), 10, 10)
        ctx.fill()

        # Annotating legend with motif sequences:

        # Variables for placement
        motif_seqs_x = 85

        # Choosing color, size and style of font
        col_fn_list(11)
        ctx.set_font_size(10)
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

        # Writing out individual motifs to their corresponding color swatches in the legend
        # Moving through system in a loop
        ctx.move_to(motif_seqs_x, (310 + x) )
        #ctx.show_text(working_motif_list[m][0])
        ctx.show_text(agenegroup.end_motifs[m].sequence)


        # Incrementing x (used in placement)
        x += 15

        ###############################################################
        ## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>