# VladTyper

A program for predicting blood types using simple SNVs.

# Inputs

1 - A definition file describing which genotype to use for predicting each antigen:
    file_genotype_definitions = "antigen.info.tsv"

2 - The path and file post/prefix name of bed/bim/fam genotype chromosome files to be used for predictions:
    snpdir = '/genotypes/'
    chr_file = 'chr'
    chr_file_postfix = '_genotypes'

3 - Paths for work directory where predictions and stats are output:
    workdir = '/work/'

3 - Paths for work directory where predictions and stats are output:
    datadir = '/data/'

5 - Optional inptuts are files containing existing blood antigen type data:
    file_abo_rh = 'abo.rh.tsv'
    file_antigens = 'antigens.tsv'

# Options

1 - To predict Sid blood group Sda antigen type set to True (this function is untested):
    PREDICT_SID = False

2 - To predict Lewis blood group Lea and Leb antigen types (these predictions are not very accurate)
    PREDICT_LEWIS = False

# Outputs

1 - Predictions file 'predictions.tsv'
    A file with 4 columns of data. First column denotes ID of a person, second column denotes antigen type, third column denotes serological antigen type status (if there was one loaded from optional files), fourth column is genotype predicted antigen status

2 - statistics file 'predictions.report.tsv'
    A file containing stats on each variant used for prediction

3 - conflicts file 'predictions.conflict.tsv'
    A file containing a list of blood antigen predictions that conflicted with serological loaded data for each person

3 - conflicts file 'predictions.badsamples.tsv'
    A file containing a list of persons with many conflicted predictions, indicating that their genotype sample might have been accidentally swapped with another person

# Notes

As per the example for input 2 the bed/bim/fam files are assumed to have names like chr1_genotypes.bed, chr2_genotypes.bed etc. If there there is nothing after the cromosome indentifier in the file name e.g chr1.bed, chr2.bed, just set chr_file_postfix = ''

The optional input files containing existing blood types can be used to calculate accuracy of predictions. Formatting of each file type refer to the included exmaple files.

The antigen status valuse in the predictions are a negative value for a negative antigen status and a positive value for positive antigen status. Values below -1 and above 1 indicate that the status was confirmed more than once, which is the case if several serological results exist for the same antigen for the same person, or if the antigen type is typed using more than one genetic variant.