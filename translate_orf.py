#! /usr/bin/env python3

import sys
import re
#import find_orf
#import translate

#1 Passed 6/6 test :) 
def vet_nucleotide_sequence(sequence):
    """
    Return None if `sequence` is a valid RNA or DNA sequence, else raise exception. 

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if sequence is valid, otherwise raise an
        exception.

    Examples
    --------
    >>> vet_nucleotide_sequence('ACGTACGT') == None
    True

    >>> vet_nucleotide_sequence('not a valid sequence')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'not a valid sequence'

    Don't allow mixing of DNA and RNA!
    >>> vet_nucleotide_sequence('AUTGC')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'AUTGC'

    Don't allow whitespace (or other characters) before, within, or after!
    >>> vet_nucleotide_sequence(' ACGT ACGT ')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: ' ACGT ACGT '

    But, an empty string should be deemed valid
    >>> vet_nucleotide_sequence('') == None
    True
    """
    ##########################################################################
    ############################ EDIT CODE BELOW #############################
    # `rna_pattern_str` and `dna_pattern_str` need to be regular expressions
    # that will match any string of zero or more RNA and DNA bases,
    # respectively (and only strings of zero or more RNA and DNA bases).
    # Currently, `rna_pattern_str` and `dna_pattern_str` are strings of literal
    # characters.
    # These are valid regular expressions, but they will only match their
    # respective strings exactly.
    # Change `rna_pattern_str` and `dna_pattern_str` so that they will match
    # any valid RNA and DNA sequence strings, respectively (and only strings of
    # RNA and DNA bases).
    # Read the docstring above for additional clues.
    sequence = sequence.upper()
    rna_pattern_str = r'^[AUGC]+$|^$'
    dna_pattern_str = r'^[ATGC]+$|^$'
    ##########################################################################

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))

#test2 Pass 5/5 :) 

def vet_codon(codon):
    """
    Return None if `codon` is a valid RNA codon, else raise an exception. 

    Parameters
    ----------
    codon : str
        A string representing a codon (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if codon is valid, otherwise raise an
        exception.

    Examples
    --------
    Valid codon
    >>> vet_codon('AUG') == None
    True

    lower-case is also vaild 
    >>> vet_codon('aug') == None
    True

    DNA is not valid
    >>> vet_codon('ATG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'ATG'

    A codon must be exactly 3 RNA bases long
    >>> vet_codon('AUGG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'AUGG'
    """
    ##########################################################################
    ############################ EDIT CODE BELOW #############################
    # `codon_pattern_str` needs to be a regular expression that will match any
    # codon (but only a string that is one codon).
    # Currently, `codon_pattern_str` is only a string of literal characters.
    # This is a valid regular expression, but it will only match 'AUG' exactly.
    # Change `codon_pattern_str` so that it will match any valid codons, and
    # only valid codons.
    # Read the docstring above for additional clues.
    codon =codon.upper()
    codon_pattern_str = r'^([AUGC]{3})$'
    
	##########################################################################

    codon_pattern = re.compile(codon_pattern_str)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))

#test3 11/13 passed :) 

def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    """
    Return the first open-reading frame in the DNA or RNA `sequence`.

    An open-reading frame (ORF) is the part of an RNA sequence that is
    translated into a peptide. It must begin with a start codon, followed by
    zero or more codons (triplets of nucleotides), and end with a stop codon.
    If there are no ORFs in the sequence, an empty string is returned.

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)
    start_codons : list of strings
        All possible start codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.
    stop_codons : list of strings
        All possible stop codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.

    Returns
    -------
    str
        An uppercase string of the first ORF found in the `sequence` that
        starts with any one of the `start_codons` and ends with any one of the
        `stop_codons`. If no ORF is found an empty string is returned.

    Examples
    --------
    When the whole RNA sequence is an ORF:
    >>> find_first_orf('AUGGUAUAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When the whole DNA sequence is an ORF:
    >>> find_first_orf('ATGGTATAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When there is no ORF:
    >>> find_first_orf('CUGGUAUAA', ['AUG'], ['UAA'])
    ''

    When there is are bases before and after ORF:
    >>> find_first_orf('CCAUGGUAUAACC', ['AUG'], ['UAA'])
    'AUGGUAUAA'
    """
    # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    # Get copies of everything in uppercase
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    # Make sure seq is RNA
    seq = seq.replace('T', 'U')

    ##########################################################################
    ############################ EDIT CODE BELOW #############################
    # `orf_pattern_str` needs to be a regular expression that will match an
    # open reading frame within a string of RNA bases. At this point we know
    # the string only contains uppercase A, C, G, and U.
    # I recommend starting by hardcoding the standard start and stop codons
    # (the ones listed as defaults for this function) into the regular
    # expression. After you get that working, then try generalizing it to work
    # for any start/stop codons.
    # Currently, `orf_pattern_str` is only a string of literal characters. This
    # is a valid regular expression, but it will only match 'AUGGUAUAA'
    # exactly. Change `orf_pattern_str` so that it will match any open reading
    # frame.
    # Read the docstring above for additional clues.
    orf_pattern_str = r'AUG(...)*U(AA|AG|GA)'
    ##########################################################################

    # Create the regular expression object
    orf_pattern = re.compile(orf_pattern_str)
    # Search the sequence
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''


#By jaime oaks
#my one is to big and failed one test I think my one will work as well. :) 
#passed 10/10 tests
def pop_next_codon(sequence):
	codon =sequence[0:3]
	remaining_seq = sequence[3:]
	return codon, remaining_seq

def translate_sequence(rna_sequence, genetic_code):
	rna_sequence = rna_sequence.upper()
	amino_acid_list = []
	while True:
		if len(rna_sequence) < 3:
			break
		codon, remaining_seq = pop_next_codon(rna_sequence)
		rna_sequence = remaining_seq
		aa = genetic_code[codon]
		if aa == "*":
			break
		amino_acid_list.append(aa)
	return "".join(amino_acid_list)


def main():
    import argparse

    # Create a command-line parser object
    parser = argparse.ArgumentParser()

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']
#    default_genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}

    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(" ".join(default_stop_codons))))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    # Check to see if start/stop codons were provided by the caller. If not,
    # use the defaults.
    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons
#    if not args.genetic_code:
#        args.genetic_code = default_genetic_code

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)
    sys.stdout.write('{}\n'.format(orf))

    protein = translate_sequence(rna_sequence = orf, genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'})
    sys.stdout.write('{}\n'.format(protein))

if __name__ == '__main__':
    main()
if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")


