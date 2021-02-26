#! /usr/bin/env python3

import sys

#1 DONE!!! Passed 9/10

def translate_sequence(rna_sequence, genetic_code):
    pass
# Read and get the RNA string
    rna = rna_sequence.upper()
    print ("\n \n RNA String: ", rna)
    x=len(rna)
    if x < 3:
        return ''
# RNA codon table(make sure you have it)
    protein_string = ""
# Generate protein string
    for i in range(0, len(rna),3):
        if genetic_code[rna[i:i+3]] == "*" :
            break
        protein_string += genetic_code[rna[i:i+3]]
    return protein_string
#    return protein_string
	# Print the protein string
    print ("\n \n Protein String: ", protein_string)

	# End

#2 Passed 2/4 All giving protein but not passing 
def get_all_translations(rna_sequence, genetic_code):
    pass

# Read and get the RNA string
    DNA = rna_sequence.upper()

    print ("\n \n RNA String1: ", DNA)

    if (DNA.find('AUG') != -1):
        pass
#        print ("Contains given substring ")
    else:
        return []
#        print ("Doesn't contains given substring")

    start = DNA.find('AUG')

    protein_string1 = ""

    if start!= -1:
        while start+2 < len(DNA):
            codon = DNA[start:start+3]
            if genetic_code[codon] == "*":
                break
            protein_string1 += genetic_code[codon]
            return [protein_string1]
            start+=3
#        print ("\n \n Protein String1: ", protein_string1)

    DNA2= DNA[1:]
    print ("\n \n RNA2 String: ", DNA2)


    if (DNA2.find('AUG') != -1):
        pass
#        print ("Contains given substring ") 
    else: 
        return []
#        print ("Doesn't contains given substring") 

    start = DNA2.find('AUG')
    protein_string2 = ""

    if start!= -1:
        while start+2 < len(DNA2):
            codon = DNA2[start:start+3]
            if genetic_code[codon] == "*":
                break
            protein_string2 += genetic_code[codon]
            return [protein_string2]
            start+=3
#        print ("\n \n Protein String2: ", protein_string2)

    DNA3= DNA[2:]
    print ("\n \n RNA3 String: ", DNA3)

    if (DNA3.find('AUG') != -1):
        pass
#        print ("Contains given substring ") 
    else: 
        return []
#        print ("Doesn't contains given substring") 

    start = DNA3.find('AUG')
    protein_string3 = ""

    if start!= -1:
        while start+2 < len(DNA3):
            codon = DNA3[start:start+3]
            if genetic_code[codon] == "*":
                break
            protein_string3 += genetic_code[codon]
            return [protein_string3]
            start+=3
#        print ("\n \n Protein String3: ", protein_string3)

#3 DONE Passed All

def get_reverse(sequence):
    pass
    sequence = sequence.upper()
    re = []
    x = len(sequence)
    for i in sequence:
        x = x - 1
        re.append(sequence[x])
    return ''.join(re)

#4 DONE Passed All

def get_complement(sequence):
    pass

    sequence = sequence.upper()
    com = []
    for i in sequence:
        if i == "U":
            com.append("A")
        if i == "A":
            com.append("U")
        if i == "G":
            com.append("C")
        if i == "C":
            com.append("G")

    return ''.join(com)

#5 DONE Passed All

def reverse_and_complement(sequence):
    pass

    sequence = sequence.upper()
    re = []
    x = len(sequence)
    for i in sequence:
        x = x - 1
        re.append(sequence[x])
#    return ''.join(re)
    com = []
    for i in re:
        if i == "U":
            com.append("A")
        if i == "A":
            com.append("U")
        if i == "G":
            com.append("C")
        if i == "C":
            com.append("G")

    return ''.join(com)

#6

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    pass

# Read and get the RNA string
    DNA = rna_sequence.upper()
    print ("\n \n RNA String1: ", DNA)
    if (DNA.find('AUG') != -1):
        pass
    else:
        return ""
    DNA2= DNA[1:]
#    print ("\n \n RNA2 String: ", DNA2)
    if (DNA2.find('AUG') != -1):
        pass
    else:
        return ""


    DNA3= DNA[2:]
#    print ("\n \n RNA3 String: ", DNA3)

    if (DNA3.find('AUG') != -1):
        pass
    else:
        return ""


    sequence = DNA[1:]

    # Find orf 2
    # Find all AUG indexs
    start_position = 1
    start_indexs = []
    stop_indexs = []
    for i in range(1, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(1, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)
    orf = []
    mark = 0
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                mark = stop_indexs[j]+3
                break
#    return orf
    orf2 = orf
    print(orf2)
    pass
    sequence = DNA[2:]
	# Find all ATG indexs
    start_position = 2
    start_indexs = []
    stop_indexs = []
    for i in range(2, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(2, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)

    orf = []
    mark = 0
    start_position = {}
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                start_position[len(sequence[start_indexs[i]:stop_indexs[j]+3])] = start_indexs[i]
                mark = stop_indexs[j]+3
                break
#    return orf
    orf3 = orf
    print(orf3)
    pass
    sequence = DNA[0:]
    start_position = 0
    start_indexs = []
    stop_indexs = []
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(0, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)

    orf = []
    mark = 0
    start_position = {}
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                start_position[len(sequence[start_indexs[i]:stop_indexs[j]+3])] = start_indexs[i]
                mark = stop_indexs[j]+3
                break
#    return orf
    orf1 = orf
    print(orf1)
    pass
    print ("\n \n RNA1 String: ", orf2)


#list2str
    orf2_str =''.join([str(elem) for elem in orf2])
    print("\n orf2_str", orf2_str)

    start = orf2_str.find('AUG')
#    print(start)
    protein_string = ""
    if start!= -1:
        while start+2 < len(orf2_str):
            codon = orf2_str[start:start+3]
            if genetic_code[codon] == "*":
                break
            protein_string += genetic_code[codon]
#            return protein_string
#            print ("\n \n Protein String: ", protein_string)
            start+=3
    orf2_protein=protein_string
    print ("\n longest Forward peptide:", orf2_protein)


###Reverse sequence longest peptide
    DNA = rna_sequence.upper()
    sequence = DNA
    re = []
    x = len(sequence)
    for i in sequence:
        x = x - 1
        re.append(sequence[x])
#    return ''.join(re)
    com = []
    for i in re:
        if i == "U":
            com.append("A")
        if i == "A":
            com.append("U")
        if i == "G":
            com.append("C")
        if i == "C":
            com.append("G")

#    return ''.join(com)
    rDNA=''.join(com)
    print ("\n \n RNA String_Reverse1: ", rDNA)
    if (rDNA.find('AUG') != -1):
        pass
    else:
        pass
#        return ""


    rDNA2= rDNA[1:]
    print ("\n \n RNA String_reverse2: ", rDNA2)
    if (rDNA2.find('AUG') != -1):
        pass
    else:
#        return ""
        pass

    rDNA3= rDNA[2:]
    print ("\n \n RNA3 String_reverse3: ", rDNA3)

    if (rDNA3.find('AUG') != -1):
        pass
    else:
#        return ""
        pass

    sequence = rDNA[1:]

    # Find orf 2
    # Find all AUG indexs
    start_position = 1
    start_indexs = []
    stop_indexs = []
    for i in range(1, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(1, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)
    orf = []
    mark = 0
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                mark = stop_indexs[j]+3
                break
#    return orf
    reverse_orf2 = orf
    print("\nORF2reverse: ",reverse_orf2)
    pass
    sequence = rDNA[2:]
	# Find all ATG indexs
    start_position = 2
    start_indexs = []
    stop_indexs = []
    for i in range(2, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(2, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)

    orf = []
    mark = 0
    start_position = {}
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                start_position[len(sequence[start_indexs[i]:stop_indexs[j]+3])] = start_indexs[i]
                mark = stop_indexs[j]+3
                break
#    return orf
    reverse_orf3 = orf
    print("\nORF3reverse: ",reverse_orf3)
    pass
    sequence = rDNA[0:]
    start_position = 0
    start_indexs = []
    stop_indexs = []
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] == "AUG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(0, len(sequence), 3):
        stops =["UAA", "UGA", "UAG"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)

    orf = []
    mark = 0
    start_position = {}
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                start_position[len(sequence[start_indexs[i]:stop_indexs[j]+3])] = start_indexs[i]
                mark = stop_indexs[j]+3
                break
#    return orf
    reverse_orf1 = orf
    print("\nORF1reverse: ",reverse_orf1)
    pass
#    print ("\n \n RNA1 String_reverse: ", orf2)

#list2str
    orf2_str =''.join([str(elem) for elem in orf2])
    print("\n orf2_str", orf2_str)

    start = orf2_str.find('AUG')
#    print(start)
    protein_string = ""
    if start!= -1:
        while start+2 < len(orf2_str):
            codon = orf2_str[start:start+3]
            if genetic_code[codon] == "*":
                break
            protein_string += genetic_code[codon]
#            return protein_string
#            print ("\n \n Protein String: ", protein_string)
            start+=3
    orf2_protein=protein_string
    print ("\n longest Reverese peptide:", orf2_protein)


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
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
