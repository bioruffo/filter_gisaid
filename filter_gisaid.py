#!/usr/bin/python3

from collections import defaultdict

def read_metadata(metafile):
    # Return a dict of {strain_id: lineage}
    header = ''
    lineage = 0
    strain = 0
    
    metadata = dict() # ex. {'Brazil/ldfd/2021': 'P.1'}
    
    print("Reading metadata...")
    for i, line in enumerate(open(metafile, "r", encoding="utf-8")):
        line = line.strip().split('\t')
        if line[:2] == ["strain","virus"]:
            header = line
            lineage = header.index('pangolin_lineage')
            strain = header.index('gisaid_epi_isl')
        else:
            assert line[strain] not in metadata
            metadata[line[strain]] = line[lineage]
    return metadata
        


def read_seqs(seqfile, strain_names, metadata):
    strains = {name: defaultdict(int) for name in strain_names}
    next_is_header = True
    
    rejected = 0
    total = 0
    for k, line in enumerate(open(seqfile, 'r', encoding="utf-8")):
        line = line.strip()
        if line.startswith(">") and next_is_header:
            total += 1
            idstring = line.split("|")[3]
            # we only want the strains we selected
            this_strain = metadata.get(idstring, None)
            if this_strain is None:
                rejected += 1
            next_is_header = False
        else:
            if not next_is_header:
                if this_strain in strains:
                    if line.count('X')/len(line) < max_x:
                        strains[this_strain][line] += 1
                next_is_header = True
            else:
                print("oops!")
                break
    print("Rejected {}/{} sequences".format(rejected, total))
    return strains
            

def output_seqs(strains, appendseq, minseq):
    print("Output...")
    if minseq > 0:
        print("Excluding sequences with {} count(s) or less".format(minseq))
    for strain_name, strain_var in strains.items():
        filename = strain_name + "_summ.fasta"
        print("Writing:", filename)
        with open(filename, "w", encoding="utf-8") as f:
            f.write(appendseq)
            suffix = 1
            for seq, count in sorted(strain_var.items(), key = lambda x: x[1], reverse=True):
                if count > minseq:
                    f.write(">({}).{}\n{}\n".format(count, suffix, seq))
                    suffix += 1
        


if __name__ == '__main__':
    
    # Metadata file
    metafile = 'data/metadata.tsv'
    seqfile = 'data/spikeprot0308.fasta'
    
    # Strain names
    strain_names = ['P.1', 'B.1.1.7', 'B.1.351']
    
    # Maximum number of indeterminate residues (X)
    max_x = 0.01
    
    # Minimum count of sequences to be written to output; set to 0 for all sequences
    minseq = 1
    
    # Append anything to the files?
    appendseq = '\n'.join(x.strip() for x in '''
    >gi|1796318598|ref|YP_009724390.1|:1-1273 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]
    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
    '''.strip().split('\n'))+'\n'
    
    metadata = read_metadata(metafile)
    strains = read_seqs(seqfile, strain_names, metadata)
    output_seqs(strains, appendseq, minseq)
   
    
