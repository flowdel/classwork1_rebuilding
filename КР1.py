from Bio import SeqIO
import argparse
import seaborn as sns

def GC_content(file, q):
    
    file = '/Users/Adele/Desktop/file.fastq'
    q = 20
    content = []
    for record in SeqIO.parse(file, 'fastq'):
        GC, AT = 0, 0
        read = record.seq
        qual = record.letter_annotations['phred_quality']
        for i in range(len(qual)):
            if qual[i] < q:
                continue
            elif qual[i] >= q:
                if read[i] == 'G' or read[i] == 'C':
                    GC += 1
                if read[i] == 'A' or read[i] == 'T':
                    AT += 1
                    
        content.append(round((GC/(AT+GC))*100))
    

    ax = sns.distplot(content)
    fig = ax.get_figure()
    fig.savefig('GC_content.png')


        
if __name__ == "__main__":
       
    parser = argparse.ArgumentParser(description='Kmer spectrum')
    
    parser.add_argument('-f', '--file', help='fastq file', type=str, required=True)
    parser.add_argument('-q', '--quality', help='min quality of each nucleotide in fatsq file', type=int, default=20)
    
    args = parser.parse_args()
    file = args.file
    q = args.quality
    
    GC_content(file, q)