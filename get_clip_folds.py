import gtf
import fasta
import rnafold
import bed

genes_parser = gtf.ParseGtf('/Users/martin/Dropbox/utils/genes/hg38/Homo_sapiens.GRCh38.80toy.gtf')
genome_parser = fasta.ParseFasta('/Users/martin/Dropbox/utils/gNome/GRCh38_2.fa')
del_parser = bed.ParseBed('/Users/martin/Dropbox/projects/year2/joel_clip/fold/prova.bed')
rnafold_handler = rnafold.HandleRNAfold()

print 'getting coordinates'
trans_exons = genes_parser.get_trans_exon()

print 'getting sequences'
trans_seqs = genome_parser.get_trans_seqs(trans_exons)

print 'getting foldings'
trans_folds = rnafold_handler.get_plfolds(trans_seqs, wind_size=80)
trans_seq = None

print 'getting folding positions'
pos_fold = {}
for trans in trans_exons.keys():
    exons = trans_exons[trans]
    trans_pos = 0
    for exon in exons:
        direc = int(exon[3] + '1')
        for base in range(int(exon[1]), int(exon[2]) + 1)[::direc]:
            pos = '_'.join([exon[0], str(base), exon[3]])
            if pos in pos_fold:
                pos_fold[pos].append(trans_folds[trans][trans_pos]) #checkthis!!!
            else:
                pos_fold[pos] = [trans_folds[trans][trans_pos]]
            trans_pos += 1

for pos, folds in pos_fold.items():
    pos_fold[pos] = sum(folds) / len(folds)