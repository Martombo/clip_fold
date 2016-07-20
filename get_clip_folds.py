import parsers as ps
import handlers as hn
import functions as fn

def msg(msg_text):
    with open('msg', 'a') as fmsg:
        fmsg.write(msg_text + '\n')

genes_parser = ps.Gtf('Homo_sapiens.GRCh38.81.gtf')
genome_parser = ps.Fasta('../GRCh38_2.fa')
del_parser_wt = ps.Bed('wt_deletions.bed')
del_parser_mut = ps.Bed('mut_deletions.bed')
out_wt = 'results_wt'
out_mut = 'results_mut'
rnafold_handler = hn.RnaFold(msg = msg)

msg('getting coordinates')
trans_exons = genes_parser.get_trans_exon()

msg('getting sequences')
trans_seqs = genome_parser.get_trans_seqs(trans_exons)

msg('getting foldings')
trans_folds = rnafold_handler.trans_plfolds(trans_seqs, wind_size=80)
trans_seqs = None

msg('converting coordinates to numpy')
exons_trans = genes_parser.trans_exon2np(trans_exon = trans_exons)

msg('getting wt clip sites')
dels_wt = del_parser_wt.get_first()

msg('intersecting wt')
counter = fn.Intersecter.FoldsCounter(trans_folds, trans_exons, span = 200)
inter = fn.Intersecter(dels_wt, exons_trans, counter, strand = False)
inter.intersect()
results_wt = counter.get_results()

with open(out_wt, 'w') as fout:
    for site in results_wt:
        fout.write(' '.join([str(x) for x in [site[1]] + site[0]]) + '\n')

msg('getting wt clip sites')
dels_mut = del_parser_mut.get_first()

msg('intersecting wt')
counter.restart()
inter = fn.Intersecter(dels_mut, exons_trans, counter, strand = False)
inter.intersect()
results_wt = counter.get_results()

with open(out_mut, 'w') as fout:
    for site in results_wt:
        fout.write(' '.join([str(x) for x in [site[1]] + site[0]]) + '\n')
