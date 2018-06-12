import gj
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
sns.set(style="darkgrid")
sns.set_context("poster")
import sys
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta

def read_fa(fa=None):
	if fa is None:
		fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def read_tmp_out(tmp_out=None,file_str=None,fa=None):
	""" caculate base enriched ratio
	"""
	file_str = file_str if file_str is not None else tmp_out.replace(".out", ".bass_count.txt")
	gj.printFuncRun('read_tmp_out')
	gj.printFuncArgs()
	fa_dict = read_fa(fa)
	tx_base_pos_dict = nested_dict(2, list) # {tx:{'A':[pos1,pos2],'T':[]}}
	base_enrich_dict = nested_dict(1, int)
	with open(tmp_out, 'r') as TMP_OUT:
		for line in TMP_OUT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split('\t')
			transcript_id = arr[0]
			transcript_len = int(arr[1])
			if transcript_len != len(fa_dict[transcript_id]):
				print "transcirpt length not conistent with reference: %s, tmp_out len: %s, reference len: %s"%(transcript_id, transcript_len, len(fa_dict[transcript_id]))
				sys.exit()
			for n,base_enrichment_score in enumerate(arr[4:]):
				score = base_enrichment_score.split(',')[0]
				#if score != "NULL" and float(score) != 0 and float(score) >= 0.3:
				if score != "NULL" and float(score) != 0:
					base = fa_dict[transcript_id][n]
					tx_base_pos_dict[transcript_id][base].append(n)
					base_enrich_dict[base] += 1
	print base_enrich_dict

	#val_ls = [base_enrich_dict[i] for i in ['A','T','C','G']]
	#gj.plot_ls_pie(labels=['A','T','C','G'],val=val_ls,dic="",title_str="",file_str=file_str)
	TXT = open(file_str, 'w')
	for i,j in base_enrich_dict.items():
		print >>TXT,i+'\t'+str(j)
	TXT.close()

	gj.printFuncRun('read_tmp_out')

if __name__ == "__main__":
	read_tmp_out(tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.tmp.out', fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa')
	#read_fa()
