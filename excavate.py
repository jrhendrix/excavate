'''
FILE:   excavate.py
AUTHOR: Jo Hendrix
URL:    http://stronglab.org
DESC:   Parse megalodon output
		lists positions with modifications
'''

# IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

class stat:
	""" Object for basic stats """

	def __init__(self, arr):
		self.num = len(arr)

		if self.num < 1:
			self.avg = 'NA'
			self.med = 'NA'
			self.std = 'NA'
			self.min = 'NA'
			self.max = 'NA'
			return

		self.avg = np.average(arr)
		self.med = np.median(arr)
		self.std = np.std(arr)
		self.min = np.min(arr)
		self.max = np.max(arr)


class base:
	""" Object for each potentially modified base """

	def __init__(self, pos, contig):
		self.pos = pos
		self.contig = contig
		#self.methcalls = []
		#self.modprobs = []
		#self.canprobs = []

		
		self.calls_hpos = []
		self.calls_hneg = []
		self.calls_mpos = []
		self.calls_mneg = []
		self.calls_cpos = []
		self.calls_cneg = []
		# Alt
		#self.calls_pos = []
		#self.calls_neg = []

		self.probs_hpos = []
		self.probs_hneg = []
		self.probs_mpos = []
		self.probs_mneg = []
		self.probs_cpos = []
		self.probs_cneg = []
		
		self.nr = None
		self.hp = None
		self.hn = None
		self.mp = None
		self.mn = None
		self.cp = None
		self.cn = None

		self.hc = None
		self.mc = None


	# Add methylation status calls
	def add_call_hpos(self, value):
		return self.calls_hpos.append(value)
	def add_call_hneg(self, value):
		return self.calls_hneg.append(value)
	def add_call_mpos(self, value):
		return self.calls_mpos.append(value)
	def add_call_mneg(self, value):
		return self.calls_mneg.append(value)
	def add_call_cpos(self, value):
		return self.calls_cpos.append(value)
	def add_call_cneg(self, value):
		return self.calls_cneg.append(value)
	# Alt
	'''
	def add_call_pos(self, value):
		return self.calls_pos.append(value)
	def add_call_neg(self, value):
		return self.calls_neg.append(value)
	'''

	# Add probability calls
	def add_prob_hpos(self, value):
		return self.probs_hpos.append(value)
	def add_prob_hneg(self, value):
		return self.probs_hneg.append(value)
	def add_prob_mpos(self, value):
		return self.probs_mpos.append(value)
	def add_prob_mneg(self, value):
		return self.probs_mneg.append(value)
	def add_prob_cpos(self, value):
		return self.probs_cpos.append(value)
	def add_prob_cneg(self, value):
		return self.probs_cneg.append(value)

	def set_stats(self):
		self.nrp = len(self.calls_mpos)
		self.nrn = len(self.calls_mneg)
		self.chp = stat(self.calls_hpos)
		self.chn = stat(self.calls_hneg)
		self.cmp = stat(self.calls_mpos)
		self.cmn = stat(self.calls_mneg)
		self.ccp = stat(self.calls_cpos)
		self.ccn = stat(self.calls_cneg)

		self.php = stat(self.probs_hpos)
		self.phn = stat(self.probs_hneg)
		self.pmp = stat(self.probs_mpos)
		self.pmn = stat(self.probs_mneg)
		self.pcp = stat(self.probs_cpos)
		self.pcn = stat(self.probs_cneg)

		#self.hc = stat(self.h_calls)
		#self.mc = stat(self.m_calls)

		#self.pchp = (self.chp.count(1)/self.nrp)*100
		#self.pchn = (self.chn.count(1)/self.nrn)*100

		#self.pcmp = (self.cmp.count(1)/self.nrp)*100
		#self.pcmn = (self.cmn.count(1)/self.nrn)*100
		#self.prm = (self.mc.count(1)/self.nr)*100


def tabMS(key, arr_h, arr_m):

	issig = False

	num1_h = arr_h.count(1)
	num0_h = arr_h.count(0)
	numn_h = arr_h.count(-1)

	num1_m = arr_m.count(1)
	num0_m = arr_m.count(0)
	numn_m = arr_m.count(-1)

	alln_h = num1_h + num0_h + numn_h
	alln_m = num1_m + num0_m + numn_m

	if alln_h != alln_m:
		print('do not add up: ', alln_h, ' 5hmC vs. ', alln_m, ' 5mC')
		exit()

	if alln_m < 1:
		prc_h = 0
		prc_m = 0
	else:
		prc_h = round((num1_h/alln_h)*100, 2)
		prc_m = round((num1_m/alln_m)*100, 2)

	if prc_h >= args.threshold or prc_m >= args.threshold:
		issig = True

	entry = '\t'.join((str(alln_m), 
		str(num1_h), str(num0_h), str(numn_h), str(prc_h),
		str(num1_m), str(num0_m), str(numn_m), str(prc_m)))

	return entry, issig

def prosDB(args, key, ob):
	sig = False
	ob.set_stats()
	c = ob.contig

	pos, sig_pos = tabMS(key, ob.calls_hpos, ob.calls_mpos)
	neg, sig_neg = tabMS(key, ob.calls_hneg, ob.calls_mneg)

	meth_stat = '\t'.join((key, str(c), pos, neg))
	meth_stat = ''.join((meth_stat, '\n'))

	prob_lists = 'NA'
	if sig_pos | sig_neg:
		sig = True

	return meth_stat, prob_lists, sig

def make_mod_entry(seqid, count, start, end, score, strand, mod):
	
	# Set seqID
	## Left pad the string with 0 to make it's length 4
	count = str(count).zfill(4)
	featid = '_'.join(('mod', str(count)))

	# Set types
	source 	= 'Megalodon'
	feat	= 'mod'
	
	# Phase does not apply
	phase	= '.'

	# Set attributes
	iname 	= '='.join(('ID', featid))
	aname 	= '='.join(('Name', mod)) # Add modification mark
	attr  	= ';'.join((iname, aname))

	data = (seqid, source, feat, str(start), str(end), str(score), strand, phase, attr)
	entry = '\t'.join(data) + '\n'

	return entry

def export(args):

	# OPEN INPUT TABLE
	try:
		f1 = open(args.input_file, 'r')
	except:
		print('ERROR: Could not open input file. Exit.')
		exit()

	# CONFIGURE OUTPUT
	#try:
	outdir = '/'.join((args.p, args.output_directory))
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	prefix = '/'.join((outdir, args.savename))
	outgff = '_'.join((prefix, 'significant_mods.gff'))
	# Summary full - avg

	# open files
	f2 = open(outgff, 'w')
	header = "##gff-version 3\n"
	f2.write(header)

	count = 1
	f1.readline()	# skip header line
	for line in f1:

		# INPUT INFORMATION
		line = line.strip()
		l = line.split('\t')
		position 	= int(l[0])
		contig 		= l[1]
		pos_num 	= int(l[2])
		pos_meth_h	= int(l[3])
		pos_na_h	= int(l[4])
		pos_con_h	= int(l[5])
		pos_prc_h	= float(l[6])
		pos_meth_m	= int(l[7])
		pos_na_m	= int(l[8])
		pos_con_m	= int(l[9])
		pos_prc_m	= float(l[10])
		neg_num		= int(l[11])
		neg_meth_h	= int(l[12])
		neg_na_h	= int(l[13])
		neg_con_h	= int(l[14])
		neg_prc_h	= float(l[15])
		neg_meth_m	= int(l[16])
		neg_na_m	= int(l[17])
		neg_con_m	= int(l[18])
		neg_prc_m	= float(l[19])

		# SET OUTPUT INFORMATION
		# Sequence ID
		seqid 	= contig
		
		# Set positions
		start	= position
		end		= position #+ 1

		# Calculate score
		pos_h = int(pos_meth_h) + int(pos_na_h) + int(pos_con_h)
		pos_m = int(pos_meth_m) + int(pos_na_m) + int(pos_con_m)
		neg_h = int(neg_meth_h) + int(neg_na_h) + int(neg_con_h)
		neg_m = int(neg_meth_m) + int(neg_na_m) + int(neg_con_m)

		# CREATE ENTRY
		if pos_h > 0:
			score = (pos_meth_h/pos_h)*100
			if score >= args.threshold:
				entry = make_mod_entry(seqid, count, start, end, score, '+', 'h')
				f2.write(entry)
				count = count + 1
		if pos_m > 0:
			score = (pos_meth_m/pos_m)*100
			if score >= args.threshold:
				entry = make_mod_entry(seqid, count, start, end, score, '+', 'm')
				f2.write(entry)
				count = count + 1
		if neg_h > 0:
			score = (neg_meth_h/neg_h)*100
			if score >= args.threshold:
				entry = make_mod_entry(seqid, count, start, end, score, '-', 'h')
				f2.write(entry)
				count = count + 1
		if neg_m > 0:
			score = (neg_meth_m/neg_m)*100
			if score >= args.threshold:
				entry = make_mod_entry(seqid, count, start, end, score, '-', 'm')
				f2.write(entry)
				count = count + 1

	f1.close()
	f2.close()

def readDB(args):

	# OPEN INPUT DATABASE FILE - table format
	try:
		f1 = open(args.input_file, 'r')
	except:
		print('ERROR: Could not open input file. Exit.')
		exit()

	# CONFIGURE OUTPUT
	try:
		outdir = '/'.join((args.p, args.output_directory))
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		prefix = '/'.join((outdir, args.savename))
		outmsf = '_'.join((prefix, 'mod_stat_full.tsv')) # mod stats - full
		outmss = '_'.join((prefix, 'mod_stat_sig.tsv')) # mod stats - sig
		#outgff = '_'.join((prefix, 'mod _sig.gff'))

		# open files
		f2 = open(outmsf, 'w') # mod stat full
		f3 = open(outmss, 'w') # mod stat sig
		#f4 = open(outgff, 'w') # mod gff file for sig

	except:
		print('ERROR: Could not configure output. Exit.')
		exit()

	# WRITE OUTPUT HEADER
	headerMS = '\t'.join(('position', 'contig',
		'pos_strand_total', 
		'pos_meth_h', 'pos_na_h', 'pos_con_h', 'pos_prc_h', 
		'pos_meth_m', 'pos_na_m', 'pos_con_m', 'pos_prc_m', 
		'neg_strand_total', 
		'neg_meth_h', 'neg_na_h', 'neg_con_h', 'neg_prc_h', 
		'neg_meth_m', 'neg_na_m', 'neg_con_m', 'neg_prc_m'))
	headerMS = ''.join((headerMS, '\n'))
	f2.write(headerMS)
	f3.write(headerMS)

	db = {}
	f1.readline() # skip header line
	curpos = 0
	count = 0
	found = 0
	for l in f1:
		count = count + 1

		# READ DATA
		line = l.strip().split('\t')
		readname			= line[0]
		chrm				= line[1]
		strand				= line[2]
		pos					= int(line[3])
		mod_prob			= np.exp(float(line[4]))*100
		can_prob			= np.exp(float(line[5]))*100
		modname				= line[6]

		# ASSIGN MOD STATE
		## NOTE: Megalodon gives a probability, NOT a call
		metcall = 0 					# Set to ambiguous
		if mod_prob >= args.threshold:
			methcall = 1 				# Set to modified
		if can_prob >= args.threshold:
			methcall = -1 				# Set to cononical

		# ASSIGN KEY
		key = str(pos)

		# ADJUST CURRENT POSITION
		if count == 1:
			curpos = pos

		# PROCESS PREVIOUS POSITION - if moving to next position
		if pos != curpos:
			#print('ONTO NEXT POS')
			# Note: using '>' ignores contigs
			if len(db) < 1:
				print('\tdb was empty. next')
				curpos = pos
			else:
				# PROCESS DB
				prevkey = str(curpos)
				meth_stat, prob_lists, issig = prosDB(args, prevkey, db[prevkey])
				f2.write(meth_stat)

				if issig:
					if chrm != '1':
						print('!! SIGNIFICANT !!')
						print(meth_stat)
					f3.write(meth_stat)

			# RESET DB
			db = {}
			curpos = pos

		# HANDLING CURRENT POSITION
		if key not in db:
			db[key] = base(pos, chrm)
		if modname == 'm':
			if strand == '+':
				db[key].add_call_mpos(methcall)
				db[key].add_prob_mpos(mod_prob)

				db[key].add_prob_cpos(can_prob)
			elif strand == '-':
				db[key].add_call_mneg(methcall)
				db[key].add_prob_mneg(mod_prob)

				db[key].add_prob_cneg(can_prob)
			else:
				print('....strand not found for m.....')
		elif modname == 'h':
			if strand == '+':
				db[key].add_call_hpos(methcall)
				db[key].add_prob_hpos(mod_prob)
			elif strand == '-':
				db[key].add_call_hneg(methcall)
				db[key].add_prob_hneg(mod_prob)
			else:
				print('.... strand not found for h.....')
		else:
			print('MODNAME NOT FOUND')

	# HANDLE LAST POSITION
	prevkey = str(curpos)
	meth_stat, prob_lists, issig = prosDB(args, prevkey, db[prevkey])
	f2.write(meth_stat)

	if issig:
		f3.write(meth_stat)

	# CLEAN UP
	f1.close()
	f2.close()
	f3.close()


	# EXPORT GFF FILE
	if not args.do_not_export_gff:
		print('Creating gff file.')
		args.input_file = outmss # Change input file
		export(args)

	print('DONE')		

def main(args):
	print('Running command: %s' % ' '.join(sys.argv))
	args.func(args)


if __name__ ==  "__main__":

	parser = argparse.ArgumentParser(description='excavate')
	subparsers = parser.add_subparsers(title="tool", dest="tool", help='available actions')
	subparsers.required = True


	# PARSER : ROOT
	__version__ = "0.0.0"
	parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))


	cwd = os.getcwd()

	# DEFINE SUBPARSERS
	# Todo: make parent parser for inheritance
	parser_readDB = subparsers.add_parser('readDB')
	parser_export = subparsers.add_parser('export')

	parser_readDB.set_defaults(func=readDB)
	parser_export.set_defaults(func=export)

	# PARSER : readDB
	# Read a Megalodon database
	parser_readDB.add_argument('-i', '--input_file', help='Path to input file.', required = True)
	parser_readDB.add_argument('-g', '--do_not_export_gff', default=False, action='store_true')
	parser_readDB.add_argument('-o', '--output_directory', help='Name of output directory', default='excavate_out', type=str)
	parser_readDB.add_argument('-p', default=cwd, help='Path to output', type=str)
	parser_readDB.add_argument('-s', '--savename', help='Name of output file', default='excavate', type=str)
	parser_readDB.add_argument('-t', '--threshold', default=60, help='Probability threashold for methylation call (%)', type=float)


	# PARSER : EXPORT
	parser_export.add_argument('-i', '--input_file', help='File with readDB data', required = True)
	parser_export.add_argument('-o', '--output_directory', help='Name of output directory', default='out', type=str)
	parser_export.add_argument('-p', default=cwd, help='Path to output', type=str)
	parser_export.add_argument('-s', '--savename', default='excavate')
	parser_export.add_argument('-t', '--threshold', default=60, help='Probability threashold for methylation call (%)', type=float)

	args = parser.parse_args()

	main(args)



