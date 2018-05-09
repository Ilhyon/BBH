#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
from pprint import pprint

#----------------------------------------------------------------------#
def importFile(filename):
	""" 
		Parameters
	    ----------
	    Returns
	    -------
	"""
	DicoFile1 = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: # browse all line
			if l and l != lines[0]:
				words=l.split('\t')
				idBBH = words[0]
				coupleHomologue = words[1]
				if words[2] :
					localisation = words[2]
				else :
					localisation = ""
				if words[3] :
					biotype = words[3]
				else :
					biotype = words[3]
				if idBBH not in DicoFile1 :
					DicoFile1[idBBH] = {coupleHomologue : {}}
					DicoFile1[idBBH][coupleHomologue].update({"Localisations" : localisation, "Biotypes" : biotype})
				elif idBBH in DicoFile1 and coupleHomologue not in DicoFile1[idBBH] :
					DicoFile1[idBBH] = {coupleHomologue : {}}
					DicoFile1[idBBH][coupleHomologue].update({"Localisations" : localisation, "Biotypes" : biotype})
				#~ else :
					#~ print idBBH
					#~ print l
					#~ pprint(DicoFile1[idBBH])
					#~ print "---------------------------------------------"
	return(DicoFile1)
#----------------------------------------------------------------------#
def GetInfo(dico1, dico2, dico3, dico4):
	""" 
		Parameters
	    ----------
	    Returns
	    -------
	"""
	
	for BBH in dico4 :
		if BBH not in dico2 :
			if BBH not in dico1:
				if BBH not in dico4 :
					print "BBH in dico1 but not in dico2"
			
		#~ if BBH in dico2:
			#~ print "BBH in dico1 and in dico2"
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'InfoBBH')
	parser.add_argument ('-p', '--path', default = '/home/local/USHERBROOKE/vana2406/Documents/Data/BBH/')
	parser.add_argument ('-sp1', '--specie1', default = 'HS')
	parser.add_argument ('-sp2', '--specie2', default = 'MM')
	return parser
#----------------------------------------------------------------------#
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the file with BBH and rG4s
	specie1=arg.specie1	# first specie to analyse
	specie2=arg.specie2	# first specie to analyse
	
	DicoOrthology = importFile(path+"InfoOrthology_Transcript_BBHrG4.txt")
	DicoHomology = importFile(path+"InfoHomology_Transcript_BBHrG4.txt")
	DicoNoOrth = importFile(path+"InfoNo_orthology_Transcript_BBHrG4.txt")
	dicoPara = importFile(path+"InfoParalogy_Transcript_BBHrG4.txt")

	GetInfo(DicoOrthology, DicoHomology, DicoNoOrth, dicoPara)
	#~ GetInfo(DicoHomology, DicoOrthology)

main()
