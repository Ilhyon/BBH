#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
import os
import argparse
from pprint import pprint
from Bio.Blast import NCBIXML
	
#----------------------------------------------------------------------#
def ParseFile(fileContent):
	fileParsed = fileContent
	return fileParsed
#----------------------------------------------------------------------#
def parserBold(element):
	fileparse = []
	e1 = r"(?i)(?P<run>g{3,})(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)(.{1,7}?)(?P=run)"
	#~ e1 = r"0\t[0-9]{2}nt,\s>ENSMUST[0-9]{1,15}.{7}\s.\n>"
	if(re.search(e1,element)): # si tracks de G
		return re.findall(e1,element)
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
				idBBH = words[0].rstrip() 
				idBBH = idBBH.split("|")[1]
				sequence = words[1].rstrip() 
				if idBBH not in DicoFile1 :
					DicoFile1[idBBH] = sequence
	return(DicoFile1)
#----------------------------------------------------------------------#
def ImportBBH(filename):
	""" Create a dictonary that contain all the regions G4s and there BBH
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene
	     of both specie as BBH
	    Returns
	    -------
	     DicoBBH : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	"""
	DicoBBH = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l and l != lines[0]:
				words=l.split('\t')
				idG4Sp1 = words[0].rstrip()  # the first word correspond to the id of the G4 from the first specie
				idG4Sp2 = words[1].rstrip()  # the second word correspond to the G4 from the second specie
				DicoBBH[idG4Sp1]= idG4Sp2  # we add the two id of the G4 in a dictionary
	return(DicoBBH)
#----------------------------------------------------------------------#
def SearchCanoG4s(dicoBBH, dicoSp1, dicoSp2):
	""" 
		Parameters
	    ----------
	    Returns
	    -------
	"""
	for idG4Sp1 in dicoBBH:
		idG4Sp2 = dicoBBH[idG4Sp1]
			
		if parserBold(dicoSp1[idG4Sp1]) :
			if parserBold(dicoSp2[idG4Sp2]) :
				if len(dicoSp1[idG4Sp1]) >= 15 or len(dicoSp2[idG4Sp2]) >= 15 :
					print "canonical\t"+str(dicoSp1[idG4Sp1].count("G"))+"\t"+str(dicoSp2[idG4Sp2].count("G"))
		if parserBold(dicoSp1[idG4Sp1]) :
			if not parserBold(dicoSp2[idG4Sp2]) :
				if len(dicoSp1[idG4Sp1]) >= 15 or len(dicoSp2[idG4Sp2]) >= 15 :
					print "canonical-non_canonical\t"+str(dicoSp1[idG4Sp1].count("G"))+"\t"+str(dicoSp2[idG4Sp2].count("G"))
		if not parserBold(dicoSp1[idG4Sp1]) :
			if parserBold(dicoSp2[idG4Sp2]) :
				if len(dicoSp1[idG4Sp1]) >= 15 or len(dicoSp2[idG4Sp2]) >= 15 :
					print "canonical-non_canonical\t"+str(dicoSp1[idG4Sp1].count("G"))+"\t"+str(dicoSp2[idG4Sp2].count("G"))
		if not parserBold(dicoSp1[idG4Sp1]) :
			if not parserBold(dicoSp2[idG4Sp2]) :
				if len(dicoSp1[idG4Sp1]) >= 15 or len(dicoSp2[idG4Sp2]) >= 15 :
					print "non_anonical\t"+str(dicoSp1[idG4Sp1].count("G"))+"\t"+str(dicoSp2[idG4Sp2].count("G"))
#----------------------------------------------------------------------#
def SearchCanoHit(dicoBBH, dicoSp1, dicoSp2, pathRes):
	""" 
		Parameters
	    ----------
	    Returns
	    -------
	"""
	bestHit = {}
	for idG4Sp1 in dicoBBH:
		idG4Sp2 = dicoBBH[idG4Sp1]
		fileSp1 = str(pathRes+"resultatHvsM/"+idG4Sp1+".blastresult.xml")
		fileSp2 = str(pathRes+"resultatMvsH/"+idG4Sp2+".blastresult.xml")
		
		xml_file = open(fileSp1) # opening of the file
		blast_out = NCBIXML.parse(xml_file) # parsing of the file
		for record in blast_out: # browse all records of the blast file
			bestHit[idG4Sp1] = ['', 0] # initialisation of the list for the new query
			for alignment in record.alignments: # browse all hit
				targetID = alignment.hit_def # retrieval of the target name
				if len(targetID.split("|"))>1:
					targetID = targetID.split("|")[1]
				maxscore = 0
				for hsp in alignment.hsps:
					# filter by e-value
					if hsp.expect < 0.01 and hsp.align_length > 6:
						#~ print hsp.score
						maxscore = max(maxscore, hsp.score)
				if maxscore > bestHit[idG4Sp1][1]:
					bestHit[idG4Sp1]=[targetID, maxscore, hsp.query, hsp.sbjct]
					#~ print hsp.query
					#~ print hsp.match
					#~ print hsp.sbjct
					#~ print "--------------------------------------------"
		xml_file.close()

		if parserBold(bestHit[idG4Sp1][2]) :
			if parserBold(bestHit[idG4Sp1][3]) :
				if len(bestHit[idG4Sp1][2]) >= 15 or len(bestHit[idG4Sp1][3]) >= 15 :
					print "canonical\t"+str(bestHit[idG4Sp1][2].count("G"))+"\t"+str(bestHit[idG4Sp1][3].count("G"))
		if parserBold(bestHit[idG4Sp1][2]) :
			if not parserBold(bestHit[idG4Sp1][3]) :
				if len(bestHit[idG4Sp1][2]) >= 15 or len(bestHit[idG4Sp1][3]) >= 15 :
					print "canonical-non_canonical\t"+str(bestHit[idG4Sp1][2].count("G"))+"\t"+str(bestHit[idG4Sp1][3].count("G"))
		if not parserBold(bestHit[idG4Sp1][2]) :
			if parserBold(bestHit[idG4Sp1][3]) :
				if len(bestHit[idG4Sp1][2]) >= 15 or len(bestHit[idG4Sp1][3]) >= 15 :
					print "canonical-non_canonical\t"+str(bestHit[idG4Sp1][2].count("G"))+"\t"+str(bestHit[idG4Sp1][3].count("G"))
		if not parserBold(bestHit[idG4Sp1][2]) :
			if not parserBold(bestHit[idG4Sp1][3]) :
				if len(bestHit[idG4Sp1][2]) >= 15 or len(bestHit[idG4Sp1][3]) >= 15 :
					print "non_anonical\t"+str(bestHit[idG4Sp1][2].count("G"))+"\t"+str(bestHit[idG4Sp1][3].count("G"))
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'Canonique')
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/Blast/')
	parser.add_argument ('-sp1', '--specie1', default = 'HS')
	parser.add_argument ('-sp2', '--specie2', default = 'MM')
	parser.add_argument ('-pB', '--pathBlast', default = '/media/anais/08c4bb0b-43b3-4183-88a4-9ca53d7ce1e8/home/anais/Documents/Data/Blast/')
	parser.add_argument ('-pQ', '--pathQuery', default = '/media/anais/08c4bb0b-43b3-4183-88a4-9ca53d7ce1e8/home/anais/Documents/Data/Blast/')
	parser.add_argument ('-o', '--outPut', default = '/media/anais/08c4bb0b-43b3-4183-88a4-9ca53d7ce1e8/home/anais/Documents/Data/Blast/')
	parser.add_argument ('-s', '--score', default = True)
	parser.add_argument ('-q', '--query', default = False)
	parser.add_argument ('-m', '--match', default = False)
	parser.add_argument ('-e', '--evalue', default = 0.0001)
	parser.add_argument ('-specie', '--specie', default = 'MM')
	return parser
#----------------------------------------------------------------------#
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	path=arg.path	# directory which contain all the file with BBH and rG4s
	specie1=arg.specie1	# first specie to analyse
	specie2=arg.specie2	# first specie to analyse
	pathOutput=arg.outPut # directory of the output
	score=arg.score # if the user want the score or not in the output
	query=arg.query # if the user want the query or not in the output
	match=arg.match # if the user want the score or not in the output
	evalue=arg.evalue # evalue to filter the result
	
	DicoBBHSp1 = importFile(path+"G4_BBH_HS.txt")
	DicoBBHSp2 = importFile(path+"G4_BBH_MM.txt")
	
	
	DicoBBH  = ImportBBH(path+"BBH.txt")
	
	#~ SearchCanoG4s(DicoBBH, DicoBBHSp1, DicoBBHSp2)
	SearchCanoHit(DicoBBH, DicoBBHSp1, DicoBBHSp2, pathOutput)
#----------------------------------------------------------------------#
main()
