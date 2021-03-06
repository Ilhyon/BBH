#!/usr/bin/env python
# -*- coding: utf-8 -*-:

from Bio.Blast import NCBIXML
import os
import argparse
			
#----------------------------------------------------------------------#
def ParseFile(fileContent):
	fileParsed = fileContent
	return fileParsed
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'parse_blast_xml')
	parser.add_argument ('-pB', '--pathBlast', default = '/home/anais/Documents/Data/Blast/Data/')
	parser.add_argument ('-pQ', '--pathQuery', default = '/home/anais/Documents/Data/Blast/Data/')
	parser.add_argument ('-o', '--outPut', default = '/home/anais/Documents/Data/Blast/Data/')
	parser.add_argument ('-s', '--score', default = True)
	parser.add_argument ('-q', '--query', default = False)
	parser.add_argument ('-m', '--match', default = False)
	parser.add_argument ('-e', '--evalue', default = 0.0001)
	parser.add_argument ('-specie', '--specie', default = 'MM')
	return parser
#----------------------------------------------------------------------#
def GetBestHit(path):
	""" Return a dictionary containing for each query file the best score and the id of the arget sequence
	    Parameters
	    ----------
	    path : str, path of the directory containing all the blast file
		
	    Returns
	    -------
	    bestHit : dictionary of list, 
	    exp --> {Query_1 : [Targetx, score]}
	"""
	filelist = sorted(os.listdir(path))# recuperation of all file name in the directory
	bestHit = {} # initialisation of the dictionary 
	cpt = 0
	bool10 = False
	bool20 = False
	bool30 = False
	bool40 = False
	bool50 = False
	bool60 = False
	bool70 = False
	bool80 = False
	bool90 = False
	bool100 = False
	for blastFile in filelist: # browse all blastfile output
		cpt += 1 
		#~ print dicoNameQuery
		namequery = blastFile.split('.')[0]
		if len(namequery.split("|"))>1:
			namequery = namequery.split("|")[1]
			#print namequery.split("|")
		blastFile = str(path+blastFile)
		xml_file = open(blastFile) # opening of the file
		blast_out = NCBIXML.parse(xml_file) # parsing of the file
		for record in blast_out: # browse all records of the blast file
			bestHit[namequery] = ['', 0] # initialisation of the list for the new query
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
				if maxscore > bestHit[namequery][1]:
					bestHit[namequery]=[targetID, maxscore]
		xml_file.close()
		stat = round((float(cpt)/float(len(filelist)))*100,1)
		#if "E" in namequery :
			#print namequery
		if stat == 10 and bool10 == False:
			bool10 = True
			print "-10%"
		elif stat == 20 and bool20 == False:
			bool20 = True
			print "--20%"
		elif stat == 30 and bool30 == False:
			bool30 = True
			print "---30%"
		elif stat == 40 and bool40 == False:
			bool40 = True
			print "----40%"
		elif stat == 50 and bool50 == False:
			bool50 = True
			print "-----50%"
		elif stat == 60 and bool60 == False:
			bool60 = True
			print "------60%"
		elif stat == 70 and bool70 == False:
			bool70 = True
			print "-------70%"
		elif stat == 80 and bool80 == False:
			bool80 = True
			print "--------80%"
		elif stat == 90 and bool90 == False:
			bool90 = True
			print "---------90%"
		elif stat == 100 and bool100 == False:
			bool100 = True
			print "----------100%"
	return bestHit
#----------------------------------------------------------------------#
def compareBH(dicoHuman, dicoMouse):
	BBH = []
	cpt = 0
	bool10 = False
	bool20 = False
	bool30 = False
	bool40 = False
	bool50 = False
	bool60 = False
	bool70 = False
	bool80 = False
	bool90 = False
	bool100 = False
	query_unknown = 0
	for queryHuman in dicoHuman:
		cpt += 1 
		TargetMouse = dicoHuman[queryHuman][0]
		if TargetMouse != '':
			if dicoMouse[TargetMouse] :
				if dicoMouse[TargetMouse][0] == queryHuman:
					# Bingo ! ~^w^~
					BBH.append(queryHuman+"\t"+TargetMouse.split("|")[0])
			else:
				query_unknown += 1
		stat = round((float(cpt)/float(len(dicoHuman)))*100,1)
		if stat == 10 and bool10 == False:
			bool10 = True
			print "-10%"
		elif stat == 20 and bool20 == False:
			bool20 = True
			print "--20%"
		elif stat == 30 and bool30 == False:
			bool30 = True
			print "---30%"
		elif stat == 40 and bool40 == False:
			bool40 = True
			print "----40%"
		elif stat == 50 and bool50 == False:
			bool50 = True
			print "-----50%"
		elif stat == 60 and bool60 == False:
			bool60 = True
			print "------60%"
		elif stat == 70 and bool70 == False:
			bool70 = True
			print "-------70%"
		elif stat == 80 and bool80 == False:
			bool80 = True
			print "--------80%"
		elif stat == 90 and bool90 == False:
			bool90 = True
			print "---------90%"
		elif stat == 100 and bool100 == False:
			bool100 = True
			print "----------100%"
	#print query_unknown
	return BBH
#----------------------------------------------------------------------#
def main():
	parser = build_arg_parser()
	arg = parser.parse_args()
	pathBlast=arg.pathBlast	# directory which contain the blast file to parse
	pathOutput=arg.outPut # directory of the output
	score=arg.score # if the user want the score or not in the output
	query=arg.query # if the user want the query or not in the output
	match=arg.match # if the user want the score or not in the output
	evalue=arg.evalue # evalue to filter the result
	
	pathBlastH = str(pathBlast+"resultHvsM/")
	print "Dico BH Human :"
	dicoBestHitHumanvsMouse = GetBestHit(pathBlastH)
	#~ for i in dicoBestHitHumanvsMouse:
		#~ mouse = dicoBestHitHumanvsMouse[i][0].split("|")[0]
		#~ print str(i +"\t"+ mouse)
	#~ print "\n"
	
	pathBlastM = str(pathBlast+"resultMvsH/")
	print "Dico BH Mouse :"
	dicoBestHitMousevsHuman = GetBestHit(pathBlastM)
	#print dicoBestHitMousevsHuman
	
	#~ for i in dicoBestHitMousevsHuman:
		#~ mouse = i.split("|")[0]
		#~ print str(mouse +"\t"+ dicoBestHitMousevsHuman[i][0])
	print "Dico BBH :"
	BBH = compareBH(dicoBestHitHumanvsMouse, dicoBestHitMousevsHuman)
	output= open(pathOutput+"/BBH.txt","w") ## file opening
	output.write('\n'.join(str(e) for e in BBH))
	#~ print '\n'.join(str(e) for e in BBH)
#----------------------------------------------------------------------#

main()
