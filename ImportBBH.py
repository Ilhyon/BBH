#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re
from pprint import pprint

#----------------------------------------------------------------------#
def ImportBBH(filename):
	""" Create a dictonary that contain all the regions G4s and there BBH
	Parameters
	    ----------
	     filename : string, name of the file containing all the gene
	     of both specie as BBH
	    Returns
	    -------
	     DicoBBH : {idG4Sp1 : idG4Sp2}
	"""
	DicoBBH = {}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			words=l.split('\t')
			idG4Sp1 = words[0] # the first word correspond to the id of the G4 from the first specie
			idG4Sp2 = words[1] # the second word correspond to the G4 from the second specie
			DicoBBH[idG4Sp1]= idG4Sp2  # we add the two id of the G4 in a dictionary
	return(DicoBBH)
#----------------------------------------------------------------------#
def ProteinToTranscript(elem, DicoGTP):
	""" This function allow us to find if eleme is an id protein and to 
		return the id transcript of this id protein
		Parameters
	    ----------
	    elem : string, correspond to the element we want to find out if
				it's a portein id or a transcript id
	    DicoGTP :	{"Transcript" : {idTranscript : {idGene : idProtein}}
					 "Protein" : {idProtein : idTranscript}}
	    Returns
	    -------
	    transcriptID : transcript id of the element
	"""
	
	if re.search("SP", elem): # the element 1 is a protein
		if elem in DicoGTP["Protein"] :
			transcriptID = DicoGTP["Protein"][elem][0]
			if type(transcriptID) is list :
				print transcriptID
		else:
			transcriptID = "No Transcript ID"
			#~ print elem
	else:
		transcriptID = elem

	return(transcriptID)
#----------------------------------------------------------------------#
def ImportIDGeneTranscriptProteins(filename):
	""" Create a dictionary with for each gene we will have the transcript ID and
	the protein ID, if there is no protein ID in ENSEMBL there will be a label "No protein"
	Parameters
	    ----------
	     filename : string, name of the file that contain for each gene the transcript id and the protein id
	    Returns
	    -------
	     DicoGeneTranscriptProtein : {"Transcript" : {idTranscript : {idGene : idProtein}}
									  "Protein" : {idProtein : idTranscript}}
	"""
	DicoGeneTranscriptProtein = {"Transcript" : {}, "Protein" : {}}
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l : 
				words=l.split('\t')
				idGene = words[0]
				idTranscrit = words[1]
				if words[2] : # there is a protein id 
					idProtein = words[2]
					if idTranscrit not in DicoGeneTranscriptProtein["Transcript"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Transcript"][idTranscrit] = {idGene : [idProtein]}
					elif idTranscrit in DicoGeneTranscriptProtein["Transcript"] and idProtein not in DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene].append(idProtein)
					if idProtein not in DicoGeneTranscriptProtein["Protein"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Protein"][idProtein] = [idTranscrit]
					elif idProtein in DicoGeneTranscriptProtein["Protein"] and idTranscrit not in DicoGeneTranscriptProtein["Protein"][idProtein] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Protein"][idProtein].append(idTranscrit)
				else : # there is no protein id
					idProtein = "No protein"
					if idTranscrit not in DicoGeneTranscriptProtein["Transcript"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Transcript"][idTranscrit] = {idGene : [idProtein]}
					elif idTranscrit in DicoGeneTranscriptProtein["Transcript"] and idProtein not in DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Transcript"][idTranscrit][idGene].append(idProtein)
					if idProtein not in DicoGeneTranscriptProtein["Protein"] : # if the transcript has not been encounter yet
						# we add it to the dictionary with the first protein encounter
						DicoGeneTranscriptProtein["Protein"][idProtein] = [idTranscrit]
					elif idProtein in DicoGeneTranscriptProtein["Protein"] and idTranscrit not in DicoGeneTranscriptProtein["Protein"][idProtein] :
						# if the transcript is already in the dictionary but not the idprot
						DicoGeneTranscriptProtein["Protein"][idProtein].append(idTranscrit)
	return(DicoGeneTranscriptProtein)
#----------------------------------------------------------------------#
def MajOrthology(GeneSp1, TranscriptSp1, GeneSp2, TranscriptSp2, DicoHomology):
	""" When reading the file of orthology of the specie 2 we add new informations
		if their are some
		Parameters
	    ----------
	    GeneSp1 : string, gene id from the specie 1
	    TranscriptSp1 : string, transcript id from the specie 1
	    GeneSp2 : string, gene id from the specie 2
	    TranscriptSp2 : string, transcript id from the specie 2
	    DicoHomology : {"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	    Returns
	    -------
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	"""
	
	# we start with the gene level
	if GeneSp1 in DicoHomology["Gene"]["Orthology"] :
		if GeneSp2 not in DicoHomology["Gene"]["Orthology"][GeneSp1] :
			#the gene from the specie 1 is known but there is a new orthologues form him
			DicoHomology["Gene"]["Orthology"][GeneSp1].append(GeneSp2)
	else :
		#the gene from the specie 1 is not known so we add it
		DicoHomology["Gene"]["Orthology"][GeneSp1] = []
		DicoHomology["Gene"]["Orthology"][GeneSp1].append(GeneSp2)
	
	# then the transcript level
	if TranscriptSp1 in DicoHomology["Transcript"]["Orthology"] :
		if TranscriptSp2 not in DicoHomology["Transcript"]["Orthology"][TranscriptSp1] :
			#the transcript from the specie 1 is known but there is a new orthologues form him
			DicoHomology["Transcript"]["Orthology"][TranscriptSp1].append(TranscriptSp2)
	else :
		#the transcript from the specie 1 is not known so we add it
		DicoHomology["Transcript"]["Orthology"][TranscriptSp1] = []
		DicoHomology["Transcript"]["Orthology"][TranscriptSp1].append(TranscriptSp2)
	return(DicoHomology)
#----------------------------------------------------------------------#
def AddGeneHomology(elem1, elem2, DicoHomologyLevel, homologyType):
	""" Add to the dictionary of homology (DicoHomology) the relation of orthology/paralogy,
		if those relations are not in the dictionary or even if they are 
		already in the dictionary, new information about the Gene/transcript are added
	Parameters
	    ----------
	     elem1 : string, correspond to a gene of the specie 1
	     elem2 : string, correspond to a gene of the specie 1 (for paralogy)
				 or to a gene of the specie 2 (for orthology)
	     DicoHomologyLevel : the subdictionary of DicoHomology with only 
							the level of gene or transcript
	     homologyType : string, can be "Paralogues" or "Orthologues"
	    Returns
	    -------
	     DicoHomologyLevel[homologyType][elem1] : {homologyType: {elem1 : [elem2]}}
	"""
	if elem1 not in DicoHomologyLevel: # if the gene from the specie 1 is not
		# in the dictionary, we add it and we create the list of homologs genes
		DicoHomologyLevel[homologyType][elem1] = []
		DicoHomologyLevel[homologyType][elem1].append(elem2) # add of the "first" homologue gene
	elif elem1 in DicoHomologyLevel[homologyType] and elem2 not in DicoHomologyLevel[homologyType][elem1] :
	 # the the gene from specie one is already in the dictionary but he 
	 # got a new homologue gene 
	 # (it is possible that the same pairs of gene homologue appear in the file, that's why this step is needed)
		DicoHomologyLevel[homologyType][elem1].append(elem2)
	else :
		print "Patate 6 : problem"
	return(DicoHomologyLevel)
#----------------------------------------------------------------------#
def ImportOrthology1(filename, DicoGTPSp1, DicoGTPSp2, DicoHomology):
	""" Import into a dictionnary, the orthology links between the genes 
		and transcript from the specie 1 vs the specie 2
		Parameters
	    ----------
	    filename : string, name of the file containg the information of orthology of the Specie 1 vs the Specie 2
	    DicoGTPSp1 : DicoGTP : {"Transcript" : {idTranscript : {idGene : idProtein}}
								"Protein" : {idProtein : idTranscript}}
	    DicoGTPSp2 : DicoGTP : {"Transcript" : {idTranscript : {idGene : idProtein}}
								"Protein" : {idProtein : idTranscript}}
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	    Returns
	    -------
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	"""
	homologyType = "Orthology"
	
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				
				GeneSp1 = words[0]
				TranscriptSp1 = words[1]
				GeneSp2 = words[2]
				TranscriptSp2 = words[3]
				
				TranscriptSp1 = ProteinToTranscript(TranscriptSp1, DicoGTPSp1)
				TranscriptSp2 =  ProteinToTranscript(TranscriptSp2, DicoGTPSp2)
				
				DicoHomology["Gene"].update(AddGeneHomology(GeneSp1, GeneSp2, DicoHomology["Gene"], homologyType))
				DicoHomology["Transcript"].update(AddGeneHomology(TranscriptSp1, TranscriptSp2, DicoHomology["Transcript"], homologyType))
	return(DicoHomology)
#----------------------------------------------------------------------#
def ImportOrthology2(filename, DicoGTPSp1, DicoGTPSp2, DicoHomology):
	""" Import into a dictionnary, the orthology links between the genes 
		and transcript from the specie 2 vs the specie 1
		Parameters
	    ----------
	    filename : string, name of the file containg the information of orthology of the Specie 2 vs the Specie 1
	    DicoGTPSp1 : {"Transcript" : {idTranscript : {idGene : idProtein}}
					  "Protein" : {idProtein : idTranscript}}
	    DicoGTPSp2 : {"Transcript" : {idTranscript : {idGene : idProtein}}
					  "Protein" : {idProtein : idTranscript}}
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	    Returns
	    -------
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	"""
	homologyType = "Orthology"
	
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				
				GeneSp2 = words[0]
				TranscriptSp2 = words[1]
				GeneSp1 = words[2]
				TranscriptSp1 = words[3]
				
				TranscriptSp1 = ProteinToTranscript(TranscriptSp1, DicoGTPSp1)
				TranscriptSp2 =  ProteinToTranscript(TranscriptSp2, DicoGTPSp2)
				
				DicoHomology.update(MajOrthology(GeneSp1, TranscriptSp1, GeneSp2, TranscriptSp2, DicoHomology))
	return(DicoHomology)
#----------------------------------------------------------------------#
def ImportParalogy(filename, DicoGTP, specie, DicoHomology):
	""" Create two dictonary that contain Homologie between the two specie,
		one with the homology between genes of the two specie and the other one
		between the transcripts of the two specie
		Parameters
	    ----------
	    filename : string, name of the file containing all the gene and transcript
	    that are homologue between two specie
	    DicoGTP : {"Transcript" : {idTranscript : {idGene : idProtein}}
				   "Protein" : {idProtein : idTranscript}}
	    specie : string, name of the specie (HS or MM for exemple)
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	    Returns
	    -------
	    DicoHomology :	{"Gene" : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								   "Paralogie Sp1" : {idGeneSp1 : [idGeneSp1]}
								   "Paralogie Sp2" : {idGeneSp2 : [idGeneSp2]}}
						 "Transcript" : {"Orthology" : {idTranscriptSp1 : [idTranscriptSp2]}
										 "Paralogie Sp1" : {idTranscriptSp1 : [idTranscriptSp1]}
										 "Paralogie Sp2" : {idTranscriptSp2 : [idTranscriptSp2]}}}
	"""
	paralogy = "Paralogy " + str(specie)
	
	with open(filename) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				Gene = words[0]
				Transcript = words[1]
				ParaGene = words[2]
				ParaTranscript = words[3]
				
				Transcript = ProteinToTranscript(Transcript, DicoGTP)
				ParaTranscript =  ProteinToTranscript(ParaTranscript, DicoGTP)
				
				DicoHomology["Gene"].update(AddGeneHomology(Gene, ParaGene, DicoHomology["Gene"], paralogy))
				DicoHomology["Transcript"].update(AddGeneHomology(Transcript, ParaTranscript, DicoHomology["Transcript"], paralogy))
				
	return(DicoHomology)
#----------------------------------------------------------------------#
def importInfoG4Transcript(filenameG4InTranscript, DicoInfoG4):
	""" Create a dictonary for each G4 with his localisations and biotype at the 
		transcript level
	Parameters
	    ----------
	    filenameG4InTranscript : string, name of the file containing all 
	    the G4 In Transcript with their localisations and biotype (transcript level)
	    DicoInfoG4 : {"Gene": {}, "Transcript" : {}}
	    Returns
	    -------
	    DicoInfoG4["Transcript"] :	{idG4 : {idTranscript : {"Localisation" : [localisation]
															 "Biotype" : [Biotype]}}}
	"""
	# from the file G4InTranscript we will extract the data of transcript : localisation and biotype
	with open(filenameG4InTranscript) as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				ids = words[0]
				ids = ids.split('|')
				idTranscript = ids[0]
				G4 = ids[1]
				localisation = words[5]
				biotype = words[6]
				if G4 not in DicoInfoG4["Transcript"]: # the g4 is not yet encounter
					DicoInfoG4["Transcript"][G4] = {}
					DicoInfoG4["Transcript"][G4][idTranscript] = {"Localisations" : []}
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
				elif G4 in DicoInfoG4["Transcript"] and idTranscript not in DicoInfoG4["Transcript"][G4] :
					# the G4 is encounter but not in this transcript
					DicoInfoG4["Transcript"][G4][idTranscript] = {"Localisations" : []}
					DicoInfoG4["Transcript"][G4][idTranscript]["Localisations"] = [localisation]
					DicoInfoG4["Transcript"][G4][idTranscript]["Biotype"] = [biotype]
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def importInfoG4Gene(DicoGTP, DicoInfoG4):
	""" Create a dictonary for each G4 with his localisations and biotype at the 
		gene level
	Parameters
	    ----------
	    DicoGTP : {"Transcript" : {idTranscript : {idGene : idProtein}}
				   "Protein" : {idProtein : idTranscript}}
		DicoInfoG4 : {"Gene": {}, "Transcript" : {}}
	    Returns
	    -------
	    DicoInfoG4 : {idG4 : {idGene : {"Localisation" : [localisation]
										"Biotype" : [Biotype]}}}
	"""
	# with the new informations and the dicoGTP we can get the missing ones :
	# for genes level the localisation and biotypes
	for G4 in DicoInfoG4["Transcript"]: # browse all G4
		G4PerTranscript = DicoInfoG4["Transcript"][G4] # get the sub dictionary of transcript where this G4 is
		for Transcrit in G4PerTranscript : # browse all the transcript of this G4
			if Transcrit in DicoGTP["Transcript"] : 
				gene = DicoGTP["Transcript"][Transcrit].keys() # get only the name of the transcript's gene
			else:
				gene = [""]
			gene = gene[0]
			#~ print gene
			if G4 not in DicoInfoG4["Gene"] : # if the G4 is not already in the dictionary
				# we add it with the first gene encounter and his relative informations
				DicoInfoG4["Gene"][G4] = {}
				DicoInfoG4["Gene"][G4][gene] = {}
				DicoInfoG4["Gene"][G4][gene] = G4PerTranscript[Transcrit]
			elif G4 in DicoInfoG4["Gene"] and gene not in DicoInfoG4["Gene"][G4] :
				# if the G4 is already in the dictionary but not the gene
				DicoInfoG4["Gene"][G4][gene] = {}
				DicoInfoG4["Gene"][G4][gene] = G4PerTranscript[Transcrit]
			elif G4 in DicoInfoG4["Gene"] and gene in DicoInfoG4["Gene"][G4] :
				# if the G4 is already in the dictionary but and the gene also
				# but for some other transcript we may have other informations 
				# we add it after verifying if we have new informations
				localisationGene = DicoInfoG4["Gene"][G4][gene]["Localisations"]
				biotypeGene = DicoInfoG4["Gene"][G4][gene]["Biotype"]
				localisationTranscript = G4PerTranscript[Transcrit]["Localisations"][0]
				biotypeTranscript = G4PerTranscript[Transcrit]["Biotype"][0]
				# 2 if because there can those two condition can be verified
				if localisationTranscript not in localisationGene :
					DicoInfoG4["Gene"][G4][gene]["Localisations"].append(localisationTranscript)
				if biotypeTranscript not in biotypeGene : 
					DicoInfoG4["Gene"][G4][gene]["Biotype"].append(biotypeTranscript)
	return(DicoInfoG4)
#----------------------------------------------------------------------#
def importInfoG4(filenameG4InTranscript, DicoGTP):
	""" Create a dictonary for each G4 with his localisation
	Parameters
	    ----------
	    filenameG4InTranscript : string, name of the file containing all 
	    the G4 In Transcript with their localisations and biotype (transcript level)
	    DicoGTP : keys (string) -> rG4s from specie 1; value (string) -> rG4s from specie 2
	    Returns
	    -------
	    DicoInfoG4 :	{"Transcript" : {idG4 : {idTranscript : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idG4 : {idGene : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
	"""
	DicoInfoG4 = {"Gene": {}, "Transcript" : {}}
	
	DicoInfoG4.update(importInfoG4Transcript(filenameG4InTranscript, DicoInfoG4))
	DicoInfoG4.update(importInfoG4Gene(DicoGTP, DicoInfoG4))
	
	return(DicoInfoG4)
