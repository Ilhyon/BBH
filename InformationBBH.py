#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import argparse
import re
import ImportBBH as imp
from pprint import pprint

#----------------------------------------------------------------------#
def GetCommonInformation(LocalisationsSp1, LocalisationsSp2, BiotypesSp1, BiotypesSp2):
	""" Get the commons informations between two list for the Biotypes and the Localisations,
		if there is no commons informations, the list will stay empty
		Parameters
	    ----------
	    LocalisationsSp1 :	list of Localisations from the specie 1
		LocalisationsSp2 :	list of Localisations from the specie 2
		BiotypesSp1 : list of Biotypes from the specie 1
		BiotypesSp2 : list of Biotypes from the specie 2
		
	    Returns
	    -------
	    commonLocalisations : list of commons Localisations of a G4 in a gene/transcript between the 2 specie
	    commonBiotypes : list of commons Biotypes of a G4 in a gene/transcript between the 2 specie
	"""
	commonLocalisations = []
	commonBiotypes = []
	if list(set(LocalisationsSp1).intersection(set(LocalisationsSp2))) :
		commonLocalisations = list(set(LocalisationsSp1).intersection(LocalisationsSp2))
	if list(set(BiotypesSp1).intersection(BiotypesSp2)) :
		commonBiotypes = list(set(BiotypesSp1).intersection(set(BiotypesSp2)))
	return(commonLocalisations, commonBiotypes)
#----------------------------------------------------------------------#
def GetLocalisationsAndBiotype(element1, element2, Dico):
	""" Get the Localisations and the Biotypes in a dictionary of dictionary 
		thanks to 2 elements
		Parameters
	    ----------
	    element1 :	string, element that bust be present in the first level of the dictionary
		element2 :	string, element that bust be present in the second level of the dictionary
		Dico :	{idTranscript/Gene : {idG4 : {"Localisation" : [localisation]
											  "Biotype" : [Biotype]}}}
	    Returns
	    -------
	    Localisations : list of Localisations of a G4 in a gene/transcript
	    Biotypes : list of Biotypes of a G4 in a gene/transcript
	"""
	if element1 in Dico :
		if element2 in Dico[element1] :
			Localisations = Dico[element1][element2]["Localisations"]
			Biotypes = Dico[element1][element2]["Biotype"]
		else : # if a gene/transcript don't have any G4s
			Localisations = []
			Biotypes = []
	else :
		Localisations = []
		Biotypes = []
	return(Localisations, Biotypes)
#----------------------------------------------------------------------#
def AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrtho, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All 
		infromation about the BBH if there is some orthology 
		Parameters
	    ----------
	    DicoInfoBBHLevel : {idBBH : {typeInfo : {couple : {"Common Localisations" : [commonLocalisations]
														   "Common Biotypes" : commonBiotypes}}}
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		coupleOrtho : string, id of the couple like geneSp1|geneSp2 or transcriptSp1|transcriptSp2
		commonLocalisations : list of string, contain the common localisation 
							  between the gene/transcripts that are in the list coupleOrtho.
							  This list could be empty
		commonBiotypes : list of string, contain the common Biotypes 
						 between the gene/transcripts that are in the list coupleOrtho.
						 This list could be empty
		typeInfo : string, here it is "Orthology"
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleOrtho : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it 
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho not in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary but 
		# not the couple there homologues's levelS
		DicoInfoBBHLevel[idBBH][coupleOrtho] = {}
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and coupleOrtho in DicoInfoBBHLevel[idBBH] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new Localisations we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new Biotypes we add them
			DicoInfoBBHLevel[idBBH][coupleOrtho]["Common Biotypes"].append(commonBiotypes)
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def AddInfoNoOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, newList, commonLocalisations, commonBiotypes, typeInfo):
	""" Adds the informations to the final dictionary containing All infromation 
		about the BBH if there is No_orthology links
		Parameters
	    ----------
	    DicoInfoBBHLevel : {idBBH : {typeInfo : {couple : {"Common Localisations" : [commonLocalisations]
														   "Common Biotypes" : commonBiotypes}}}
		idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		newList : list of string, id of the gene that are not orthologues
						but they still have common informations in common
		commonLocalisations : list of string, contain the common localisation 
							  between the gene/transcripts that are in the list coupleOrtho.
							  This list could be empty
		commonBiotypes : list of string, contain the common Biotypes 
						 between the gene/transcripts that are in the list coupleOrtho.
						 This list could be empty
		typeInfo : string, here it is "No_orthology"
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {"Orthology" : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
																	   "Common Biotypes" : commonBiotypes}}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	newList = tuple(newList)
	if idBBH not in DicoInfoBBHLevel :
		# First time the couple of BBH is encounter
		# so we create all information relative to it
		DicoInfoBBHLevel[idBBH] = {typeInfo : {}}
		DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and typeInfo not in DicoInfoBBHLevel[idBBH] : 
		DicoInfoBBHLevel[idBBH].update({typeInfo : {}})
		DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
		DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and newList not in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary but 
		# not the id that are not orthologues, this can be du to many reasons
		# maybe there just one id in plus the list or maybe the list doesn't exist
		keys = DicoInfoBBHLevel[idBBH][typeInfo].keys()
		for existingList in keys :
			if set(existingList).issubset(set(newList)) :
				# the existingList is a subset of newList
				newComLoca = list(set(commonLocalisations) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Localisations"]))
				newComBiot = list(set(commonBiotypes) ^ (DicoInfoBBHLevel[idBBH][typeInfo][existingList]["Common Biotypes"]))
				DicoInfoBBHLevel[idBBH][typeInfo][newList] = dictionary.pop(DicoInfoBBHLevel[idBBH][typeInfo][existingList])
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"].append(newComLoca)
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"].append(newComBiot)
			else :
				# the list of id is new
				DicoInfoBBHLevel[idBBH][typeInfo][newList] = {}
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"] = commonLocalisations
				DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"] = commonBiotypes
	elif idBBH in DicoInfoBBHLevel and newList in DicoInfoBBHLevel[idBBH][typeInfo] :
		# the couple of BBH is already in the dictionary and also the couple of homologue
		# yet we still have a chance to get some new informations so we need to check it
		if list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"]) ^ set(commonLocalisations)) :
			# if there is new Localisations we add them
			DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"].append(list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Localisations"]) ^ set(commonLocalisations)))
		if list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"]) ^ set(commonBiotypes)) :
			# if there is new Biotypes we add them
			DicoInfoBBHLevel[idBBH][idNoOrthology][typeInfo][newList]["Common Biotypes"].append(list(set(DicoInfoBBHLevel[idBBH][typeInfo][newList]["Common Biotypes"]) ^ set(commonBiotypes)))
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel):
	""" Here we test if there is common information between gene that hasn't orthology links
		Parameters
	    ----------
	    nonHomologueslevelSp2 : list of string, list of the gene that are not orthologs to the gene of 
								the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleNoOrth : {"Common Localisations" : [commonLocalisations]
													 "Common Biotypes" : commonBiotypes}}}
	"""
	idNoOrthology = [levelSp1]
	commonLocalisationsNoOrth = []
	commonBiotypesNoOrth = []
	for levelSp2 in nonHomologueslevelSp2 : 
		# retrivial of the Localisations and biotype of the two specie
		# from the gene/transcript that are homologues
		LocalisationsSp1WithoutH, BiotypesSp1WithoutH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
		LocalisationsSp2WithoutH, BiotypesSp2WithoutH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
		commonLocalisationsWithoutH, commonBiotypesWithoutH = GetCommonInformation(LocalisationsSp1WithoutH, LocalisationsSp1WithoutH, BiotypesSp1WithoutH, BiotypesSp1WithoutH)
		if commonLocalisationsWithoutH or commonBiotypesWithoutH :
			idNoOrthology.append(levelSp2)
			commonLocalisationsNoOrth = commonLocalisationsNoOrth + commonLocalisationsWithoutH
			commonBiotypesNoOrth = commonBiotypesNoOrth + commonBiotypesWithoutH
	commonLocalisationsNoOrth = list(set(commonLocalisationsNoOrth))
	commonBiotypesNoOrth =  list(set(commonBiotypesNoOrth))
	DicoInfoBBHLevel = AddInfoNoOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, idNoOrthology, commonLocalisationsNoOrth, commonBiotypesNoOrth, "No_orthology")
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetInfoOfParalogues(listParalogues, idG4s,DicoHomologyLevel, DicoInfoG4Level):
	""" The gene are orthologues, so we will search if they have also paralogues
		Parameters
	    ----------
	    listParalogues : list of string, list that contain all the paralogues of a gene/transcript
	    idG4 : string, id of the G4s that we are trying to analise
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
	    Returns
	    -------
	    paraloguesG4 : {"Localisations" : {Loc1 : [paralogues]}
						"Biotypes" : {Bt1 : [paralogues]}}
	"""

	paraloguesG4 = {"Localisations" : {}, "Biotype" : {}}
	for paralogue in listParalogues :
		if paralogue in DicoInfoG4Level[idG4s] :
			for loc in DicoInfoG4Level[idG4s][paralogue]["Localisations"]:
				if loc not in paraloguesG4["Localisations"] :
					paraloguesG4["Localisations"][loc] = []
					paraloguesG4["Localisations"][loc].append(paralogue)
				elif loc in paraloguesG4["Localisations"] and paralogue not in paraloguesG4["Localisations"][loc] : 
					paraloguesG4["Localisations"][loc].append(paralogue)
			for Bt in DicoInfoG4Level[idG4s][paralogue]["Biotype"] : 
				if Bt not in paraloguesG4["Biotype"] :
					paraloguesG4["Biotype"][Bt] = []
					paraloguesG4["Biotype"][Bt].append(paralogue)
				elif Bt in paraloguesG4["Biotype"] and paralogue not in paraloguesG4["Biotype"][Bt] : 
					paraloguesG4["Biotype"][Bt].append(paralogue)
	return(paraloguesG4)
#----------------------------------------------------------------------#
def GetCommonInfoOfParalogues(paraloguesG4Sp1, paraloguesG4Sp2):
	""" The gene are orthologues, so we will search if they have also paralogues
		Parameters
	    ----------
	    paraloguesG4Sp1/Sp2 : {"Localisations" : {Loc1 : [paralogues]}
							   "Biotypes" : {Bt1 : [paralogues]}}
	    Returns
	    -------
	    commonParalogues : list of pralogues that have common information
	    commonloc : list of common Localisations
	    commonBt : list of common biotype
	"""
	commonloc = []
	commonBt = []
	commonParalogues = []
	for loc in paraloguesG4Sp1["Localisations"] :
		if loc in paraloguesG4Sp2["Localisations"] :
			# the localisation is present in both species's paralogues
			# so we add those information to the common lists
			commonloc.append(loc)
			if len(paraloguesG4Sp2["Localisations"][loc]) > 1 :
				for i in range(len(paraloguesG4Sp2["Localisations"][loc])) :
					#~ print paraloguesG4Sp2["Localisations"][loc][i]
					commonParalogues.append(paraloguesG4Sp2["Localisations"][loc][i])
			if len(paraloguesG4Sp1["Localisations"][loc]) > 1 :
				for i in range(len(paraloguesG4Sp1["Localisations"][loc])) :
					commonParalogues.append(paraloguesG4Sp1["Localisations"][loc][i])
					#~ print paraloguesG4Sp1["Localisations"][loc][i]
	for Bt in paraloguesG4Sp1["Biotype"] :
		if Bt in paraloguesG4Sp2["Biotype"] :
			# the biotype is present in both species's paralogues
			# so we add those information to the common lists
			commonBt.append(Bt)
			if len(paraloguesG4Sp2["Biotype"][Bt]) > 1 :
				for i in range(len(paraloguesG4Sp2["Biotype"][Bt])) :
					commonParalogues.append(paraloguesG4Sp2["Biotype"][Bt][i])
					#~ print paraloguesG4Sp2["Biotype"][Bt][i]
			if len(paraloguesG4Sp1["Biotype"][Bt]) > 1 :
				for i in range(len(paraloguesG4Sp1["Biotype"][Bt])) :
					commonParalogues.append(paraloguesG4Sp1["Biotype"][Bt][i])
					#~ print paraloguesG4Sp1["Biotype"][Bt][i]
	return(commonParalogues, commonloc, commonBt)
#----------------------------------------------------------------------#
def GetInfoParalogy(levelSp1, levelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevelHom, DicoHomologyLevel, specie1, specie2):
	""" The gene are orthologues, so we will search if they have also paralogues
		Parameters
	    ----------
	    levelSp1 : string, transcript/gene from the specie 1
	    levelSp2 : string, transcript/gene from the specie 2
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		
	    Returns
	    -------
	    DicoHomologyLevelHom : {commonParalogues : {"Common Localisations" : [commonLocalisations]
													"Common Biotypes" : [commonBiotypes]}}}
	"""
	if levelSp2 in DicoHomologyLevel["Paralogy "+str(specie2)] :
		paraloguesSp2 = DicoHomologyLevel["Paralogy "+str(specie2)][levelSp2]
		paraloguesSp2.append(levelSp2)
	else :
		paraloguesSp2 = []
	if levelSp1 in DicoHomologyLevel["Paralogy "+str(specie1)] :
		paraloguesSp1 = DicoHomologyLevel["Paralogy "+str(specie1)][levelSp1]
		paraloguesSp1.append(levelSp1)
	else:
		paraloguesSp1 = []
	paraloguesG4Sp1 = GetInfoOfParalogues(paraloguesSp1, idG4Sp1,DicoHomologyLevel, DicoInfoG4Sp1Level)
	paraloguesG4Sp2 = GetInfoOfParalogues(paraloguesSp2, idG4Sp2,DicoHomologyLevel, DicoInfoG4Sp2Level)
	
	commonParalogues, commonloc, commonBt = GetCommonInfoOfParalogues(paraloguesG4Sp1, paraloguesG4Sp2)
	
	if tuple(commonParalogues) not in DicoInfoBBHLevelHom :
		DicoInfoBBHLevelHom[tuple(commonParalogues)] = {"Common Localisations" : [], "Common Biotypes" : []}
		for i in range(len(commonloc)):
			DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Localisations"].append(commonloc[i])
		for i in range(len(commonBt)):
			DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Biotypes"].append(commonBt[i])
	elif tuple(commonParalogues) in DicoInfoBBHLevelHom :
		if set(commonloc) ^ set(DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Localisations"]) :
			if len(list(set(commonloc) ^ set(DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Localisations"]))) > 1 :
				newlist = list(set(commonloc) ^ set(DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Localisations"]))
				for i in range(len(newlist)) :
					DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Localisations"].append(newlist[i])
		if set(commonBt) ^ set(DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Biotypes"]) :
			newlist = list(set(commonBt) ^ set(DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Biotypes"]))
			for i in range(len(newlist)) :
				DicoInfoBBHLevelHom[tuple(commonParalogues)]["Common Biotypes"].append(newlist[i])
	return(DicoInfoBBHLevelHom)
#----------------------------------------------------------------------#
def ParalogyNoHomology(levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevelHom, DicoHomologyLevel, specie1, specie2):
	""" We test if there is paralogy
		Parameters
	    ----------
		levelSp1 : string, transcript/gene from the specie 1
	    levelsSp2 : list of string, transcript/gene from the specie 2
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
	    Returns
	    -------
	    DicoTmp : {commonParalogues : {"Common Localisations" : [commonLocalisations]
									   "Common Biotypes" : [commonBiotypes]}}}
	"""
	if levelSp1 in DicoHomologyLevel["Paralogy "+str(specie1)] :
		for levelSp2 in levelsSp2 :
			if levelSp2 in DicoHomologyLevel["Paralogy "+str(specie2)] :
				DicoTmp = GetInfoParalogy(levelSp1, levelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevelHom, DicoHomologyLevel, specie1, specie2)
			else :
				DicoTmp = {"No Paralogy" : {"Common Localisations" : [], "Common Biotypes" : []}} 
				#~ print "Sp2 : " + str(idG4Sp2) + "\t" + str(levelSp2)
	else :
		DicoTmp = {"No Paralogy" : {"Common Localisations" : [], "Common Biotypes" : []}} 
		#~ print "Sp1 : " + str(idG4Sp1) + "\t" + str(levelSp1)
	return(DicoTmp)
#----------------------------------------------------------------------#
def GetInfoHomology(HlevelSp1, levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, DicoHomologyLevel, specie1, specie2):
	""" Here we test many conditions (with or without homology), then we 
		retrieve the informations relative to them
		Parameters
	    ----------
	    HlevelSp1 : list of string, list of the gene that are homologue to the gene of 
					the specie 1 where there is the G4 BBH we are looking at for
		levelSp1 : string, transcript/gene from the specie 1
		levelsSp2 : list of string, list of transcript/gene from the specie 2
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {idGeneSpecie1 : [idGeneSpecie2]}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
	"""
	idBBH = str(idG4Sp1+"|"+idG4Sp2)
	if len(HlevelSp1) > 1 : # there is orthologues
		if list(set(HlevelSp1).intersection(levelsSp2)) :
			# if there an intersection between the homologues of the gene/transcript of specie 1
			# and with the levelS of the specie 2 where the G4 is
			
			for levelSp2 in list(set(HlevelSp1).intersection(levelsSp2)) :
				# retrivial of the Localisations and biotype of the two specie
				# from the gene/transcript that are homologues
				LocalisationsSp1WithH, BiotypesSp1WithH = GetLocalisationsAndBiotype(idG4Sp1, levelSp1, DicoInfoG4Sp1Level)
				LocalisationsSp2WithH, BiotypesSp2WithH = GetLocalisationsAndBiotype(idG4Sp2, levelSp2, DicoInfoG4Sp2Level)
				commonLocalisationsWithH, commonBiotypesWithH = GetCommonInformation(LocalisationsSp1WithH, LocalisationsSp2WithH, BiotypesSp1WithH, BiotypesSp2WithH)
				coupleOrthologue = str(levelSp1+"|"+levelSp2)
				DicoInfoBBHLevel.update(AddInformationsOrthology(DicoInfoBBHLevel, idG4Sp1, idG4Sp2, coupleOrthologue, commonLocalisationsWithH, commonBiotypesWithH, "Orthology"))
				if idBBH not in DicoInfoBBHLevel[idBBH] :
					DicoInfoBBHLevel[idBBH] = {"Orthology" : {}, "Homology" : {}, "No_orthology" : {}, "Paralogy" : {}}
				DicoInfoBBHLevel[idBBH]["Homology"].update(GetInfoParalogy(levelSp1, levelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel[idBBH]["Homology"], DicoHomologyLevel, specie1, specie2))
				
			if list(set(HlevelSp1) ^ set(levelsSp2)) : 
			# retrivial of the list of G4 that are not homologues (unique in the two lists)
				if idBBH not in DicoInfoBBHLevel :
					DicoInfoBBHLevel[idBBH] = {"Orthology" : {}, "Homology" : {}, "No_orthology" : {}, "Paralogy" : {}}
				nonHomologueslevelSp2 = list(set(HlevelSp1) ^ set(levelsSp2))
				DicoInfoBBHLevel.update(GetInfoNoOrthology(levelSp1, nonHomologueslevelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel))
				DicoInfoBBHLevel[idBBH]["Paralogy"].update(ParalogyNoHomology(levelSp1, nonHomologueslevelSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel[idBBH]["Paralogy"], DicoHomologyLevel, specie1, specie2))
		else :
			# there is no element in comon => no orthologues between the gene/transcript
			# where the G4 BBH are
			if idBBH not in DicoInfoBBHLevel :
				DicoInfoBBHLevel[idBBH] = {"Orthology" : {}, "Homology" : {}, "No_orthology" : {}, "Paralogy" : {}}
			DicoInfoBBHLevel.update(GetInfoNoOrthology(levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel))
			DicoInfoBBHLevel[idBBH]["Paralogy"].update(ParalogyNoHomology(levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel[idBBH]["Paralogy"], DicoHomologyLevel, specie1, specie2))
	else : # there is no orthologues
		if idBBH not in DicoInfoBBHLevel :
			DicoInfoBBHLevel[idBBH] = {"Orthology" : {}, "Homology" : {}, "No_orthology" : {}, "Paralogy" : {}}
		DicoInfoBBHLevel.update(GetInfoNoOrthology(levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel))
		DicoInfoBBHLevel[idBBH]["Paralogy"].update(ParalogyNoHomology(levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel[idBBH]["Paralogy"], DicoHomologyLevel, specie1, specie2))
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level, specie1, specie2):
	""" Create a dictonary for each G4 with his localisation and biotype,
		the informations that are retrieve depend on the level transcript/gene
		Parameters
	    ----------
	    idG4Sp1 : string, id of a G4 from the specie 1 which is a BBH with BBHsp2
	    idG4Sp2 : string, id of a G4 from the specie 2 which is a BBH with BBHsp1
		DicoHomologyLevel : {"Orthology" : {idGeneSp1 : [idGeneSp2]}
							 "Paralogues Sp1" : {idGeneSp1 : [idGeneSp1]}
							 "Paralogues Sp2" : {idGeneSp2 : [idGeneSp2]}}
		DicoInfoG4Sp1Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4Sp2Level : {idTranscript : {idG4 : {"Localisation" : [localisation]
													  "Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
	"""
	# first we retrieve the Dico of informations from a level (transcript/gene)
	# where the G4 of each specie is
	if idG4Sp1 in DicoInfoG4Sp1Level : # retrievial of informations for the specie 1
		levelsSp1 =  DicoInfoG4Sp1Level[idG4Sp1]
		
	if idG4Sp2 in DicoInfoG4Sp2Level : # retrievial of informations for the specie 2
		levelsSp2 =  DicoInfoG4Sp2Level[idG4Sp2]
			
	# Then we browse all genes where the G4 from the specie 1 is to get
	# all the informations we want (Paralogy  Sp1/Sp2, Orthology or not)
	for levelSp1 in levelsSp1 : # browse all genes where the G4 from the specie 1 is
		DicoOrthologyLevel = DicoHomologyLevel["Orthology"]
		
		# Orthology
		if levelSp1 in DicoOrthologyLevel : # if this gene have some orthologues
			OrtLevelSp1 = DicoOrthologyLevel[levelSp1] # retrevial of the 
			# gene/transcript from the specie 2 that are orthologues to levelSp1
		else : # there is no orthologues for this gene
			OrtLevelSp1 = []
		DicoInfoBBHLevel.update(GetInfoHomology(OrtLevelSp1, levelSp1, levelsSp2, idG4Sp1, idG4Sp2, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, DicoHomologyLevel, specie1, specie2))
		
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def GetDicoInfo(DicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBHLevel, level, specie1, specie2):
	""" Create a dictonary for each G4 with his localisation and biotype
		the informations that are retrieve depend on the level transcript/gene
	Parameters
	    ----------
	     DicoBBH : dictionnary of BBH, see ImportBBH doc
	     DicoHomology : {"Gene": {"Orthology" : {idGeneSp1 : [idGeneSp2]}
								  "Paralogues Sp1" : {idGeneSp1 : [idGeneSp1]}
								  "Paralogues Sp2" : {idGeneSp2 : [idGeneSp2]}}
						"Transcript" : {"Orthology" : {idTranscririptSp1 : [idTranscriptSp2]}}
										"Paralogues Sp1" : {{idTranscririptSp1 : [idTranscriptSp1]}}
										"Paralogues Sp2" : {{idTranscririptSp2 : [idTranscriptSp2]}} 
		 DicoInfoG4Sp1 :{"Transcript" : {idTranscript : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idGene : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
		 DicoInfoG4Sp2 :{"Transcript" : {idTranscript : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
						 "Gene" : {idGene : {idG4 : {"Localisation" : [localisation]
																"Biotype" : [Biotype]}}}
		DicoInfoG4 : empty dictionnary for the part of gene level, se FindInfoHomology
		level : string, to choose if we want the gene level or the transcript level
		
	    Returns
	    -------
	    DicoInfoBBHLevel : {idBBH : {coupleHomologue : {"Common Localisations" : [commonLocalisations]
														"Common Biotypes" : commonBiotypes}}}
							Yet, the coupleHomologue must be equal to "No homologues"
	"""
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
	DicoInfoG4Sp1Level = DicoInfoG4Sp1[level]
	DicoInfoG4Sp2Level = DicoInfoG4Sp2[level]
	DicoHomologyLevel = DicoHomology[level]
	cpt = 0
	lentot = len(DicoBBH)
	
	for idG4Sp1 in DicoBBH : # browse all BBH couple
		idG4Sp2 = DicoBBH[idG4Sp1]
		
		DicoInfoBBHLevel = GetDicoInfoBBH(idG4Sp1, idG4Sp2, DicoHomologyLevel, DicoInfoG4Sp1Level, DicoInfoG4Sp2Level, DicoInfoBBHLevel, level, specie1, specie2)
		
		cpt += 1
		
		stat = round((float(cpt)/float(len(DicoBBH)))*100,1)
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
	
	return(DicoInfoBBHLevel)
#----------------------------------------------------------------------#
def WriteFile(DicoInfoBBH, path, level, homologyType):
	output = open(path+str(homologyType)+'_'+level+'_BBHrG4.txt',"w") # file opening for reading
	output.write("BBH\tcouple Homologue\tLocalisations\tBiotypes")
	for BBH in DicoInfoBBH[level] :
		if homologyType in DicoInfoBBH[level][BBH] :
			for coupleHomology in DicoInfoBBH[level][BBH][homologyType] :
				if DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Localisations"] :
					if DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Biotypes"] :
						output.write(str(BBH)+"\t"+str(coupleHomology)+"\t"+"|".join(DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Localisations"])+"\t"+"|".join(DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Biotypes"])+"\n")
					else :
						output.write(str(BBH)+"\t"+str(coupleHomology)+"\t"+"|".join(DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Localisations"])+"\t\n")
				else :
					if DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Biotypes"] :
						output.write(str(BBH)+"\t"+str(coupleHomology)+"\t\t"+"|".join(DicoInfoBBH[level][BBH][homologyType][coupleHomology]["Common Biotypes"])+"\n")
					else :
						output.write(str(BBH)+"\t"+str(coupleHomology)+"\t\t\n")
	output.close()
#----------------------------------------------------------------------#
def importData(path, specie1, specie2):
	""" import of all the data we need 
	Parameters
	    ----------
	     filename : strin, name of the path for all file
	    Returns
	    -------
	     	DicoBBH
	     	DicoHomology
	     	DicoGeneTranscriptProteinSpecie2
	     	DicoInfoG4Specie1
	     	DicoInfoG4Specie2

	"""
	filenameBBH = path+'BBH.txt'
	filenameGeneTranscriptProteinSpecie1 = path+specie1+"_All_GTP.txt"
	filenameGeneTranscriptProteinSpecie2 = path+specie2+"_All_GTP.txt"
	filenameParalogySp1 = path+specie1+"_paralogy.txt"
	filenameParalogySp2 = path+specie2+"_paralogy.txt"
	filenameOrthologySp1 = path+specie1+"_orthology_"+specie2+".txt"
	filenameOrthologySp2 = path+specie2+"_orthology_"+specie1+".txt"
	filenameInfoSpecie1 = path+specie1+"_All_G4InTranscript.txt"
	filenameInfoSpecie2 = path+specie2+"_All_G4InTranscript.txt"
	DicoHomology = {"Gene" : {"Paralogy "+str(specie1) : {}, "Paralogy "+str(specie2) : {}, "Orthology" : {}}, "Transcript" : {"Paralogy "+str(specie1) : {}, "Paralogy "+str(specie2) : {}, "Orthology" : {}}}
	
	DicoBBH = imp.ImportBBH(filenameBBH)
	print "Import of BBH -> Done"
	DicoGTPSp1 = imp.ImportIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie1)
	print "Import of GTP from specie 1 -> Done"
	DicoGTPSp2 = imp.ImportIDGeneTranscriptProteins(filenameGeneTranscriptProteinSpecie2)
	print "Import of GTP from specie 2 -> Done"
	DicoInfoG4Specie1 = imp.importInfoG4(filenameInfoSpecie1, DicoGTPSp1)
	print "Import of informations from specie 1 -> Done"
	DicoInfoG4Specie2 = imp.importInfoG4(filenameInfoSpecie2, DicoGTPSp2)
	print "Import of informations from specie 2 -> Done"
	
	DicoHomology.update(imp.ImportParalogy(filenameParalogySp1, DicoGTPSp1, specie1, DicoHomology))
	print "Import of Paralogy for specie 1 -> Done"
	DicoHomology.update(imp.ImportParalogy(filenameParalogySp2, DicoGTPSp2, specie2, DicoHomology))
	print "Import of Paralogy for specie 2 -> Done"
	DicoHomology.update(imp.ImportOrthology1(filenameOrthologySp1, DicoGTPSp1, DicoGTPSp2, DicoHomology))
	print "Import of Orthology links from Sp1 to Sp 2 -> Done"
	DicoHomology.update(imp.ImportOrthology2(filenameOrthologySp2, DicoGTPSp1, DicoGTPSp2, DicoHomology))
	print "Update of Orthology -> Done"

	
	#~ pprint(DicoGeneTranscriptProteinSpecie2)
	#~ pprint(DicoHomology)
	#~ pprint(DicoInfoG4Specie2)
	return(DicoBBH, DicoHomology, DicoInfoG4Specie1, DicoInfoG4Specie2)
#----------------------------------------------------------------------#
def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'InfoBBH')
	parser.add_argument ('-p', '--path', default = '/home/anais/Documents/Data/Blast/')
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
	
	
	DicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2 = importData(path, specie1, specie2)
	#~ importData(path, specie1, specie2, idTypePerSpecie)
	
	DicoInfoBBH = {"Gene" : {}, "Transcript" : {}}
	DicoInfoBBH["Gene"] = GetDicoInfo(DicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBH, "Gene", specie1, specie2)
	DicoInfoBBH["Transcript"] = GetDicoInfo(DicoBBH, DicoHomology, DicoInfoG4Sp1, DicoInfoG4Sp2, DicoInfoBBH, "Transcript", specie1, specie2)
	
	WriteFile(DicoInfoBBH, path, "Gene", "Orthology")
	WriteFile(DicoInfoBBH, path, "Gene", "Homology")
	WriteFile(DicoInfoBBH, path, "Gene", "No_orthology")
	WriteFile(DicoInfoBBH, path, "Gene", "Paralogy")
	print "File for the gene level done"
	WriteFile(DicoInfoBBH, path, "Transcript", "No_orthology")
	WriteFile(DicoInfoBBH, path, "Transcript", "Orthology")
	WriteFile(DicoInfoBBH, path, "Transcript", "Homology")
	WriteFile(DicoInfoBBH, path, "Transcript", "Paralogy")
	print "File for the transcript level done"
#----------------------------------------------------------------------#	

main()
