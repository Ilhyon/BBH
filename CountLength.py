#!/usr/bin/env python
# -*- coding: utf-8 -*-:v

import re

Dico = {"Gene" : {"Total" : 0}, "Transcript" : {"Total" : 0}}

with open("/home/local/USHERBROOKE/vana2406/Documents/Data/Human/All/HS_length_Genome_Transcriptome.txt") as f: # file opening
		content = f.read()
		lines = content.split('\n')
		for l in lines: #parcour de toute les lignes
			if l :
				words=l.split('\t')
				idGene = words[0]
				idTranscript = words[1]
				startGene = int(words[2])
				endGene = int(words[3])
				startTranscript = int(words[4])
				endTranscript = int(words[5])
				strand = int(words[6])
				if strand == 1 :
					lengthGene =  startGene - endGene
					lengthTranscript = startTranscript - endTranscript
					if idGene not in Dico["Gene"] :
						Dico["Gene"][idGene] = lengthGene
						Dico["Gene"]["Total"] += lengthGene
					if idTranscript not in Dico["Transcript"] :
						Dico["Transcript"][idGene] = lengthTranscript
						Dico["Transcript"]["Total"] += lengthTranscript
				else :
					lengthGene =  endGene - startGene
					lengthTranscript = endTranscript - startTranscript
					if idGene not in Dico["Gene"] :
						Dico["Gene"][idGene] = lengthGene
						Dico["Gene"]["Total"] += lengthGene
					if idTranscript not in Dico["Transcript"] :
						Dico["Transcript"][idGene] = lengthTranscript
						Dico["Transcript"]["Total"] += lengthTranscript

print Dico["Gene"]["Total"], Dico["Transcript"]["Total"]
