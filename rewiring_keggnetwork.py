# -*- coding: utf-8 -*-
"""
Created on Thu May 14 10:28:45 2020
@author: Swapnil Keshari
Summary: This file rewires the network form the Dict 2 which is imported from dataprocessing_1 file
"""

import networkx as nx
import matplotlib.pyplot as plt
import networkConstructor as nc
from dataprocessing_1 import Dict2
G = nx.DiGraph()
# Removed Lines 150, 151, 152 from networkConstructor
#Creates a network
code='04060'
aliasDict, dict1, dict2={}, {}, {} # set up dicts for reading KEGG files
nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
coder=str('ko'+code) 
nc.uploadKEGGcodes([coder], G, dict2)
coder=str('hsa'+code)
nc.uploadKEGGcodes_hsa([coder], G,dict1, dict2)

#Check for common Genes in the graph and the data (CSV) genes
Gene = set(Dict2.keys())
Gene1 = set (G.nodes())
Gene2 = Gene1.intersection(Gene)

nx.write_gml(G, "04060_initial.gml")

removeNodeList=[]
#Creating a list of gene which has phosporylation position attached to it
for gene in Gene2:
    for seq in Dict2[gene].keys():
        if Dict2[gene][seq].keys() != None:
            removeNodeList.append(gene)
            break
#Rewiring program starts
for gene in removeNodeList:
    print gene
    #All the different sequence will be in parallel to each other
    for seq in Dict2[gene].keys():
        print seq
        monoList=[]#Used to store monophosphorylations
        diList=[]#Used to store diphosphorylations
        triList=[]#Used to store triphosphorylations
        removeDiList=[]#Used to store diphosphorylations node which will be removed from diList
        removeTriList=[]#Used to store monophosphorylations which will be removed from triList
        for phospho in Dict2[gene][seq].keys():
            if '1x' in phospho:
                monoList.append(str(gene)+"_"+str(seq)+"_"+str(phospho))
            elif '2x' in phospho:
                diList.append(str(gene)+"_"+str(seq)+"_"+str(phospho))            
            elif '3x' in phospho:
                triList.append(str(gene)+"_"+str(seq)+"_"+str(phospho))
        G.add_nodes_from(monoList)
        G.add_nodes_from(diList)
        G.add_nodes_from(triList)
        #Checking if monophosphorylations exist
        if monoList !=[]:
            #All mono are attached to the source nodes
            for mono in monoList:
                for start in G.predecessors(gene):
                    edge1=G.get_edge_data(start,gene)["signal"]
                    if edge1=="i":
                        G.add_edge(start,mono,signal="i")
                    else:
                        G.add_edge(start,mono,signal="a")
                new1=mono.split('_1x')
                temp=[]
                temp.append(new1[1])
                v1=set(temp)
                if diList!=[]:
                    isSubsetofDi=0
                    for di in diList:
                        new2=di.split('_2x')
                        v2=set(new2[1].split('-'))
                        if v1.issubset(v2):#Checks if the mono are present in di or not
                            G.add_edge(mono,di,signal="a")
                            isSubsetofDi=isSubsetofDi+1
                            if triList!=[]:
                                isSubsetofTri=0
                                for tri in triList:
                                    new3=tri.split('_3x')
                                    v3=set(new3[1].split('-'))
                                    if v2.issubset(v3):#Checks if the same di is present in the tri
                                        G.add_edge(di,tri,signal="a")
                                        isSubsetofTri=isSubsetofTri+1
                                        for finish in G.successors(gene):
                                            edge2=G.get_edge_data(gene,finish)["signal"]
                                            if edge2=="i":
                                                G.add_edge(tri,finish,signal="i")
                                            else:
                                                G.add_edge(tri,finish,signal="a")
                                        removeTriList.append(tri)
                                if isSubsetofTri==0:#If di is not the subset of tri we join it to finish
                                    for finish in G.successors(gene):
                                        edge2=G.get_edge_data(gene,finish)["signal"]
                                        if edge2=="i":
                                            G.add_edge(di,finish,signal="i")
                                        else:
                                            G.add_edge(di,finish,signal="a")
                                    removeDiList.append(di)#Add it to the remove node list so that is not added to source later
                            else:#If there is no tri list
                                for finish in G.successors(gene):
                                    edge2=G.get_edge_data(gene,finish)["signal"]
                                    if edge2=="i":
                                        G.add_edge(di,finish,signal="i")
                                    else:
                                        G.add_edge(di,finish,signal="a")
                                removeDiList.append(di)#Add it to the remove node list so that is not added to source later
                    if isSubsetofDi==0:#Even though diList is there but mono isn't the subset of it
                        if triList!=[]:#We check that is it present in tri list
                            isSubsetofTri=0
                            for tri in triList:
                                new3=tri.split('_3x')
                                v3=set(new3[1].split('-'))
                                if v1.issubset(v3):#If mono is the subset of tri
                                    G.add_edge(mono,tri,signal="a")
                                    isSubsetofTri=isSubsetofTri+1
                                    for finish in G.successors(gene):
                                        edge2=G.get_edge_data(gene,finish)["signal"]
                                        if edge2=="i":
                                            G.add_edge(tri,finish,signal="i")
                                        else:
                                            G.add_edge(tri,finish,signal="a")
                                    removeTriList.append(tri)#Add it to the remove node list so that is not added to source later
                            if isSubsetofTri==0:#Even if tri is present but mono is not the subset of tri
                                for finish in G.successors(gene):
                                    edge2=G.get_edge_data(gene,finish)["signal"]
                                    if edge2=="i":
                                        G.add_edge(mono,finish,signal="i")
                                    else:
                                        G.add_edge(mono,finish,signal="a")
                        else:#if it is not in the subset of di but tri doen't exist
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(mono,finish,signal="i")
                                else:
                                    G.add_edge(mono,finish,signal="a")                                    
                elif triList!=[]:#if there is no di list but tri list exist
                    isSubsetofTri=0
                    for tri in triList:
                        new3=tri.split('_3x')
                        v3=set(new3[1].split('-'))
                        if v1.issubset(v3):#If mono is in the subset of tri
                            G.add_edge(mono,tri,signal="a")
                            isSubsetofTri=isSubsetofTri+1
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(tri,finish,signal="i")
                                else:
                                    G.add_edge(tri,finish,signal="a")
                            removeTriList.append(tri)#Add it to the remove node list so that is not added to source later
                    if isSubsetofTri==0:#Even if tri exist but mono isn't in the subset of tri
                        for finish in G.successors(gene):
                            edge2=G.get_edge_data(gene,finish)["signal"]
                            if edge2=="i":
                                G.add_edge(mono,finish,signal="i")
                            else:
                                G.add_edge(mono,finish,signal="a")
                else:#Only mono list exist
                    for finish in G.successors(gene):
                        edge2=G.get_edge_data(gene,finish)["signal"]
                        if edge2=="i":
                            G.add_edge(mono,finish,signal="i")
                        else:
                            G.add_edge(mono,finish,signal="a")
            diList=[x for x in diList if x not in removeDiList] #Remove all the diList elements which have mono as source 
        #if no monophosphrylation then checking existance of diphosphorylation
        elif diList!=[]:#Elements which don't have mono as the source present
            for di in diList:
                for start in G.predecessors(gene):
                    edge1=G.get_edge_data(start,gene)["signal"]
                    if edge1=="i":
                        G.add_edge(start,di,signal="i")
                    else:
                        G.add_edge(start,di,signal="a")
                new2=di.split('_2x')
                v2=set(new2[1].split('-'))
                if triList!=[]:#If trilist exist
                    isSubsetofTri=0
                    for tri in triList:
                        new3=tri.split('_3x')
                        v3=set(new3[1].split('-'))
                        if v2.issubset(v3):#If di exist in tri
                            G.add_edge(di,tri,signal="a")
                            isSubsetofTri=isSubsetofTri+1
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(tri,finish,signal="i")
                                else:
                                    G.add_edge(tri,finish,signal="a")
                            removeTriList.append(tri)#Add it to the remove node list so that is not added to source later
                    if isSubsetofTri==0:#Even though tri exist but still it is not in any tri
                        for finish in G.successors(gene):
                            edge2=G.get_edge_data(gene,finish)["signal"]
                            if edge2=="i":
                                G.add_edge(mono,finish,signal="i")
                            else:
                                G.add_edge(mono,finish,signal="a")                                
                else:#If tri doesn't exist
                    for finish in G.successors(gene):
                        edge2=G.get_edge_data(gene,finish)["signal"]
                        if edge2=="i":
                            G.add_edge(di,finish,signal="i")
                        else:
                            G.add_edge(di,finish,signal="a")        
            triList=[x for x in triList if x not in removeTriList]   #Remove all the trList elements which have mono or di as source           
        #Else triphosphorylation shall exist since we have some type of phosphorylation attached
        elif triList!=[]:#Elements which donot have and di or mono as source
            for tri in triList:
                for start in G.predecessors(gene):
                    edge1=G.get_edge_data(start,gene)["signal"]
                    if edge1=="i":
                        G.add_edge(start,tri,signal="i")
                    else:
                        G.add_edge(start,tri,signal="a")
                for finish in G.successors(gene):
                    edge2=G.get_edge_data(gene,finish)["signal"]
                    if edge2=="i":
                        G.add_edge(tri,finish,signal="i")
                    else:
                        G.add_edge(tri,finish,signal="a")

G.remove_nodes_from(removeNodeList)                        
nx.write_gml(G, "04060_final.gml")
