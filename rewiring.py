# -*- coding: utf-8 -*-
"""
Created on Thu May 14 10:28:45 2020

@author: 98kes
"""

import networkx as nx
import matplotlib.pyplot as plt
from dataprocessing_1 import Dict2
G = nx.DiGraph()

Gene = Dict2.keys()
G.add_nodes_from(Gene)

for index in range(len(Gene)):
    if index != len(Gene)-1:
        G.add_edges_from([(Gene[index],Gene[index+1])],signal='i')
#    else:
#        G.add_edges_from([(Gene[index],Gene[0])],signal='i')

nx.write_gml(G, "test.gml") 
       
#G.add_edges_from([(1,3)],signal='i')
#G.add_edges_from([(2,3),(3,4)],signal='a')
#nx.write_gml(G, "test.gml")

#GeneDict={1:None,2:["3P","6P","2P","3P-6P","3P-6P-2P"],3:["4P","5P","3P-4P"],4:None}

removeNodeList=[]

for gene in Dict2.keys():
    for seq in Dict2[gene].keys():
        if Dict2[gene][seq].keys() != None:
            removeNodeList.append(gene)
            break

#for rm in removeNodeList:
#    for value in GeneDict[rm]:
#        if '-' not in value:
#            for start in G.predecessors(rm):
#                edge1=G.get_edge_data(start,rm)["signal"]
#                if edge1=="i":
#                    G.add_edge(start,value,signal="i")
#                else:
#                    G.add_edge(start,value,signal="a")
#            for finish in G.successors(rm):
#                edge2=G.get_edge_data(rm,finish)["signal"]
#                if edge2=="i":
#                    G.add_edge(value,finish,signal="i")
#                else:
#                    G.add_edge(value,finish,signal="a")




for gene in removeNodeList:
    print gene
    for seq in Dict2[gene].keys():
        print seq
        monoList=[]
        diList=[]
        triList=[]
        removeDiList=[]
        removeTriList=[]
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
        if monoList !=[]:
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
                    flag2=0
                    for di in diList:
                        new2=di.split('_2x')
                        v2=set(new2[1].split('-'))
                        if v1.issubset(v2):
                            G.add_edge(mono,di,signal="a")
                            flag2=flag2+1
                            if triList!=[]:
                                flag3=0
                                for tri in triList:
                                    new3=tri.split('_3x')
                                    v3=set(new3[1].split('-'))
                                    if v2.issubset(v3):
                                        G.add_edge(di,tri,signal="a")
                                        flag3=flag3+1
                                        for finish in G.successors(gene):
                                            edge2=G.get_edge_data(gene,finish)["signal"]
                                            if edge2=="i":
                                                G.add_edge(tri,finish,signal="i")
                                            else:
                                                G.add_edge(tri,finish,signal="a")
                                        removeTriList.append(tri)
                                if flag3==0:
                                    for finish in G.successors(gene):
                                        edge2=G.get_edge_data(gene,finish)["signal"]
                                        if edge2=="i":
                                            G.add_edge(di,finish,signal="i")
                                        else:
                                            G.add_edge(di,finish,signal="a")
                                    removeDiList.append(di)
                            else:
                                for finish in G.successors(gene):
                                    edge2=G.get_edge_data(gene,finish)["signal"]
                                    if edge2=="i":
                                        G.add_edge(di,finish,signal="i")
                                    else:
                                        G.add_edge(di,finish,signal="a")
                                removeDiList.append(di)
                    if flag2==0:
                        if triList!=[]:
                            flag3=0
                            for tri in triList:
                                new3=tri.split('_3x')
                                v3=set(new3[1].split('-'))
                                if v1.issubset(v3):
                                    G.add_edge(mono,tri,signal="a")
                                    flag3=flag3+1
                                    for finish in G.successors(gene):
                                        edge2=G.get_edge_data(gene,finish)["signal"]
                                        if edge2=="i":
                                            G.add_edge(tri,finish,signal="i")
                                        else:
                                            G.add_edge(tri,finish,signal="a")
                                    removeTriList.append(tri)
                            if flag3==0:
                                for finish in G.successors(gene):
                                    edge2=G.get_edge_data(gene,finish)["signal"]
                                    if edge2=="i":
                                        G.add_edge(mono,finish,signal="i")
                                    else:
                                        G.add_edge(mono,finish,signal="a")
                        else:
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(mono,finish,signal="i")
                                else:
                                    G.add_edge(mono,finish,signal="a")                                    
                elif triList!=[]:
                    flag3=0
                    for tri in triList:
                        new3=tri.split('_3x')
                        v3=set(new3[1].split('-'))
                        if v1.issubset(v3):
                            G.add_edge(mono,tri,signal="a")
                            flag3=flag3+1
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(tri,finish,signal="i")
                                else:
                                    G.add_edge(tri,finish,signal="a")
                            removeTriList.append(tri)
                    if flag3==0:
                        for finish in G.successors(gene):
                            edge2=G.get_edge_data(gene,finish)["signal"]
                            if edge2=="i":
                                G.add_edge(mono,finish,signal="i")
                            else:
                                G.add_edge(mono,finish,signal="a")
                else:
                    for finish in G.successors(gene):
                        edge2=G.get_edge_data(gene,finish)["signal"]
                        if edge2=="i":
                            G.add_edge(mono,finish,signal="i")
                        else:
                            G.add_edge(mono,finish,signal="a")
            diList=[x for x in diList if x not in removeDiList]  
        elif diList!=[]:
            for di in diList:
                for start in G.predecessors(gene):
                    edge1=G.get_edge_data(start,gene)["signal"]
                    if edge1=="i":
                        G.add_edge(start,di,signal="i")
                    else:
                        G.add_edge(start,di,signal="a")
                new2=di.split('_2x')
                v2=set(new2[1].split('-'))
                if triList!=[]:
                    flag3=0
                    for tri in triList:
                        new3=tri.split('_3x')
                        v3=set(new3[1].split('-'))
                        if v2.issubset(v3):
                            G.add_edge(di,tri,signal="a")
                            flag3=flag3+1
                            for finish in G.successors(gene):
                                edge2=G.get_edge_data(gene,finish)["signal"]
                                if edge2=="i":
                                    G.add_edge(tri,finish,signal="i")
                                else:
                                    G.add_edge(tri,finish,signal="a")
                            removeTriList.append(tri)
                    if flag3==0:
                        for finish in G.successors(gene):
                            edge2=G.get_edge_data(gene,finish)["signal"]
                            if edge2=="i":
                                G.add_edge(mono,finish,signal="i")
                            else:
                                G.add_edge(mono,finish,signal="a")                                
                else:
                    for finish in G.successors(gene):
                        edge2=G.get_edge_data(gene,finish)["signal"]
                        if edge2=="i":
                            G.add_edge(di,finish,signal="i")
                        else:
                            G.add_edge(di,finish,signal="a")        
            triList=[x for x in triList if x not in removeTriList]             
        elif triList!=[]:
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
nx.write_gml(G, "test1.gml")
#nx.draw(G,with_labels=True)
#plt.savefig("graph.png")
#plt.show()






#    if addNodeList != []:
#        for add in addNodeList:
#            for element in add.split('-'):
#                G.add_edge(element,add,signal="a")
#                for start in G.predecessors(rm):
#                    edge1=G.get_edge_data(start,rm)["signal"]
#                    if edge1=="i":
#                        G.add_edge(start,element,signal="i")
#                    else:
#                        G.add_edge(start,element,signal="a")
#            for finish in G.successors(rm):
#                edge2=G.get_edge_data(rm,finish)["signal"]
#                if edge2=="i":
#                    G.add_edge(add,finish,signal="i")
#                else:
#                    G.add_edge(add,finish,signal="a")
#    else:
#        for value in GeneDict[rm]:
#            if value not in addNodeList:
#                for start in G.predecessors(rm):
#                    edge1=G.get_edge_data(start,rm)["signal"]
#                    if edge1=="i":
#                        G.add_edge(start,value,signal="i")
#                    else:
#                        G.add_edge(start,value,signal="a")
#                for finish in G.successors(rm):
#                    edge2=G.get_edge_data(rm,finish)["signal"]
#                    if edge2=="i":
#                        G.add_edge(value,finish,signal="i")
#                    else:
#                        G.add_edge(value,finish,signal="a")
#            
#
#G.remove_nodes_from(removeNodeList)
#
#nx.write_gml(G, "test1.gml")







                        
#for rm in removeNodeList:
#	for start in G.predecessors(rm):
#		edge1=G.get_edge_data(start,rm)['signal']
#		if edge1=='i':
#			for element in rm.split('-'):
#				G.add_edge(start,element,signal='i')
#		else:
#			for element in rm.split('-'):
#				G.add_edge(start,element,signal='a')
#                
#	for finish in G.successors(rm):
#		edge2=G.get_edge_data(rm,finish)['signal']		
#		if edge2=='i':
#			for element in rm.split('-'):
#				G.add_edge(element,finish,signal='i')
#		else:
#			for element in rm.split('-'):
#				G.add_edge(element,finish,signal='a')
#	G.remove_node(rm)