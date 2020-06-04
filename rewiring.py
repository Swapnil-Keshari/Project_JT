# -*- coding: utf-8 -*-
"""
Created on Thu May 14 10:28:45 2020

@author: 98kes
"""

import networkx as nx
G = nx.DiGraph()

Gene = [1,2,3,4]
G.add_nodes_from(Gene)
G.add_edges_from([(1,3)],signal='i')
G.add_edges_from([(2,3),(3,4)],signal='a')
nx.write_gml(G, "test.gml")

GeneDict={1:None,2:["3P","6P","2P","3P-6P-2P"],3:["3P","5P","3P-4P"],4:None}

removeNodeList=[]

for key in GeneDict.keys():
    if GeneDict[key] != None:
        removeNodeList.append(key)

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




for rm in removeNodeList:
    addNodeList=[]
    for value in GeneDict[rm]:
        if '-' in value:
            flag = True
            for element in value.split('-'):
                if element not in GeneDict[rm]:
                    flag= False
                    break
            if flag:
                addNodeList.append(value)
    if addNodeList != []:
        for add in addNodeList:
            for element in add.split('-'):
                G.add_edge(element,add,signal="a")
                for start in G.predecessors(rm):
                    edge1=G.get_edge_data(start,rm)["signal"]
                    if edge1=="i":
                        G.add_edge(start,element,signal="i")
                    else:
                        G.add_edge(start,element,signal="a")
            for finish in G.successors(rm):
                edge2=G.get_edge_data(rm,finish)["signal"]
                if edge2=="i":
                    G.add_edge(add,finish,signal="i")
                else:
                    G.add_edge(add,finish,signal="a")
    else:
        for value in GeneDict[rm]:
            if value not in addNodeList:
                for start in G.predecessors(rm):
                    edge1=G.get_edge_data(start,rm)["signal"]
                    if edge1=="i":
                        G.add_edge(start,value,signal="i")
                    else:
                        G.add_edge(start,value,signal="a")
                for finish in G.successors(rm):
                    edge2=G.get_edge_data(rm,finish)["signal"]
                    if edge2=="i":
                        G.add_edge(value,finish,signal="i")
                    else:
                        G.add_edge(value,finish,signal="a")
            

G.remove_nodes_from(removeNodeList)

nx.write_gml(G, "test1.gml")
                        
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