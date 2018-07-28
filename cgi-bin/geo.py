#geo.py
#Geometry lookup table
import xml.etree.ElementTree as etree
import numpy as np

def ListGeometries():
	return list(ParseXML().keys())

def GetTab(Geometry):
	return ParseXML()[Geometry]

def Lookup(tab,SeekProp,Val,ResProp):
	return tab

def ParseXML():
	root = (etree.parse('GeometryData.xml')).getroot()
	tabs={}
	#the property names
	for t in range(len(root)):
		props = list(root[t][0].attrib.keys())
		nrows = len(root[t])
		tab = {}
		for p in range(len(props)):	
			tab[(props[p])] = np.zeros(nrows)
		for r in range(nrows):	
			for p in range(len(props)):
				tab[props[p]][r] = root[t][r].attrib[props[p]]
		tabs[root[t].attrib['Name']] = tab
	return tabs
