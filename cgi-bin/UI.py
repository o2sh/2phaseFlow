#!/usr/bin/env python3
import cgi
import cgitb
import os
import re
import Graph
import sat
import geo

def main():
	Param = {'BendDiameter': 0.100,'MassFlux': 300.0,'VaporQuality': 0.1,'Tsat': 13., 'Dia': 0.014,'QFlux':7500.0,'Fluid':'R134a','P_drop':10, 'Geometry':'UBend_Horizontal_flow'}	#defaults
	ParseQStr(Param) # you write the new parameters

	#if request type POST, read form and redirect
	if os.environ['REQUEST_METHOD']=='POST': 
		RedirectViaFormData(Param) # I don't Know, really
	else:
		WritePage(Param) #calculate everything launch getImage

#read the query string to numbers
def ParseQStr(Param):
	p = re.compile(r"([A-Za-z0-9]+)=([^;&]+)(?:[;&]|$)")
	l = p.findall(os.environ['QUERY_STRING']) # you get the new value
	for i in range(len(l)):
		TryValue(Param,l[i][0],l[i][1])

#attempt to convert a string to value matching the current type in Param
def TryValue(Param,key,value):
	try:
		if(type(Param[key])==float):
			Param[key] = float(value)
		elif(type(Param[key])==str):
			Param[key] = value
	except:
		pass

def RedirectViaFormData(Param):
	form = cgi.FieldStorage()		#data from the http POST
	for key in Param:
		TryValue(Param,key,form.getfirst(key, ""))
	#issue a html redirection encoding the new data into the URI
	print("Content-Type: text/html")		
	print("Refresh: 0; url=./UI.py?"+MakeQStr(Param))
	print()# end CGI header and send

#make a query string from numbers
def MakeQStr(Param):
	return "BendDiameter={BendDiameter:+.5e};MassFlux={MassFlux:+.5e};VaporQuality={VaporQuality:+.5e};Tsat={Tsat:+.5e};Dia={Dia:+.5e};QFlux={QFlux:+.5e};Fluid={Fluid}; P_drop={P_drop:+.5e};Geometry={Geometry}".format(**Param)

def WritePage(Param):
    Param['QStr']=MakeQStr(Param) #query string
    Param['ImgPath']='img_cache/'+Param['QStr']   #You create an .png corresponding to the pattern
    Param['ErrStr']=""
    Graph.GetImage(Param)# in this function you modify param['imgPath']
    Fluids = sat.ListFluids() #get the list of fluids
    Param['FluidOpts']=OptsList(Fluids,Param['Fluid']) #build the fluid list
    Geometries = geo.ListGeometries() #get the list of fluids
    Param['GeometryOpts']=OptsList(Geometries,Param['Geometry'])
    fin = open('ong.html');
    contents = fin.read();
    fin.close()
    print("Content-Type: text/html")
    print() #end CGI header
    print(contents.format(**Param)) #render html from template file

#make a list of <options> for an html <select>, including one selected
def OptsList(List,Select):
	Opts=''
	for i in range(len(List)):
		if List[i]==Select:
			Opts += """<option value="{0}" selected="selected">{0}</option>""".format(List[i])
		else:
			Opts += """<option value="{0}" >{0}</option>""".format(List[i])
	return Opts

#execute main when called as a script, but not as Python module
if __name__ == "__main__":
	cgitb.enable()				#dump errors to the browser
	main()
