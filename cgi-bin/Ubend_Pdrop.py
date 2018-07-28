
#### Code written by: Wilhelm Steffen Hench, Michael Gullberg

import numpy as np
import Flowpatternmap_Pdrop as flow
import cmath as m
import matplotlib.pyplot as plt

def fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry):
    ThetaDry=0
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2-0.01
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2
    return 0.215*((delta/(Dia-2*delta))**0.515)*(((DensL-DensV)*g*(Dia**2)*(1/STen))**(-0.298))*(ViscL/ViscV)**0.285*(DensL*(Ul**2)/(DensV*(Uv**2)))**-0.044
    


###################################################################
############### ANNULAR FLOW PRESSURE DROPS #######################    

#################
################# IF HORIZONTAL FLOW ANNULAR

def deltapUHANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,D,DeV,L,Al):     ##Friction coeff. 
    ThetaDry=0
    fUHANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.2*(Dia/D)**(1.313)*DeV**0.358)            ##Total pressure drop
    return 4*fUHANR*(L/Dia)*DensV*(Uv**2)*0.5
#################
################# IF VERTICAL DOWNFLOW ANNULAR
def deltapUVDANR(ViscL,ViscV,DensL,DensV,g,Dia,STen,Ul,Uv,G,x,D,L,Al):
    ThetaDry=0
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Refilm=((2*G*delta)/(ViscL))*((1-x)/(1-flow.epsi(x,DensV,DensL,STen,g,G)))
    Frfilm=(Ul**2)/(g*2*delta)
    fUVDANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.5*(Dia/D)**(1.339)*(Frfilm**(-0.193))*(Refilm**(0.555)))  ## Friction Coeff
    return 4*fUVDANR*(L/Dia)*DensV*(Uv**2)*0.5                               ## Total pressure drop

################## IF VERTICAL UPFLOW ANNULAR
def deltapUVUANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,G,x,D,L,Al):
    ThetaDry=0
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Frfilm=(Ul**2)/(g*2*delta)                       ##Film Froude
    fUVUANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+8*(Dia/D)**(1.178)*(Frfilm**0.196)*((Ul**2*DensV)/(Uv**2*DensL))**0.106) ##Friction Coeff
    return (4*fUVUANR*(L/Dia)*DensV*(Uv**2)*0.5)    ## Total pressure drop
    
###################################################################
########## SLUG (SLG) AND INTERMITTENT (IMT) FLOWS ################



#################
################# IF HORIZONTAL FLOW SLG AND IMT
                                          ## all liquid single phase Dean number
def deltapUHSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al):
    ThetaDry=0
    fUHL0=fSTL0*(1+103.19*(10**3)*((Dia/D)**2.405)*(DeL0**-0.653)) ## Final Modified Friction coeff SLG/IMT Horizontal
    deltapUHL0=4*fUHL0*(L/Dia)*0.5*DensL*(UL0**2)                  ##SLG/IMT horizontal pressure drop factor

    fUHANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.2*(Dia/D)**(1.313)*DeV**0.358)     ##Friction coeff Annular horizontal.
    deltapUHANR=4*fUHANR*(L/Dia)*DensV*(Uv**2)*0.5            ##Annular pressure drop factor

    return (deltapUHL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25)+deltapUHANR*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25))

#################
################# IF VERTICAL DOWNFLOW SLG AND IMT

def deltapUVDSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al):
    ThetaDry=0
    fUVDL0=fSTL0*(1+717.7*(10**5)*((Dia/D)**1.244)*(DeL0**(-1.461))) ##Final Modified Friction coeff SLG/IMT Vertical downflow
    deltapUVDL0=4*fUVDL0*(L/Dia)*0.5*DensL*(UL0**2)                ##SLG/IMT Vertival downflow pressure drop factor

    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Refilm=((2*G*delta)/(ViscL))*((1-x)/(1-flow.epsi(x,DensV,DensL,STen,g,G)))        ##Film reynolds
    Frfilm=(Ul**2)/(g*2*delta)                           ##Film Froude
    fUVDANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.5*(Dia/D)**(1.339)*(Frfilm**(-0.193))*(Refilm**(0.555)))  ## Friction Coeff annular Vertical downflow
    deltapUVDANR=4*fUVDANR*(L/Dia)*DensV*(Uv**2)*0.5                               ## Total pressure drop annular vertical downflow
    return deltapUVDL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25)+deltapUVDANR*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25) ## Total pressure drop

#################
################# IF VERTICAL UPFLOW SLG AND IMT
def deltapUVUSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al):
    ThetaDry=0
    fUVUL0=fSTL0*(1+4.09*(10**3)*((Dia/D)**1.002)*(DeL0**-0.381)) ##Final Modified Friction coeff SLG/IMT Vertical upflow
    deltapUVUL0=4*fUVUL0*(L/Dia)*0.5*DensL*(UL0**2)                ##SLG/IMT vertical upflow pressure drop factor

    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Frfilm=(Ul**2)/(g*2*delta)                       ##Film Froude
    fUVUANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+8*(Dia/D)**(1.178)*(Frfilm**0.196)*((Ul**2*DensV)/(Uv**2*DensL))**0.106) ##Friction Coeff annular vertical upflow
    deltapUVUANR=4*fUVUANR*(L/Dia)*DensV*(Uv**2)*0.5       ## Total pressure drop Annular Vertical upflow
    return deltapUVUL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25)+deltapUVUANR*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**(0.25) ##Total pressure drop

###################################################################
##########   STRATIFIED WAVY (STW) FLOW PATTERN    ################

#################
################# IF HORIZONTAL FLOW STW

def deltapUHSTW(A,WeFrl,G, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al):
    ThetaDry=((flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)*ThetaStrat 
    thetadry1=ThetaDry/(2*np.pi)

    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUHV=fSTV*(1+5.102*(10**6)*((Dia/D)**1.109)*(DeV**(-1.222)))
    fUHSTW=thetadry1*fUHV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.2*(Dia/D)**(1.313)*DeV**0.358)

    return 4*fUHSTW*(L/Dia)*DensV*(Uv**(2))/2

#################
################# IF VERTICAL DOWNFLOW STW

def deltapUVDSTW(A,WeFrl,G, Dia, x, ViscV, D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al):
    ThetaDry=((flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)*ThetaStrat
    thetadry1=ThetaDry/(2*np.pi)
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUVDV=fSTV*(1+8.39*((Dia/D)**(1.278))*(DeV**(0.057)))

    
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Refilm=((2*G*delta)/(ViscL))*((1-x)/(1-flow.epsi(x,DensV,DensL,STen,g,G)))
    Frfilm=(Ul**2)/(g*2*delta)
    fUVDSTW=thetadry1*fUVDV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al, ThetaDry)*(1+0.5*(Dia/D)**(1.339)*(Frfilm**(-0.193))*(Refilm**(0.555)))  ## Friction Coeff

    return 4*fUVDSTW*(L/Dia)*DensV*(Uv**(2))/2

#################
################# IF VERTICAL UPFLOW STW

def deltapUVUSTW(A,WeFrl,G, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al):
    
    ThetaDry=((flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)*ThetaStrat
    thetadry1=ThetaDry/(2*np.pi)
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUVUV=fSTV*(1+47.1*10**3*(Dia/D)**2.707*DeV**(-0.507))
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Frfilm=(Ul**2)/(g*2*delta)
    fUVUSTW=thetadry1*fUVUV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+8*(Dia/D)**(1.178)*(Frfilm**0.196)*((Ul**2*DensV)/(Uv**2*DensL))**0.106) ##Friction Coeff
    return 4*fUVUSTW*(L/Dia)*DensV*(Uv**(2))/2

###################################################################
######### SLUG STRATIFIED WAVY (SSW) FLOW PATTERN #################

#################
################# IF HORIZONTAL FLOW SSW

def deltapUHSSW(A,epsilonIA, fSTL0, DensL, UL0,WeFrl, Dia, x, ViscV, D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G):
    fUHL0=fSTL0*(1+103.19*(10**3)*((Dia/D)**2.405)*(DeL0**-0.653)) ## Final Modified Friction coeff SLG/IMT Horizontal
    deltapUHL0=4*fUHL0*(L/Dia)*0.5*DensL*(UL0**2)                  ##SLG/IMT horizontal pressure drop factor
    
    ThetaDry=(ThetaStrat)*(x/flow.xia(DensV,DensL,ViscL,ViscV))*((flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)  ## This is the change##It is important to note that in this case, flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen) =flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)(flow.xia(DensV,DensL,ViscL,ViscV)) and flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G) = flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)(flow.xia(DensV,DensL,ViscL,ViscV))
    thetadry1=ThetaDry/(2*np.pi)

    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUHV=fSTV*(1+5.102*(10**6)*((Dia/D)**1.109)*(DeV**(-1.222)))
    fUHSTW=thetadry1*fUHV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al, ThetaDry)*(1+0.2*(Dia/D)**(1.313)*DeV**0.358)
    deltapUHSTW1=4*fUHSTW*(L/Dia)*DensV*(Uv**(2))/2
    return deltapUHL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25+deltapUHSTW1*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25

################# IF VERTICAL DOWNFLOW SSW
def deltapUVDSSW(A,epsilonIA, fSTL0, DensL, UL0,WeFrl,  Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G):
    fUVDL0=fSTL0*(1+717.7*(10**5)*((Dia/D)**1.244)*(DeL0**(-1.461)))
    deltapUVDL0=4*fUVDL0*(L/Dia)*0.5*DensL*(UL0**2) 
	
    
    ThetaDry=(x/flow.xia(DensV,DensL,ViscL,ViscV))*((flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)*(ThetaStrat)  ## This is the change##It is important to note that in this case, flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen) =flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)(flow.xia(DensV,DensL,ViscL,ViscV)) and flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G) = flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)(flow.xia(DensV,DensL,ViscL,ViscV))
    thetadry1=ThetaDry/(2*np.pi)
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Refilm=((2*G*delta)/(ViscL))*((1-x)/(1-flow.epsi(x,DensV,DensL,STen,g,G)))
    Frfilm=(Ul**2)/(g*2*delta)

    fUVDV=fSTV*(1+8.39*(Dia/D)**1.278*DeV**0.057)
    fUVDSTW=thetadry1*fUVDV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al, ThetaDry)*(1+0.5*(Dia/D)**(1.339)*(Frfilm**(-0.193))*(Refilm**(0.555)))  ## Friction Coeff
    deltapUVDSTW1=4*fUVDSTW*(L/Dia)*DensV*(Uv**(2))/2
    return deltapUVDL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25+deltapUVDSTW1*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25


################# IF VERTICAL UPFLOW SSW

def deltapUVUSSW(A,epsilonIA, fSTL0, DensL, UL0, WeFrl, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G):
    fUVUL0=fSTL0*(1+4.09*(10**3)*((Dia/D)**1.002)*(DeL0**-0.381)) ##Final Modified Friction coeff SLG/IMT Vertical upflow
    deltapUVUL0=4*fUVUL0*(L/Dia)*0.5*DensL*(UL0**2) 

    
    ThetaDry=(x/flow.xia(DensV,DensL,ViscL,ViscV))*((flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-G)/(flow.wavy(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensL,DensV,g,WeFrl,G,STen)-flow.strat(flow.xia(DensV,DensL,ViscL,ViscV),A,Dia,DensV,DensL,ViscL,g,STen,G)))**(0.61)*(ThetaStrat)  ## This is the change##It is important to note that in this case, flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen) =flow.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)(flow.xia(DensV,DensL,ViscL,ViscV)) and flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G) = flow.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)(flow.xia(DensV,DensL,ViscL,ViscV))
    thetadry1=ThetaDry/(2*np.pi)
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Frfilm=(Ul**2)/(g*2*delta)

    fUVUV=fSTV*(1+47.1*10**3*(Dia/D)**2.707*DeV**(-0.507))
    fUVUSTW=thetadry1*fUVUV+(1-thetadry1)*fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+8*(Dia/D)**(1.178)*(Frfilm**0.196)*((Ul**2*DensV)/(Uv**2*DensL))**0.106) ##Friction Coeff
    deltapUVUSTW1= 4*fUVUSTW*(L/Dia)*DensV*(Uv**(2))/2
	
    return deltapUVUL0*(1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25+deltapUVUSTW1*(flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25


###################################################################
######### Mist (MST) Flow Pattern #################


################# IF HORIZONTAL FLOW MIST
def deltapUHMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L):
    epsilonH=(1+((1-x)/x)*(DensV/DensL))**(-1)
    DensMist=DensL*(1-epsilonH)+DensV*epsilonH
    ViscMist=ViscL*(1-x)+ViscV*x
    ReMist=(G*Dia)/ViscMist
    DeMST=(G*Dia/ViscMist)*(Dia/D)**(0.5)
    fSTMST=0.079/(ReMist**0.25)
    fUHMST=fSTMST*(1+300.35*(Dia/D)**1.623*DeMST**(-0.189))
    return 4*fUHMST*(L/Dia)*((G**2)/(2*DensMist))

################# IF VERTICAL DOWNFLOW

def deltapUVDMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L):
    epsilonH=(1+((1-x)/x)*(DensV/DensL))**(-1)
    DensMist=DensL*(1-epsilonH)+DensV*epsilonH
    ViscMist=ViscL*(1-x)+ViscV*x
    ReMist=(G*Dia)/ViscMist
    fSTMST=0.079/(ReMist**0.25)
    return 4*fSTMST*(L/Dia)*((G**2)/(2*DensMist))


################# IF Vertical Upflow MIST 
   
def deltapUVUMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L):

    epsilonH=(1+((1-x)/x)*(DensV/DensL))**(-1)
    DensMist=DensL*(1-epsilonH)+DensV*epsilonH
    ViscMist=ViscL*(1-x)+ViscV*x
    ReMist=(G*Dia)/ViscMist
    DeMST=((G*Dia)/ViscMist)*(Dia/D)**(0.5)
    fSTL0=0.079/(ReMist**0.25)

    fUVUMST=fSTL0*(1+24.491*(Dia/D)**1.090*DeMST**0.011)
    return 4*fUVUMST*(L/Dia)*G**2/(2*DensMist)


###################################################################
######### Dryout Flow #################


####Horizontal Flow Dryout####

def deltapUHDYT(A,ThetaStratdi,ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,qcrit,q):

    if patternbef==1:
	
        deltapUHDYT1=deltapUHANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi)-(x-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))/(flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))*(deltapUHANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi)-deltapUHMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q),D, L)) #deltapUHANR has to be calculated at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit) and deltapUHMST at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
		
    elif patternbef==2:
        
        deltapUHDYT1=deltapUHSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV,  D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-(x-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))/(flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))*(deltapUVUSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV,  D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-deltapUHMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q),D, L)) #deltapUHANR has to be calculated at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit) and deltapUHMST at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
    
    
    return deltapUHDYT1

### vertical downflow dryout####

def deltapUVDDYT(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,ThetaStratdi,q,qcrit,A):
	
    if patternbef==1:
        
        deltapUVDDYT1=deltapUVDANR(ViscL,ViscV,DensL,DensV,g,Dia,STen,Uldi,Uvdi,G,flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit),D,L,Aldi)-\
            ((x-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))/(flow.xde(G,Dia,DensV,DensL,STen,g,q,qcrit)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)))*\
            (deltapUVDANR(ViscL,ViscV,DensL,DensV,g,Dia,STen,Uldi,Uvdi,G,flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit),D,L,Aldi)-\
             deltapUVDMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,DensV,DensL,STen,g,q,qcrit),D, L))
             
            

    elif patternbef==2:
        deltapUVDDYT1=deltapUVDSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV, D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-\
            ((x- flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))/(flow.xde(G,Dia,DensV,DensL,STen,g,q,qcrit)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))) *(deltapUVDSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV, D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-\
            deltapUVDMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q),D, L))                                                                      
    return deltapUVDDYT1

### vertical upflow dryout

def deltapUVUDYT(A,ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,ThetaStratdi,q,qcrit):
    if patternbef==1:

        deltapUVUDYT1=deltapUVUANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,G,flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit),D,L,Aldi)-((x-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)))/(flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))*(deltapUVUANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,G,flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit),D,L,Aldi)-(deltapUVUMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q),D, L))) #deltapUHANR has to be calculated at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit) and deltapUHMST at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
        
    elif patternbef==2:
		
        deltapUVUDYT1=deltapUVUSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV,  D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-((x-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)))/(flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)-flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit))*(deltapUVUSTW(A,WeFrl,G, Dia, flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit), ViscV,  D, DeVdi, L, DensV, Uvdi, ThetaStratdi, ViscL, DensL, g, STen, Uldi, Aldi)-deltapUVUMST(DensV,DensL,G,Dia,ViscL,ViscV,flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q),D, L)) #deltapUHANR has to be calculated at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit) and deltapUHMST at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
    
    return deltapUVUDYT1
		
###################################################################
######### Stratified Flow Pattern #################	

####horizontal Flow Stratified#####

##if x>flow.xia(DensV,DensL,ViscL,ViscV)
def deltapUHSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L):
    fSTV=0.079/(Rev**0.25)
    fUHV=fSTV*(1+5.102*(10**6)*((Dia/D)**1.109)*(DeV**(-1.222)))
    ThetaDry=ThetaStrat
	
    fUHANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.2*(Dia/D)**(1.313)*DeV**0.358)  
	
    fUHSTD=ThetaStrat*fUHV+(1-ThetaStrat)*fUHANR
    return 4*fUHSTD*L/Dia*DensV*(Uv**2)/2

##if x<flow.xia(DensV,DensL,ViscL,ViscV)
def deltaUHSTDsmaller(x,G,fSTL0,DeL0,UL0,epsilonIA,ViscV,ViscL,DensL,DensV,g,STen,Ul,Uv,Al,Dia,D,DeV,ThetaStrat,Rev,L):

    fUHL0=fSTL0*(1+103.19*(10**3)*((Dia/D)**2.405)*(DeL0**-0.653)) ## Final Modified Friction coeff SLG/IMT Horizontal
    deltapUHL0=4*fUHL0*(L/Dia)*0.5*DensL*(UL0**2)                  ##SLG/IMT horizontal pressure drop factor

    return deltapUHL0*((1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)+deltapUHSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L)*((flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)


####Vertical Downflow Stratified#####	
### if x>flow.xia(DensV,DensL,ViscL,ViscV)
def deltapUVDSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x):
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUVDV=fSTV*(1+8.39*((Dia/D)**(1.278))*(DeV**(0.057)))
	
    ThetaDry=ThetaStrat
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Refilm=((2*G*delta)/(ViscL))*((1-x)/(1-flow.epsi(x,DensV,DensL,STen,g,G)))
    Frfilm=(Ul**2)/(g*2*delta)
    fUVDANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+0.5*(Dia/D)**(1.339)*(Frfilm**(-0.193))*(Refilm**(0.555)))  ## Friction Coeff
	
    fUVDSTD=ThetaStrat*fUVDV+(1-ThetaStrat)*fUVDANR
    return 4*fUVDSTD*L/Dia*DensV*(Uv**2)/2

###if x<flow.xia(DensV,DensL,ViscL,ViscV)
def deltapUVDSTDsmaller(ThetaStrat,ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,G,x,Rev,D,DeV,epsilonIA,DeL0,UL0,fSTL0,L):
	
    fUVDL0=fSTL0*(1+717.7*(10**5)*((Dia/D)**1.244)*(DeL0**(-1.461))) ##Final Modified Friction coeff SLG/IMT Vertical downflow
    deltapUVDL0=4*fUVDL0*(L/Dia)*0.5*DensL*(UL0**2) 
	
    return deltapUVDL0*((1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)+deltapUVDSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x)*((flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)
	
	
####Vertical Upflow Stratified####

##if x>flow.xia(DensV,DensL,ViscL,ViscV)
def deltapUVUSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x):
    ThetaDry=ThetaStrat
	
    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
    fSTV=0.079/(Rev**0.25)
    fUVUV=fSTV*(1+47.1*10**3*(Dia/D)**2.707*DeV**(-0.507))
	
    delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-ThetaDry))
    if(delta.imag != 0):
        delta=Dia/2
    delta=float(delta.real)
    if(delta > (Dia/2)):
            delta=Dia/2  
    Frfilm=(Ul**2)/(g*2*delta)                       ##Film Froude
    fUVUANR=fSTANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,ThetaDry)*(1+8*(Dia/D)**(1.178)*(Frfilm**0.196)*((Ul**2*DensV)/(Uv**2*DensL))**0.106) ##Friction Coeff

    fUVUSTD=ThetaStrat*fUVUV+(1-ThetaStrat)*fUVUANR
    return 4*fUVUSTD*L/Dia*DensV*(Uv**2)/2

##  if x<flow.xia(DensV,DensL,ViscL,ViscV)

def deltapUVUSTDsmaller(ThetaStrat,ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,epsilonIA,Rev,D,DeV,fSTL0,DeL0,L,UL0,G,x):

    fUVUL0=fSTL0*(1+4.09*(10**3)*((Dia/D)**1.002)*(DeL0**-0.381)) ##Final Modified Friction coeff SLG/IMT Vertical upflow
    deltapUVUL0=4*fUVUL0*(L/Dia)*0.5*DensL*(UL0**2)                ##SLG/IMT vertical upflow pressure drop factor
	
    return deltapUVUL0*((1-flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)+deltapUVUSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x)*((flow.epsi(x,DensV,DensL,STen,g,G)/epsilonIA)**0.25)

def computePgrad(G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L,D,orientation):

    #orientation 2=horizontal; 3=Vertical Downflow; 4=Vertical Upflow
    #D is bend diameter
    g=9.81
    qcrit=0.131*DensV**0.5*(HV-HL)*(g*(DensL-DensV)*STen)**0.25


    vaporQualities=list(map(lambda x: x/500.0, range(1, 500)))
    Pgradient=[]
    patternloopbef=0
    patternbef=0
    for x in vaporQualities:
        pattern = flow.computePattern(x,G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q)
        (Pdrop,patternloopbef,patternbef) = computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,D,g,G,qcrit,q,HV,HL,orientation,patternloopbef,patternbef)
        Pgrad = Pdrop/L
        Pgradient.append(Pgrad)
    return (Pgradient)

def computePgradOne(G,x,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L,D,orientation):
    g=9.81
    hlv=HV-HL
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25
    patternloopbef=0
    patternbef=0
    pattern = flow.computePattern(x,G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q)
    (Pdrop,patternloopbef,patternbef) = computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,D,g,G,qcrit,q,HV,HL,orientation,patternloopbef,patternbef)
    Pgrad = Pdrop/L                                                                                                                                                #PGRADIENT OR PDROP ???
    return int(Pgrad)

def computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,D,g,G,qcrit,q,HV,HL,orientation,patternloopbef,patternbef):

    epsilonIA=flow.epsi(flow.xia(DensV,DensL,ViscL,ViscV),DensV,DensL,STen,g,G)
    Ul=(G*(1-x))/(DensL*(1-flow.epsi(x,DensV,DensL,STen,g,G)))                        ##Liquid speed
    Uv=(G*x)/(DensV*flow.epsi(x,DensV,DensL,STen,g,G))                                ##Vapor speed
    WeL=(1/STen)*Dia*DensL*(Ul**2)                            ##Liquid Weber
    Al=((np.pi*Dia**2)/4*(1-flow.epsi(x,DensV,DensL,STen,g,G)))                             ##Liquid Cross-sectional area
    ThetaStrat = (2*np.pi)-2*(np.pi*(1-flow.epsi(x,DensV,DensL,STen,g,G))+((np.pi*3/2)**(1/3))*(1-2*(1-flow.epsi(x,DensV,DensL,STen,g,G))+(1-flow.epsi(x,DensV,DensL,STen,g,G))**(1/3)-flow.epsi(x,DensV,DensL,STen,g,G)**(1/3))-(1/200)*(1-flow.epsi(x,DensV,DensL,STen,g,G))*flow.epsi(x,DensV,DensL,STen,g,G)*(1-2*(1-flow.epsi(x,DensV,DensL,STen,g,G)))*(1+4*((1-flow.epsi(x,DensV,DensL,STen,g,G))**2+flow.epsi(x,DensV,DensL,STen,g,G)**2)))#Biberg (1999)
    DeV=((G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G)))*(Dia/D)**0.5                 ##Dean number, explains inertial to viscous forces
    UL0=(G/DensL)           ##All-liquid single-phase velocity
    ReL0=(G*Dia)/(ViscL)      ##Liquid single phase reynolds number
    fSTL0=0.079/(ReL0**0.25)##Friction coeff not modified
    DeL0=ReL0*(Dia/D)**0.5     ##Dean number all liquid single phase


	###### values used for dryout calculations#####


    Rev=(G*Dia*x)/(ViscV*flow.epsi(x,DensV,DensL,STen,g,G))
                         #Gravitational acceleration [m/s2]
    A=np.pi*Dia**2/4            #Cross-sectional area [m2]
    hlv=HV-HL                   #Latent heat of vaporisation [J/kg]
    Wel=G**2*Dia/(DensL*STen)   #Liquid Weber number [...]
    Frl=G**2/(DensL**2*g*Dia)   #liquid Froude number [...]
    WeFrl=(g*Dia**2*DensL)/STen #(Wel/Frl)
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25

    xdi=flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)
    xde=flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
    
    Uldi=(G*(1-xdi))/(DensL*(1-flow.epsi(xdi,DensV,DensL,STen,g,G)))                        ##Liquid speed at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)
    Uvdi=(G*xdi)/(DensV*flow.epsi(xdi,DensV,DensL,STen,g,G))                                ##Vapor speed at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)
    #Ulde=(G*(1-xde))/(DensL*(1-flow.epsi(xde,DensV,DensL,STen,g,G)))                        ##Liquid speed at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
    Uvde=(G*xde)/(DensV*flow.epsi(xde,DensV,DensL,STen,g,G))                                ##Vapor speed at flow.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
    DeVdi=((G*Dia*xdi)/(ViscV*flow.epsi(xdi,DensV,DensL,STen,g,G)))*(Dia/D)**0.5                 ##Dean number at flow.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)
    DeVde=((G*Dia*xde)/(ViscV*flow.epsi(xdi,DensV,DensL,STen,g,G)))*(Dia/D)**0.5
    Aldi=((np.pi*Dia**2)/4*(1-flow.epsi(xdi,DensV,DensL,STen,g,G)))

    ThetaStratdi = (2*np.pi)-2*(np.pi*(1-flow.epsi(xdi,DensV,DensL,STen,g,G))+((np.pi*3/2)**(1/3))*(1-2*(1-flow.epsi(xdi,DensV,DensL,STen,g,G))+(1-flow.epsi(xdi,DensV,DensL,STen,g,G))**(1/3)-flow.epsi(xdi,DensV,DensL,STen,g,G)**(1/3))-(1/200)*(1-flow.epsi(xdi,DensV,DensL,STen,g,G))*flow.epsi(xdi,DensV,DensL,STen,g,G)*(1-2*(1-flow.epsi(xdi,DensV,DensL,STen,g,G)))*(1+4*((1-flow.epsi(xdi,DensV,DensL,STen,g,G))**2+flow.epsi(xdi,DensV,DensL,STen,g,G)**2)))#Biberg (1999)

    ## orientation decides if straight tube(1) or horizontal flow(2), vertical downflow(3) or vertical upflow(4)

    if(pattern!=patternloopbef and pattern==5):
        patternbef=patternloopbef
    
    if(pattern==1):
        if(orientation==2):
            DeltaP=deltapUHANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,D,DeV,L,Al)
        elif(orientation==3):
            DeltaP=deltapUVDANR(ViscL,ViscV,DensL,DensV,g,Dia,STen,Ul,Uv,G,x,D,L,Al)
        elif(orientation==4):
            DeltaP=deltapUVUANR(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,G,x,D,L,Al)
    elif(pattern==2):      
        if(orientation==2):
            DeltaP=deltapUHSTW(A,WeFrl,G, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al)
        elif(orientation==3):
            DeltaP=deltapUVDSTW(A,WeFrl,G, Dia, x, ViscV, D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al)
        elif(orientation==4):
            DeltaP=deltapUVUSTW(A,WeFrl,G, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, DensL, g, STen, Ul, Al)
    elif(pattern==8):
        
        if(orientation==2):
            DeltaP=deltapUHSSW(A,epsilonIA, fSTL0, DensL, UL0,WeFrl, Dia, x, ViscV, D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G)
        elif(orientation==3):
            DeltaP=deltapUVDSSW(A,epsilonIA, fSTL0, DensL, UL0,WeFrl,  Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G)
        elif(orientation==4):
            DeltaP=deltapUVUSSW(A,epsilonIA, fSTL0, DensL, UL0, WeFrl, Dia, x, ViscV,  D, DeV, L, DensV, Uv, ThetaStrat, ViscL, g, STen, Ul, Al, DeL0, G)
        
    elif(pattern==6 or pattern==7):
        if(orientation==2):
            DeltaP=deltapUHSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al)
        elif(orientation==3):
            DeltaP=deltapUVDSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al)
        elif(orientation==4):
            DeltaP=deltapUVUSLGIMT(UL0,ReL0,DeL0,fSTL0,ViscV,ViscL,DensL,DensV,DeV,g,Dia,STen,Ul,Uv,epsilonIA,G,x,D,L,Al)
    elif(pattern==4):
        if(orientation==2):
            DeltaP=deltapUHMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L)
        elif(orientation==3):
            DeltaP=deltapUVDMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L)
        elif(orientation==4):
            DeltaP=deltapUVUMST(DensV,DensL,G,Dia,ViscL,ViscV,x,D, L)
    elif(pattern==5):
        
        if(orientation==2):
            DeltaP=deltapUHDYT(A,ThetaStratdi,ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,qcrit,q)
        elif(orientation==3):
            DeltaP=deltapUVDDYT(ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,ThetaStratdi,q,qcrit,A)
        elif(orientation==4):
            DeltaP=deltapUVUDYT(A,ViscV,ViscL,DensL,DensV,g,Dia,STen,Uldi,Uvdi,D,DeVdi,L,Aldi,x,G,patternbef,WeFrl,ThetaStratdi,q,qcrit)
    elif(pattern==3):        
        if(orientation==2):
            if(x>flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltapUHSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L)
            elif(x<flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltaUHSTDsmaller(x,G,fSTL0,DeL0,UL0,epsilonIA,ViscV,ViscL,DensL,DensV,g,STen,Ul,Uv,Al,Dia,D,DeV,ThetaStrat,Rev,L)
        elif(orientation==3):
            if(x>flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltapUVDSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x)
            elif(x<flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltapUVDSTDsmaller(ThetaStrat,ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,G,x,Rev,D,DeV,epsilonIA,DeL0,UL0,fSTL0,L)
        elif(orientation==4):
            if(x>flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltapUVUSTDbigger(ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,D,DeV,ThetaStrat,Rev,L,G,x)
            elif(x<flow.xia(DensV,DensL,ViscL,ViscV)):
                DeltaP=deltapUVUSTDsmaller(ThetaStrat,ViscV,ViscL,DensL,DensV,g,Dia,STen,Ul,Uv,Al,epsilonIA,Rev,D,DeV,fSTL0,DeL0,L,UL0,G,x)

    patternloopbef=pattern
    return (DeltaP,patternloopbef,patternbef)   
