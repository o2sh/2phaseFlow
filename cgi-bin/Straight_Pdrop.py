import numpy as np
import matplotlib.pyplot as plt
import cmath as m
import Flowpatternmap_Pdrop as fpm
###### Straight tubes - Jesus Moreno ######


""" Check EpsilonIA computation, for x<xIA, stratified wavy gives complex pressure drop due to epsion bigger than epsilonIA
    Check the HAJAL delta that can become negatif
    Should use the formula for laminar reynolds too ? """

def computeDelta(DensL,DensV,ViscL,ViscV,STen,Dia,thetaDry,epsilon,Hajal=False):
    if(Hajal):
        #Computation of the film thickness using et Hajal et al. (2003b)
        A=np.pi*Dia**2/4
        Al=A*(1-epsilon)#Surface of the cross section occupied by liquid
        xIA=fpm.xia(DensV,DensL,ViscL,ViscV)        
        delta=Dia/2-m.sqrt(((Dia/2)**2)-(2*Al)/(2*np.pi-thetaDry))
        if(delta.imag != 0):
            delta=Dia/2
        delta=float(delta.real)
        if(delta > (Dia/2)):
            delta=Dia/2
    else:
        delta = (np.pi*Dia*(1-epsilon))/(2*(2*np.pi-thetaDry))

    return delta    

def computeSinglePhaseFrictionFactor(Dens,Visc,U,Dia):
    Re=(Dens*U*Dia)/Visc#Reynolds number computation
    if Re>=2300: 
        f0=0.079/(Re**0.25) #Turbulent flow(Blasius formula)
    else:
        f0=16/Re #Laminar flow (Fanning formula)
    return f0

def computeEpsilon(DensL,DensV,G,g,x,STen):
    return (x/DensV)*((1+0.12*(1-x))*((x/DensV)+((1-x)/DensL))+((1.18*(1-x)*(g*STen*(DensL-DensV))**0.25)/(G*DensL**0.5)))**-1

def computeThetaStrat(eps):
    return 2*np.pi-2*(np.pi*(1-eps)+(3*np.pi/2)**(1./3.)*(1-2*(1-eps)+(1-eps)**(1./3.)-eps**(1./3.))-(-1./120.)*(1-eps)*eps*(1-2*(1-eps))*(1+4*((1-eps)**2+eps**2)))#Biberg (1999)

def computePdropAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,Hajal=False):
    """ Hajal = Boolean value. If set to true, the film thickness will be computed using the Hajal et al.(2003b) expression which requires Gwavy and Gstrat."""
    thetaDry=0.0
    fiAnnular=computeFiAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,thetaDry,Hajal)
    #print(Hajal)
    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ug = (G/DensV)*(x/epsilon)
    return (4*fiAnnular*(L/Dia)*DensV*(Ug**2)/2)#Pressure drop

def computeFiAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,thetaDry,Hajal=False):
    """ Hajal = Boolean value. If set to true, the film thickness will be computed using the Hajal et al.(2003b) expression which requires Gwavy and Gstrat."""
    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ul = (G/DensL)*((1-x)/(1-epsilon))
    Ug = (G/DensV)*(x/epsilon)

    delta=computeDelta(DensL,DensV,ViscL,ViscV,STen,Dia,thetaDry,epsilon,Hajal)

    return(0.215*(delta/(Dia-2*delta))**0.515*((DensL-DensV)*g*Dia**2/STen)**(-0.298)*(ViscL/ViscV)**0.285*(DensL*Ul**2/(DensV*Ug**2))**(-0.044))


def computePdropSW(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,Gwavy,Gstrat,Hajal=False):
    
    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ug = (G/DensV)*(x/epsilon)
    f0g = computeSinglePhaseFrictionFactor(DensV,ViscV,Ug,Dia)#Single phase gas friction factor 
    xIA=fpm.xia(DensV,DensL,ViscL,ViscV)
    thetaStrat = computeThetaStrat(epsilon)
    if(Hajal):
        thetaDry =  x/xIA*((Gwavy-G)/(Gwavy-Gstrat))**(0.61)*thetaStrat#Wojtan et al. (2005b)
    else:    
        thetaDry =  (((Gwavy-G)/(Gwavy-Gstrat))**(0.61))*thetaStrat#Wojtan et al. (2005b)
    fiAnnular = computeFiAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,thetaDry,Hajal)#Annular flow pressure drop                                            #USE DIFFERENT THETADRY FOR PANNULAR ? IS GWAVY AT XIA ?
    thetaDryStar=thetaDry/(2*np.pi)

    fsw=thetaDryStar*f0g+(1-thetaDryStar)*fiAnnular#Friction factor
    return 4*fsw*(L/Dia)*(DensV*Ug**2)/2#Pressure drop

def computePdropSI(DensL,DensV,ViscL,ViscV,STen,x,epsilonIA,Dia,L,g,G):

    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ul0 = G/DensL
    f0l = computeSinglePhaseFrictionFactor(DensL,ViscL,Ul0,Dia)        
    DeltaP0l=4*f0l*(L/Dia)*(DensL*Ul0**2)/2 #Single phase liquid ONLY pressure drop
    DeltaPAnnular=computePdropAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G)
    DeltaPsi = DeltaP0l*((1-epsilon/epsilonIA)**0.25)+DeltaPAnnular*((epsilon/epsilonIA)**0.25)
    return DeltaPsi

def computePdropSSW(DensL,DensV,ViscL,ViscV,STen,x,epsilonIA,Dia,L,g,G,Gwavy,Gstrat):

    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ul0 = G/DensL       
    f0l = computeSinglePhaseFrictionFactor(DensL,ViscL,Ul0,Dia)
    DeltaP0l=4*f0l*(L/Dia)*(DensL*Ul0**2)/2 #Single phase liquid ONLY pressure drop
    DeltaPsw=computePdropSW(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,Gwavy,Gstrat,False)#Pressure drops in stratified wavy                                      # TO HAJAL, OR NOT ? THAT'S THE MOTHERFUCKING QUESTION
    DeltaPslugSW = (DeltaP0l*(1-epsilon/epsilonIA)**0.25)+(DeltaPsw*(epsilon/epsilonIA)**0.25)                                                               #HAJAL complex. For G=70, HAJAL ONLY SEEMS TO ENHANCE THE JUMP...
    return DeltaPslugSW

def computePdropMist(DensL,DensV,ViscL,ViscV,x,G,Dia,L):
     
    epsilonH=1/(1+((1-x)/x)*(DensV/DensL))#Homogeneous void fraction
    DensH=DensL*(1-epsilonH)+DensV*epsilonH#Homegeneous density
    ViscH=x*ViscV+(1-x)*ViscL#Homogeneous viscosity
    Rem=G*Dia/ViscH
    fm=0.079/(Rem**(1/4))   
    DeltaPm=2*fm*(L/Dia)*(G**2)/DensH
    return DeltaPm

def computePdropDryout(DensL,DensV,ViscL,ViscV,STen,x,xdi,xde,G,Dia,L,g,Gwavy,Gstrat):

    DeltaPm=computePdropMist(DensL,DensV,ViscL,ViscV,xde,G,Dia,L)
    WeFrl=(g*Dia**2*DensL)/STen
    A=np.pi*Dia**2/4
    
    if(xde<=xdi):
        return 0 #Dryout zone doesn't exist
    if(G>fpm.wavy(xdi,A,Dia,DensL,DensV,g,WeFrl,G,STen)): #Flow comes from annular  #NEED TO COMPUTE THIS GWAVY @ xdi                                        # GOOD CONDITIONS ?
           DeltaPtp= computePdropAnnular(DensL,DensV,ViscL,ViscV,STen,xdi,Dia,L,g,G)
    else:#Flow comes from stratified wavy
           DeltaPtp=computePdropSW(DensL,DensV,ViscL,ViscV,STen,xdi,Dia,L,g,G,Gwavy,Gstrat)   
    return (DeltaPtp-((x-xdi)/(xde-xdi))*(DeltaPtp-DeltaPm))
           

def computePdropStratified(DensL,DensV,ViscL,ViscV,STen,x,xIA,epsilonIA,G,Dia,L,g):
    epsilon = computeEpsilon(DensL,DensV,G,g,x,STen)
    Ug = (G/DensV)*(x/epsilon)
    Ul = (G/DensL)*((1-x)/(1-epsilon))
    thetaStrat = computeThetaStrat(epsilon)
    f0g=computeSinglePhaseFrictionFactor(DensV,ViscV,Ug,Dia)
    thetaDry=thetaStrat
    fa = computeFiAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,thetaDry)                                                                                    # ARENT WE USING FIANNULAR OUT OF ITS WORKING BOUNDARIES ???
         
    thetaStratStar=thetaStrat/(2*np.pi)
    fs= thetaStratStar*f0g+(1-thetaStratStar)*fa
    
    DeltaPstrat1=4*fs*L/Dia*DensV*(Ug**2)/2
    if(x>=xIA):
        DeltaPstrat=DeltaPstrat1
    else:
        Ul0 = G/DensL
        f0l = computeSinglePhaseFrictionFactor(DensL,ViscL,Ul0,Dia)
        DeltaP0l=4*f0l*(L/Dia)*(DensL*Ul0**2)/2 #Single phase liquid pressure drop
        DeltaPstrat=DeltaP0l*(1-epsilon/epsilonIA)**(1/4)+DeltaPstrat1*(epsilon/epsilonIA)**(1/4)
    return DeltaPstrat

def computePgrad(G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L):
    g=9.81
    hlv=HV-HL
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25
    vaporQualities=list(map(lambda x1: x1/500.0, range(1, 500)))
    Pgradient=[]
    for x in vaporQualities:
        pattern = fpm.computePattern(x,G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q)
        Pdrop = computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,qcrit,q)
        Pgrad = Pdrop/L                                                                                                                                                #PGRADIENT OR PDROP ???
        Pgradient.append(Pgrad)

    return(Pgradient)
  #  fpm.plotPattern(x,'fluid',300.,G,Dia,DensV,DensL,HV,HL,STen,ViscV,ViscL,q)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax2 = ax.twinx()
    #ax2.plot(vaporQualities, Pgradient, 'g')
    #ax.set_ylabel(r"MassFlux")
    #ax2.set_ylabel(r"Pgrad")

   #
    #plt.show()
def computePgradOne(G,x,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L):
    g=9.81
    hlv=HV-HL
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25
    pattern = fpm.computePattern(x,G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q)
    Pdrop = computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,qcrit,q)
    Pgrad = Pdrop/L                                                                                                                                                #PGRADIENT OR PDROP ???
    return int(Pgrad)

def computePdrop(pattern,DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,qcrit,q):
    #Logic
    WeFrl=(g*Dia**2*DensL)/STen
    A=np.pi*Dia**2/4
    xIA=fpm.xia(DensV,DensL,ViscL,ViscV)
    if(pattern==1):
        DeltaP=computePdropAnnular(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G)
        #print('annular')
    elif(pattern==6 or pattern==7):
        epsilonIA=fpm.epsi(xIA,DensV,DensL,STen,g,G)
        DeltaP=computePdropSI(DensL,DensV,ViscL,ViscV,STen,x,epsilonIA,Dia,L,g,G)
        #print('SI')
    elif(pattern==2):
        Gwavy=fpm.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)
        Gstrat=fpm.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)
        DeltaP=computePdropSW(DensL,DensV,ViscL,ViscV,STen,x,Dia,L,g,G,Gwavy,Gstrat)
        #print('sw')
    elif(pattern==8):       
        Gwavy=fpm.wavy(xIA,A,Dia,DensL,DensV,g,WeFrl,G,STen)
        Gstrat=fpm.strat(xIA,A,Dia,DensV,DensL,ViscL,g,STen,G)
        epsilonIA=fpm.epsi(xIA,DensV,DensL,STen,g,G)
        DeltaP=computePdropSSW(DensL,DensV,ViscL,ViscV,STen,x,epsilonIA,Dia,L,g,G,Gwavy,Gstrat)
        #print('ssw')
    elif(pattern==4):
        DeltaP=computePdropMist(DensL,DensV,ViscL,ViscV,x,G,Dia,L)
        #print('mist')
    elif(pattern==5):
        xdi1=fpm.xdi(G,Dia,DensV,DensL,STen,g,q,qcrit)
        xde1=fpm.xde(G,Dia,STen,DensV,DensL,g,qcrit,q)
        Gwavy=fpm.wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)
        Gstrat=fpm.strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)
        DeltaP=computePdropDryout(DensL,DensV,ViscL,ViscV,STen,x,xdi1,xde1,G,Dia,L,g,Gwavy,Gstrat)
        #print('dryout')
    else:# pattern ==3
        epsilonIA=fpm.epsi(xIA,DensV,DensL,STen,g,G)
        DeltaP=computePdropStratified(DensL,DensV,ViscL,ViscV,STen,x,xIA,epsilonIA,G,Dia,L,g)
        #print('s')
#    print ('DeltaP =',DeltaP)
#    print('x=',x)
    return DeltaP   
