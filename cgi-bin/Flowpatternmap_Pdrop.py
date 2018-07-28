# -*- coding: utf-8 -*-
"""
Flow pattern map code

Made by Robert Nederlof
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import Straight_Pdrop as st
import Ubend_Pdrop as ub
import geo
import sat



def computePattern(x,G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q):
    g=9.81                      #Gravitational acceleration [m/s2]
    A=np.pi*Dia**2/4            #Cross-sectional area [m2]
    hlv=HV-HL                   #Latent heat of vaporisation [J/kg]
    WeFrl=(g*Dia**2*DensL)/STen #(Wel/Frl)
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25   #Critical heat flux

    #Vertical and horizontal lines
    xia=((0.34**(1./0.875)*(DensV/DensL)**(-1./1.75)*(ViscL/ViscV)**(-1./7.))+1)**-1
    slug=wavy(xia,A,Dia,DensL,DensV,g,WeFrl,G,STen)

    #Values of parameter z

    #Annular=1
    #Stratified-Wavy=2
    #Stratified=3
    #Mist=4
    #Dryout=5
    #Intermittent=6
    #Slug=7
    #Slug+Stratified-Wavy

    if xia<=x<=xdi(G,Dia,DensV,DensL,STen,g,q,qcrit):
        if wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen)<=G<mist(x,DensV,Dia,DensL,STen,q,qcrit,g):
            z=1
        elif strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)<=G<wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen):
            z=2
        elif G<=strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G):
            z=3
        else:
            z=4

    elif x>xdi(G,Dia,DensV,DensL,STen,g,q,qcrit) and G>=strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G):
        if G<mist(x,DensV,Dia,DensL,STen,q,qcrit,g):
            z=5
        else:
            z=4

    elif x<xia:
        if G>=wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen):
            z=6
        elif slug<=G<wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen):
            z=7
        elif strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G)<=G<slug:
            z=8
        else:
            z=3
    else:
        z=3

    return z

def plotPattern(x,fluid,Tsat,G,Dia,DensV,DensL,HV,HL,STen,ViscV,ViscL,q,D,Geometry):
    Gdryout=[]
    Gmist=[]
    Gwavy=[]
    Gstrat=[]
    GG=[0,10000]
    xd=[]
    dx=0.001                    #Step size
    xx=list(map(lambda x1: x1/500.0, range(1, 500)))     #Vapor quality's
    g=9.81                      #Gravitational acceleration [m/s2]
    A=np.pi*Dia**2/4            #Cross-sectional area [m2]
    hlv=HV-HL                   #Latent heat of vaporisation [J/kg]
    WeFrl=(g*Dia**2*DensL)/STen #(Wel/Frl)
    qcrit=0.131*DensV**0.5*hlv*(g*(DensL-DensV)*STen)**0.25   #Critical heat flux

    #Vertical line
    xia=((0.34**(1./0.875)*(DensV/DensL)**(-1./1.75)*(ViscL/ViscV)**(-1./7.))+1)**-1

    #Creating lines
    for i in range(len(xx)):
        X=xx[i]

        #Startified line
        a=strat(X,A,Dia,DensV,DensL,ViscL,g,STen,G)
        if X>=xia:
            Gstrat.append(a)
        else:
            None

        #Mist line
        b=mist(X,DensV,Dia,DensL,STen,q,qcrit,g)
        Gmist.append(b)

        #Dryout line
        c=dry(X,Dia,DensV,DensL,STen,g,q,qcrit)
        if a<=c<b and X>xia:
            Gdryout.append(c)
            xd.append(X)
        else:
            None

        #Wavy line
        d=wavy(X,A,Dia,DensL,DensV,g,WeFrl,G,STen)
        if d<=c or X<xia:
            Gwavy.append(d)
        else:
            None

    #Straight lines
    XIA=[xia,xia]
    GG[0]=Gstrat[0]
    XSW=[dx,xia]
    GSW=[(Gstrat[0]),(Gstrat[0])]
    slug=wavy(xia,A,Dia,DensL,DensV,g,WeFrl,G,STen)
    GS=[slug,slug]
    #Defining letter D in plot
    if xde(500,Dia,STen,DensV,DensL,g,qcrit,q)-xdi(500,Dia,DensV,DensL,STen,g,q,qcrit)<=0.075:
        h=150
        if q<5000:
            l=0.01
        else:
            l=0.02
    elif 0.075<=(xde(500,Dia,STen,DensV,DensL,g,qcrit,q)-xdi(500,Dia,DensV,DensL,STen,g,q,qcrit))<=0.15:
        h=400
        l=0.035
    else:
        h=500
        l=0.075

    if G<700:
        p=700
    else:
        p=G+50

    if q>=400:
        k=1
    else:
        k=0

    #Plotting
    ax=plt.subplot(111)
    plt.plot(xd,Gdryout,'b')
    plt.plot(xx,Gmist,'b')
    plt.plot(xx[:(len(Gwavy))],Gwavy,'b')
    plt.plot(xx[(len(xx)-len(Gstrat)):],Gstrat,'b')
    plt.plot(XIA,GG,'b')
    plt.plot(XSW,GSW,'b')
    plt.plot(XSW,GS,'b')
    plt.plot(x,G,color='r',marker='o',markersize=15)
    plt.plot([0,1],[G,G],color='r',linestyle='--')
    #Edditing
    plt.axis([0,1,0,p])
    plt.text((xia+xdi(400,Dia,DensV,DensL,STen,g,q,qcrit))/2,400,'ANR',fontsize=20)
    plt.text(0.12,(Gstrat[0]/3),'STD',fontsize=10)
    plt.text((xia/3),(Gstrat[0]+75),'SSW',fontsize=15)
    plt.text((xia/1.5),wavy((xia/2.5),A,Dia,DensL,DensV,g,WeFrl,G,STen),'IMT',fontsize=20)
    plt.text((xia/4),((wavy((xia/4.),A,Dia,DensL,DensV,g,WeFrl,G,STen)-slug)/3+slug),'SLG',fontsize=15)
    plt.text((xia+0.2),((wavy((xia+0.15),A,Dia,DensL,DensV,g,WeFrl,G,STen)+strat((xia+0.15),A,Dia,DensV,DensL,ViscL,g,STen,G))/2),'SW',fontsize=15)
    if k==1:
        plt.text((xdi(h,Dia,DensV,DensL,STen,g,q,qcrit)+l),(h),'DYT',fontsize=15)
        plt.text(((1+xde(600,Dia,STen,DensV,DensL,g,qcrit,q))/2),(600),'MST',fontsize=15)
    plt.xlabel('Vapor Quality [-]')
    plt.ylabel('Mass Flux [kg/s/m$^2$]', color="b")
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.xaxis.grid(True,'minor')
    ax.yaxis.grid(True,'minor')
    ax.xaxis.grid(True,'major')
    ax.yaxis.grid(True,'major')
    ax.set_xlim(0,1)
    MassF=[0,p]
    plt.ylim(0, 1.05*p)
    gmax= np.ceil(int(p/100))*100
    plt.yticks(range(0, int(1.05*gmax), int(np.ceil(gmax/10))))
    plt.tick_params(axis="y", labelcolor="b", pad=8)
    L=2.0
    if Geometry ==1:
        Pgradient=st.computePgrad(G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L)
        pd=st.computePgradOne(G,x,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L)
        plt.title('%s, G=%s kg/s/m$^2$, Tsat=%s$^\circ$C, d=%s mm, q=%s kW/m$^2$'%(fluid,G,Tsat-273,(Dia*1000),(q/1000.)),fontsize=14)
    else:
        Pgradient=ub.computePgrad(G,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L,D,Geometry)
        pd=ub.computePgradOne(G,x,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,q,L,D,Geometry)
        plt.title('%s, G=%s kg/s/m$^2$, Tsat=%s$^\circ$C, d=%s m, D=%s m, q=%s kW/m$^2$'%(fluid,G,Tsat-273,Dia,D,(q/1000.)),fontsize=12)
    ax2 = ax.twinx()
    ax2.set_xlim(0,1)
    pmax=np.ceil(int(max(Pgradient)/100))*100
    plt.plot(x,pd,color='r',marker='o',markersize=15)
    ax2.set_ylabel('Pressure Gradient [Pa/m]', color="g")
    plt.ylim(0, int(1.05*pmax))
    plt.yticks(range(0, int(1.05*pmax), int(np.ceil(pmax/10))))
    plt.tick_params(axis="y", labelcolor="g", pad=8)
    ax2.plot(xx, Pgradient, 'g', linewidth=2)



    plt.show()

#Functions
def epsi(x,DensV,DensL,STen,g,G):
    eps=(x/DensV)*((1+0.12*(1-x))*((x/DensV)+((1-x)/DensL))+((1.18*(1-x)*(g*STen*(DensL-DensV))**0.25)/(G*DensL**0.5)))**-1
    return eps

def thet(eps):
    theta=2*np.pi-2*(np.pi*(1-eps)+(3*np.pi/2)**(1./3.)*(1-2*(1-eps)+(1-eps)**(1./3.)-eps**(1./3.))-(-1./120.)*(1-eps)*eps*(1-2*(1-eps))*(1+4*((1-eps)**2+eps**2)))
    return theta

def strat(x,A,Dia,DensV,DensL,ViscL,g,STen,G):
    eps=epsi(x,DensV,DensL,STen,g,G)
    Avd=A*eps/(Dia**2)
    Ald=A*(1-eps)/(Dia**2)
    a=((226.3**2*Ald*Avd**2*DensV*(DensL-DensV)*ViscL*g)/(x**2*(1-x)*np.pi**3))**(1./3.)
    return a

def mist(x,DensV,Dia,DensL,STen,q,qcrit,g):
    b=((1./0.0058)*(np.log(0.61/x)+0.57)*(Dia/(DensV*STen))**-0.38*(1./(g*Dia*DensV*(DensL-DensV)))**-0.15*(DensV/DensL)**0.09*(q/qcrit)**-0.27)**0.943
    return b

def dry(x,Dia,DensV,DensL,STen,g,q,qcrit):
    c=((1./0.235)*(np.log(0.58/x)+0.52)*(Dia/(DensV*STen))**-0.17*(1./(g*Dia*DensV*(DensL-DensV)))**-0.37*(DensV/DensL)**-0.25*(q/qcrit)**-0.70)**0.926
    return c

def wavy(x,A,Dia,DensL,DensV,g,WeFrl,G,STen):
    eps=epsi(x,DensV,DensL,STen,g,G)
    Avd=A*eps/(Dia**2)
    theta=thet(eps)
    hld=0.5*(1-np.cos((2*np.pi-theta)/2))
    d=(((16*(Avd**3)*g*Dia*DensL*DensV)/(x**2*np.pi**2*(1-(2*hld-1)**2)**0.5))*((np.pi**2/(25*hld**2))*(WeFrl)**-1+1))**0.5+50
    return d

def xdi(G,Dia,DensV,DensL,STen,g,q,qcrit):
    Wev=G**2*Dia/(DensV*STen)
    Frv=(G**2)/(DensV*(DensL-DensV)*g*Dia)
    di=0.58*np.exp(0.52-0.235*Wev**0.17*Frv**0.37*(DensV/DensL)**0.25*(q/qcrit)**0.7)
    if di>1:
        di=1
    return di

def xde(G,Dia,STen,DensV,DensL,g,qcrit,q):
    Wev=G**2*Dia/(DensV*STen)
    Frv=(G**2)/(DensV*(DensL-DensV)*g*Dia)
    de=0.61*np.exp(0.57-0.0058*Wev**0.38*Frv**0.15*(DensV/DensL)**-0.09*(q/qcrit)**0.27)
    if de>1:
        de=1
    return de

def xia(DensV,DensL,ViscL,ViscV):
    xia=((0.34**(1./0.875)*(DensV/DensL)**(-1./1.75)*(ViscL/ViscV)**(-1./7.))+1)**-1
    return xia

