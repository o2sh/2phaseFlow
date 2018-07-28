# the first two lines are required for headless linux server
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sat
import os
import traceback
import geo
import Flowpatternmap_Pdrop as flow
from matplotlib.ticker import MultipleLocator
import Ubend_Pdrop as ub
import Straight_Pdrop as st

# big function that compute the results in main
def GetImage(Param):
    #if(not(os.path.isfile(Param['IStr']+'.png'))): #do nothing if image exist
    if (True):  #no image caching
        try:
            Render(Param)
        except Exception as Err:  #display error message
            Param['ImgPath'] = "fatal_error"
            ErrStr = "RENDERING ERROR: \n{}\n\n STACK TRACE:\n{}".format(Err, traceback.format_exc())
            ErrStr = ErrStr.replace("\n", "<br />")
            Param['ErrStr'] = ErrStr


#Warning:Param is the same object, not a copy, of the argument in the calling function

# we get the parameters and make calculations then plot
def Render(Param):
    BendDiameter = Param['BendDiameter']
    MassFlux = Param['MassFlux']
    VaporQuality = Param['VaporQuality']
    Dia = Param['Dia']
    QFlux = Param['QFlux']
    Tsat = Param['Tsat'] +273

    tab = sat.GetTab(Param['Fluid'])  # extract the fluid characteristics from sat

    tab2 = geo.GetTab(Param['Geometry'])


    Lim = [np.amin(tab['T']), np.amax(tab['T'])]
    if Tsat < Lim[0] or Tsat > Lim[1]:
        Param['ImgPath'] = "fatal_error"
        Param['ErrStr'] = "ERROR: Temperature, {} K, is outside the range {} to {} K acceptable for {}".format(
            Tsat, Lim[0], Lim[1], Param['Fluid'])
        return

    Lim = [1.0e-7, 2.0e-2]
    if Dia < Lim[0] or Dia > Lim[1]:
        Param['ImgPath'] = "fatal_error"
        Param['ErrStr'] = "ERROR: Diameter, {} m, is outside the range {} to {} m".format(
            Dia, Lim[0], Lim[1])
        return


    Lim = [0.0025, 1000]
    if BendDiameter < 1.1*Dia or BendDiameter > 1000:
        Param['ImgPath'] = "fatal_error"
        Param['ErrStr'] = "ERROR: Bend Diameter, {} m, is outside the range {} to {} m".format(
            BendDiameter, 1.1*Dia, 1000)
        return


    Lim = [100, 1200]
    if MassFlux < Lim[0] or MassFlux > Lim[1]:
        Param['ImgPath'] = "fatal_error"
        Param['ErrStr'] = "ERROR: Mass Flux, {}  kg/s/m², is outside the range {} to {}  kg/s/m²".format(
            MassFlux, Lim[0], Lim[1])
        return

    ViscL = np.interp(Tsat, tab['T'], tab['ViscL'])
    DensL = np.interp(Tsat, tab['T'], tab['DensL'])
    ViscV = np.interp(Tsat, tab['T'], tab['ViscV'])
    DensV = np.interp(Tsat, tab['T'], tab['DensV'])
    STen = np.interp(Tsat, tab['T'], tab['STen'])
    HL = np.interp(Tsat, tab['T'], tab['HL'])
    HV = np.interp(Tsat, tab['T'], tab['HV'])
    Geometry = tab2['Geo']

    L=2.0
    flow.plotPattern(VaporQuality,Param['Fluid'],Tsat,MassFlux,Dia,DensV,DensL,HV,HL,STen,ViscV,ViscL,QFlux,BendDiameter,Geometry)
    if Geometry ==1:
        Param['P_drop']=st.computePgradOne(MassFlux,VaporQuality,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,QFlux,L)
    else:
        Param['P_drop']=ub.computePgradOne(MassFlux,VaporQuality,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,QFlux,L,BendDiameter,Geometry)

 #   st.computePgrad(MassFlux,Dia,DensV,DensL,STen,HV,HL,ViscV,ViscL,QFlux,L)
    plt.savefig(Param['ImgPath'] + '.png', bbox_inches='tight')  # you put the flow pattern in param['imgPath']
    plt.savefig(Param['ImgPath'] + '.svg', bbox_inches='tight')






