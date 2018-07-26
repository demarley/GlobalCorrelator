"""
Created:        17 March 2018
Last Updated:   20 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Read output data file from propagator
Data from:
 https://github.com/demarley/GlobalCorrelator/tree/tk-mu-linking/
                             pulsar_devkit/simple_examples/HLSIPs

each line in text file:
 in_eta,out_eta,in_phi,out_phi,in_rinv,inz0,
"""
import os
import sys
from math import fabs,log,sqrt


# plotting options
import matplotlib
matplotlib.use("PDF")
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='sans-serif')
import numpy as np

import matplotlib.pyplot as plt
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

# TeX labels for different features
featureNames = {
'deta':   r"$\Delta\eta$",
'dphi':   r"$\Delta\phi$",
'sinhEta':r"sinh($\eta$)",
'eta':    r"$\eta$",
'phi':    r"$\phi_{\text{0}}$",
'rinv':   r"R$_{\text{inv}}$",
'charge': r"Charge",
'z0':     r"$\text{z}_\text{0}$",
'propEta':r"Propagated $\eta$",
'propPhi':r"Propagated $\phi_{\text{0}}$"
}



def plot1D(data=None,label="",saveAs="",extra_text=""):
    """Make one-dimensional plot"""
    fig,ax = plt.subplots()

    plt.hist( data,bins=10,histtype='stepfilled' )

    plt.xlabel(label,position=(1,0),ha='right',fontsize=20)
    x_tick_labels = []
    if extra_text.lower()=="binary":
        x_tick_labels = [int(i) for i in ax.get_xticks()]
    else:
        x_tick_labels = ax.get_xticks()
    ax.set_xticklabels(x_tick_labels,fontsize=18)

    plt.ylabel("N Tracks",position=(0,1),ha='right',fontsize=20)
    ax.set_yticklabels([int(i) for i in ax.get_yticks()],fontsize=18)

    if extra_text:
        ax.text(0.97,0.97,extra_text,fontsize=18,ha="right",va="top",transform=ax.transAxes)

    plt.savefig(saveAs,format='pdf',dpi=300,bbox_inches='tight')
    plt.close()

    return


def plot2D(x_data=None,y_data=None,x_label="",y_label="",saveAs=""):
    """Make two-dimensional plot"""
    fig,ax = plt.subplots()

    plt.hist2d(x_data,y_data,bins=20,cmap='viridis')

    if "binary" in saveAs.lower():
        x_tick_labels = [int(i) for i in ax.get_xticks()]
    else:
        x_tick_labels = ax.get_xticks()

    plt.xlabel(x_label,position=(1,0),ha='right',fontsize=20)
    ax.set_xticklabels(x_tick_labels,fontsize=18)

    plt.ylabel(y_label,position=(0,1),ha='right',fontsize=20)
    ax.set_yticklabels(ax.get_yticks(),fontsize=18)

    plt.savefig(saveAs,format='pdf',dpi=300,bbox_inches='tight')
    plt.close()

    return


class Track(object):
    def __init__(self):
        self.index   = 0
        self.rinv    = 0
        self.charge  = 0
        self.simPt      = -999
        self.simPhi     = -999
        self.simEta     = -999
        self.swTrackPt      = -999
        self.swTrackPhi     = -999
        self.swTrackEta     = -999
        self.hwTrackPt      = -999
        self.hwTrackPhi     = -999
        self.hwTrackEta     = -999
        self.swPropTrackPt      = -999
        self.swPropTrackPhi     = -999
        self.swPropTrackEta     = -999

def bit2int(value):
    factor   = 1
    if value[0]=="1": 
        value = "0"+value[1:]
        factor = -1

    intvalue = factor*int(value,base=2)

    return intvalue


def file2list(filename):
    """Load text file and dump contents into a list"""
    listOfFiles = open( filename,'r').readlines()
    listOfFiles = [i.rstrip('\n') for i in listOfFiles if not i.startswith("#")]
    return listOfFiles



print " Read input files "

# Input Data
f1 = file2list("output_data.txt")

tracks = []
print " Read file "
for fw in f1:
    #print fw
    dataBinary = fw.split(',')
    t1 = Track()  # binary
    t1.simPhi = float(dataBinary[4][5:])
    t1.swTrackPhi = float(dataBinary[8][5:])
    t1.swTrackSector = float(dataBinary[11][8:])
    tracks.append(t1)

    continue

print " Make input plots "

simPhiData = [getattr(x,'simPhi') for x in tracks]
swTrackPhiData  = [getattr(x,'swTrackPhi') for x in tracks]
swTrackSectorData  = [getattr(x,'swTrackSector') for x in tracks]

swTrackPhiMult = swTrackPhiData + swTrackSectorData * (2/9.)*np.pi

plot2D(x_data=simPhiData,y_data=swTrackPhiMult,
       x_label='Sim Phi',y_label='Sw Track Phi',
       saveAs="simPhi_vs_SwTrackPhi.pdf")

# Plot some of the inputs
"""



    featureName = featureNames[feature]

    plot1D(data=floatData,label=featureName,extra_text="Floating-point",
           saveAs="tracks_float_{0}.pdf".format(feature))

    plot1D(data=binaryData,label=featureName,extra_text="Binary",
           saveAs="tracks_binary_{0}.pdf".format(feature))





# Output Data
print " Read output files "

f1 = file2list("firmware_prop.txt")
f2 = file2list("software_prop.txt")

tracks_fw = []
tracks_sw = []

for fw,sw in zip(f1,f2):

    data_fw = fw.split(',')
    data_sw = sw.split(',')


    # firmware
    t1 = Track()

    t1.eta     = float(data_fw[0])
    t1.propEta = float(data_fw[1])
    t1.phi     = float(data_fw[2])
    t1.propPhi = float(data_fw[3])
    t1.rinv    = float(data_fw[4])
    t1.z0      = float(data_fw[5])
    tracks_fw.append(t1)

    # software
    t2 = Track()

    t2.eta     = float(data_sw[0])
    t2.propEta = float(data_sw[1])
    t2.phi     = float(data_sw[2])
    t2.propPhi = float(data_sw[3])
    t2.rinv    = float(data_sw[4])
    t2.z0      = float(data_sw[5])
    tracks_sw.append(t2)


print " Make output plots "
# Plot some of the features
for feature in ['eta','phi','rinv','z0','propEta','propPhi']:
    print " > Plotting "+feature
    featureName = featureNames[feature]

    fwData = [getattr(x,feature) for x in tracks_fw]
    swData = [getattr(x,feature) for x in tracks_sw]

    plot2D(x_data=fwData,y_data=swData,
           x_label=featureName+"-fw",
           y_label=featureName+"-sw",
           saveAs="tracks_fwsw_prop_{0}.pdf".format(feature))


deltaEta = [ i.propEta-j.propEta for i,j in zip(tracks_sw,tracks_fw) ]
deltaPhi = [ i.propPhi-j.propPhi for i,j in zip(tracks_sw,tracks_fw) ]
inEta    = [ i.eta for i in tracks_sw ]
inPhi    = [ i.phi for i in tracks_sw ]
inRinv   = [ i.rinv for i in tracks_sw ]

# DeltaEta vs inEta
plot2D(x_data=inEta,y_data=deltaEta,
       x_label="Input {0}".format(featureNames['eta']),
       y_label="Propagated {0}".format(featureNames['deta']),
       saveAs="tracks_fwsw_prop_dEta-v-eta.pdf")

# PropEta vs InEta
propEta = [ i.propEta for i in tracks_sw ]
plot2D(x_data=inEta,y_data=propEta,
       x_label="Input {0}".format(featureNames['eta']),
       y_label="Propagated {0}".format(featureNames['eta']),
       saveAs="tracks_fwsw_inEta-v-propEta.pdf")

propPhi = [ i.propPhi for i in tracks_sw ]
plot2D(x_data=inPhi,y_data=propPhi,
       x_label="Input {0}".format(featureNames['phi']),
       y_label="Propagated {0}".format(featureNames['phi']),
       saveAs="tracks_fwsw_inPhi-v-propPhi.pdf")

# DeltaPhi vs inEta
plot2D(x_data=inEta,y_data=deltaPhi,
       x_label="Input {0}".format(featureNames['eta']),
       y_label="Propagated {0}".format(featureNames['dphi']),
       saveAs="tracks_fwsw_prop_dPhi-v-eta.pdf")

# DeltaPhi vs inPhi
plot2D(x_data=inPhi,y_data=deltaPhi,
       x_label="Input {0}".format(featureNames['phi']),
       y_label="Propagated {0}".format(featureNames['dphi']),
       saveAs="tracks_fwsw_prop_dPhi-v-phi.pdf")

# DeltaEta vs inRinv
plot2D(x_data=inRinv,y_data=deltaEta,
       x_label="Input {0}".format(featureNames['rinv']),
       y_label="Propagated {0}".format(featureNames['deta']),
       saveAs="tracks_fwsw_prop_dEta-v-rinv.pdf")

# DeltaPhi vs inRinv
plot2D(x_data=inRinv,y_data=deltaPhi,
       x_label="Input {0}".format(featureNames['rinv']),
       y_label="Propagated {0}".format(featureNames['dphi']),
       saveAs="tracks_fwsw_prop_dPhi-v-rinv.pdf")
"""


## THE END ##
