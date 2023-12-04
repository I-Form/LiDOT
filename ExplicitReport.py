# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2019 replay file
# Internal Version: 2018_09_24-19.41.51 157541
# Run by 16438722 on Mon May  3 12:27:22 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
import sys
file_path=sys.argv[8]
file_name=sys.argv[9]

from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=50, 
    height=50)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
import os
os.chdir(file_path)
o1 = session.openOdb(
    name=file_path + file_name)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)

#Reaction Force
odb = session.odbs[file_path + file_name]
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF3'), )), ), nodeSets=("TOP", ))
x0 = session.xyDataObjects['RF:RF3 PI: RIGID_PLATE_2 N: 1']
session.xyReportOptions.setValues(numberFormat=SCIENTIFIC)
session.writeXYReport(fileName='RFReport.rpt', xyData=(x0, ))

#Internal and kinetic energy
odb = session.odbs[file_path + file_name]
xy0 = session.XYDataFromHistory(name='ALLIE Whole Model-1', odb=odb, 
    outputVariableName='Internal energy: ALLIE for Whole Model', steps=(
    'Step-1', ), __linkedVpName__='Viewport: 1')
c0 = session.Curve(xyData=xy0)
xy1 = session.XYDataFromHistory(name='ALLKE Whole Model-1', odb=odb, 
    outputVariableName='Kinetic energy: ALLKE for Whole Model', steps=(
    'Step-1', ), __linkedVpName__='Viewport: 1')

x0 = session.xyDataObjects['ALLIE Whole Model-1']
x1 = session.xyDataObjects['ALLKE Whole Model-1']
session.writeXYReport(fileName='Energy.rpt', xyData=(x0, x1))

#Element volume
# odb = session.odbs['C:/Users/16438722/Documents/GitHub/Lattice/Matlab/abq_out/%s.odb'% file_name]
# session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES)
# session.writeFieldReport(fileName='ElVol.csv', append=ON, 
    # sortItem='Element Label', odb=odb, step=0, frame=50, 
    # outputPosition=WHOLE_ELEMENT, variable=(('EVOL', WHOLE_ELEMENT), ), 
    # stepFrame=ALL)
    
# odb = session.odbs['C:/Users/16438722/Documents/GitHub/Lattice/Matlab/abq_out/%s.odb'% file_name]
# nf = NumberFormat(numDigits=8, precision=0, format=ENGINEERING)
# session.fieldReportOptions.setValues(numberFormat=nf)
# session.writeFieldReport(fileName='ElVol.rpt', append=ON, 
    # sortItem='Element Label', odb=odb, step=0, frame=50, 
    # outputPosition=WHOLE_ELEMENT, variable=(('EVOL', WHOLE_ELEMENT), ), 
    # stepFrame=ALL)

