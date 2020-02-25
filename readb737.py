# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:27:25 2020

@author: Raven
"""
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


################ INPUT

#opening input file
inp=open('B737.inp','r')
#reading
inptext=inp.read()
inp.close()
#beginphrase to state where data starts
beginphraseinp='name=B737\n*Node\n'
#endphrase to state where data ends
endphraseinp='\n*Element, type=S4R'
#stripping the data to get coordinates with node number
stripinp=inptext[inptext.find(beginphraseinp) + len(beginphraseinp):].split(endphraseinp,1)[0]
#converting to pd dataframe
dfinp = pd.DataFrame([x.split(',') for x in stripinp.split('\n')],columns=['node','x','y','z']).astype(float)
#setting node as index
dfinp=dfinp.set_index('node')

################ INPUT END

#plotting yz plane of aileron
#df.plot(kind='scatter',x='z',y='y',color='red')
#plt.show()


################ OUTPUT

#opening output file
outp = open('B737output','r')
#reading
outptext=outp.read()
outp.close()
#beginphrase to state where data starts
beginphraseoutp='----------------------------------------------------------------------------------------\n'
#endphrase to state where data ends
endphraseoutp='\n\n\n  Minimum'
#stripping the data to get coordinates with node number
stripoutp=outptext[outptext.find(beginphraseoutp) + len(beginphraseoutp):].split(endphraseoutp,1)[0]
#stripping data
stripoutp=stripoutp.lstrip()
stripoutp=stripoutp.replace('               ',',')
stripoutp=stripoutp.replace('            ',',')
stripoutp=stripoutp.replace('     ',',')
stripoutp=stripoutp.replace('    ',',')
stripoutp=stripoutp.replace('   ',',')
#converting to pd dataframe
dfoutp = pd.DataFrame([x.split(',') for x in stripoutp.split('\n,')],columns=['node','1','smises1','smises2','ss1','ss2'])[['node','smises1','smises2','ss1','ss2']].astype(float)


#setting node as index
dfoutp=dfoutp.set_index('node')

################ OUTPUT END

#left joining input with output on key
dfmerged=dfinp.join(dfoutp)
#calculating averages
dfmerged['smises'] = dfmerged[['smises1','smises2']].mean(axis=1)
dfmerged['ss'] = dfmerged[['ss1','ss2']].mean(axis=1)
dfmerged=dfmerged[['x','y','z','smises','ss']]

'''
threedee = plt.figure().gca(projection='3d')
threedee.scatter(df.index, df['z'], df['y'])
threedee.set_xlabel('x')
threedee.set_ylabel('y')
threedee.set_zlabel('z')
plt.show()
'''


