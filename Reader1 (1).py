import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import * 

h = 0.205
r = h/2
rmsy =[]
rmsz = []
rmst = []

def readStress(head, coor):
    ''' Reads data from the CSV file, gives out file in format: index, x, y, z, avg von misses stress, avg shear stress'''

    rptBend1 = pd.read_csv('B737.rpt', header=head, nrows=5778,
                           names=['index', 'intPt', 'Mises@1', 'Mises@2', 'SS12@1', 'SS12@2'], delim_whitespace=True,
                           index_col='index')
    rptBend2 = pd.read_csv('B737.rpt', header=head + 5789, nrows=856,
                           names=['index', 'intPt', 'Mises@1', 'Mises@2', 'SS12@1', 'SS12@2'], delim_whitespace=True,
                           index_col='index')
    rptBend = rptBend1.append(rptBend2)
    rptBend['avgMises'] = rptBend[['Mises@1', 'Mises@2']].mean(axis=1)
    rptBend['avgShear'] = rptBend[['SS12@1', 'SS12@2']].mean(axis=1)
    rptBend = rptBend.drop(['intPt', 'Mises@1', 'Mises@2', 'SS12@1', 'SS12@2'], axis=1)
    case1Stress = pd.concat([coor, rptBend.reindex(coor.index)], axis=1)



    return case1Stress

def readDeflection(head):
    file = pd.read_csv('B737.rpt', header=head, nrows=6588,
                           names=['index', 'dMag', 'dX', 'dY', 'dZ'], delim_whitespace=True,
                           index_col='index')
    
   # print(file)
    
    
    file['dMag'] /=1000
    file['dX']   /=1000
    file['dY']   /=1000
    file['dZ']   /=1000


    return file

def angle(y1, dy1 , y2 , dy2, z1, dz1 , z2, dz2 ):
    y = ((y1 + dy1) - (y2 + dy2))
    z = ((z1 + dz1 ) - (z2 + dz2))
    angle = - np.arctan2(y, z)
  #  angle = angle* 180/np.pi
    return angle

# ---------------------------------------------------------------------------------------------------------
#Loading up the data
#SS12 is a shear stress

coor = pd.read_csv('B737.inp', header=None, nrows=6588, names = ['index', 'x' , 'y' , 'z'], comment='*', index_col = 'index') #Reads initial file


# Case 1 ------------------------------------------------------------------------------------
# Just bending max at
case1Stress= readStress(13, coor)
case1def = readDeflection(20010)
case1 = pd.concat([case1Stress, case1def.reindex(case1Stress.index)], axis=1)

case1['x'] /=1000
case1['y'] /=1000
case1['z'] /=1000
case1['z'] -= r 

case1['avgMises'] *=10**9
case1['avgShear'] *=10**9


# Output is:    'index'   'x'   'y'   'z'   'avgMisses'   'avgS12'   'dMag'   'dX'   'dY'   'dZ'

#print(case1)


# Case 2 -------------------------------------------------------------------------------------
# Jammed and Bent


case2Stress = readStress(6679, coor)
case2def = readDeflection(26641)
case2 = pd.concat([case2Stress, case2def.reindex(case2Stress.index)], axis=1)

case2['x'] /=1000
case2['y'] /=1000
case2['z'] /=1000
case2['z'] -= r 

case2['avgMises'] *=10**9
case2['avgShear'] *=10**9


# Case 3 ------------------------------------------------------------------------------------
# Jammed and Straight


case3Stress = readStress( 13345, coor)
case3def = readDeflection(33360-88)
case3 = pd.concat([case3Stress, case3def.reindex(case3Stress.index)], axis=1)

case3['x'] /=1000
case3['y'] /=1000
case3['z'] /=1000
case3['z'] -= r 

case3['avgMises'] *=10**9
case3['avgShear'] *=10**9

case1 = case2



#Choosing hinge line --------------------------------------------------------------------------------------

hingeline1 = case1.loc[(case1['z'] == 0) & ( case1['y'] == 0) ]
hingeline1 = hingeline1.sort_values('x')



# Calculating the angular deflection --------------------------------------------------------------------

# Finding the leading and trailing edge
leadEdge1 = case1.loc[(case1['z'] == 0 ) & ( case1['y'] == 0) ]
leadEdge1 = leadEdge1.sort_values('x') # Sorting it with ascending x

trailEdge1 = case1.loc[(case1['z'] == -0.605 ) & ( -0.000001<= case1['y'] ) & (case1['y'] <= 0.0000001) ]
trailEdge1 = trailEdge1.sort_values('x') # Sorting it with ascending x
#Calculating tangent



trailEdge1.index = leadEdge1.index

#trailEdge1['z+dz'] = trailEdge1['z'] + trailEdge1['dZ']
#leadEdge1['z+dz'] = leadEdge1['z'] + leadEdge1['dZ']
#leadEdge1['test'] = leadEdge1['dY']

trailEdge1['angle'] = angle(leadEdge1['y'],leadEdge1['dY'], trailEdge1['y'], trailEdge1['dY'] , leadEdge1['z'],leadEdge1['dZ'], trailEdge1['z'], trailEdge1['dZ'])
#print(trailEdge1)

# Adding the dx dy dz to x y z coordinates
k = case1.loc[case1['x'] == 1.43620837]

'''
plt.plot(k['z'],k['y'])
plt.xlabel('z')
plt.ylabel('y')
plt.show()


plt.scatter(k['z'],k['y'],c=k['avgMises'],cmap = 'jet')
plt.xlabel('z[m]', fontsize=16)
plt.ylabel('y[m]', fontsize=16)
#plt.legend()
c=plt.colorbar()
c.set_label('Total stress [N/m^2]')
plt.title('Von Mises stress distribution')
#plt.savefig("flow.png")
plt.show()
'''


trailEdge1.to_csv('trailEdge1.txt', columns = 'x')
leadEdge1.to_csv('leadEdge1.txt', columns = 'x')

# 2d plots -------------------------------------------------------
#plt.plot( hingeline1['x'] , hingeline1['dZ'])
#plt.show()

'''
#Case 2
xver = [0., 0.29566667, 0.59133333, 0.887   ,   1.18266667 ,1.47833333,
 1.774    ,  2.06966667 ,2.36533333, 2.661     ]

Dy = [ 1.22335503e-02,  8.62534510e-03,  5.18045648e-03,  2.22466131e-03,
        6.90571122e-05, -2.49968619e-04,  1.93018819e-03,  6.08741941e-03,
        1.15563208e-02,  1.76115727e-02]

Dz  = [-6.82770018e-03, -4.40941437e-03, -2.18506230e-03, -5.71207911e-04,
        9.56402039e-06, -5.87890869e-04, -2.05030643e-03, -4.13553902e-03,
       -6.61702541e-03, -9.27042406e-03]

Twist = [-0.0093243 , -0.00928255, -0.00915694, -0.00894666, -0.01037132,
       -0.0127212 , -0.01298129, -0.0132283 , -0.01341143, -0.01350802]
'''

xver = [0., 0.29566667, 0.59133333, 0.887   ,   1.18266667 ,1.47833333,
 1.774    ,  2.06966667 ,2.36533333, 2.661     ]

Dy = [-4.38129032e-04,  2.15365260e-04,  7.31888459e-04,  7.92551526e-04,
        4.32721625e-05, -9.99157927e-04, -1.42946386e-03, -1.26891462e-03,
       -6.79958046e-04,  1.11861897e-04] 

Dz = [-8.90683623e-05,  6.28840923e-05,  1.80589462e-04,  1.90303641e-04,
        2.31705334e-05, -1.89644098e-04, -2.63921333e-04, -2.23853290e-04,
       -1.10338424e-04,  3.51946587e-05]

Twist = [-0.00908039, -0.00907283, -0.00902937, -0.00889949, -0.01040362,
       -0.01270926, -0.01290959, -0.01309376, -0.01321462, -0.01326411]

i = 0 
for x in xver:
    params = hingeline1.iloc[(hingeline1['x']-x).abs().argsort()[:1]]

    Dyv = params.iloc[0]['dY']
    #print(Dyv)
    #print(Dy[i])
    p = (Dyv-Dy[i])**2
    #print(p)
    i = i+1
    rmsy.append(p)
    
MSEy = (sum(rmsy)/len(rmsy))**(1/2)
    
print("The MSE for the dY is equal to",MSEy)

i = 0 
for x in xver:
    params = hingeline1.iloc[(hingeline1['x']-x).abs().argsort()[:1]]

    Dzv = params.iloc[0]['dZ']
    #print(Dzv)
    #print(Dz[i])
    p = (Dzv-Dz[i])**2
    #print(p)
    i = i+1
    rmsz.append(p)
    
MSEz = (sum(rmsz)/len(rmsz))**(1/2)
    
print("The MSE for the dZ is equal to",MSEz)

i = 0 
for x in xver:
    params = trailEdge1.iloc[(trailEdge1['x']-x).abs().argsort()[:1]]

    Dtwist = params.iloc[0]['angle']
    #print(Dzv)
    #print(Dz[i])
    p = (Dtwist-Twist[i])**2
    #print(p)
    i = i+1
    rmst.append(p)
    
MSEz = (sum(rmst)/len(rmst))**(1/2)
    
print("The MSE for the dZ is equal to",MSEz)
    
   
        
''''
plt.plot( hingeline1['x'] , hingeline1['dY'])
plt.plot(xver,Dy)
plt.legend(['Validation Model','Numerical Model'])
plt.title('Deflection in the y direction',fontsize=20)
plt.xlabel('x (m)',fontsize =15) 
plt.ylabel('y (m)',fontsize=15) 
plt.show()

'''
plt.plot( trailEdge1['x'] , trailEdge1['angle'])
plt.plot(xver,twist)
plt.title('Twist',fontsize=20)
plt.xlabel('x (m)',fontsize=15)
plt.ylabel('Twist angle (rad)',fontsize=15)


'''
# Plot showing the stresses / deflection with color

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
img = ax.scatter(case2['x'], case2['y']   , case2['z'], c= case1['dMag'], cmap='bone')
fig.colorbar(img)
plt.show()



#3d projection

threedee = plt.figure().gca(projection='3d')

threedee.scatter(coor.x, coor.y, coor.z, color = 'red')
threedee.set_xlabel('x')
threedee.set_ylabel('y')
threedee.set_zlabel('z')
plt.show()
# I want to find the hinge line, as I need them for deflection

print(max(case1['avgMises']))

print("Maximum Von-Mises-----------------------------------")
c = case1.loc[case1['avgMises'].idxmax()]
print(c['x'])



print("Maximum dY---------------------------------------")
print(case1.loc[case1['dY'].idxmax()])
print(max(abs((case1['dY']))))

print("Maximum dZ----------------------------------------")
print(case1.loc[case1['dZ'].idxmin()])
print(max(abs(case1['dZ'])))

print("Maximum Twist------------------------------------------")
print(trailEdge1.loc[trailEdge1['angle'].idxmin()])
print("Minimum angle is", min(trailEdge1['angle']))

print(max(abs(trailEdge1['angle'])))

#print(case1Stress)
'''

