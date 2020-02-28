# Imports
import numpy as np
from Interpolation import patchinterpolate
from NEW_Forces_Deflections import output, locationvalue
from Reaction_Forces import reaction_forces
from Cross_sectional_properties import cross_sect_prop
from ShearCalcs import CalcTorsStif
import matplotlib.pyplot as plt
from Deflections import deflections
from AllStressCalcs import CalcShearCenter,AileronFullStress


# Input parameters
Ca    = 0.505       #Chord Length Aileron
la    = 1.611       #Span of the Aileron
x1    = 0.125       #x-location of hinge 1
x2    = 0.498       #x-location of hinge 2
x3    = 1.494       #x-location of hinge 3
xa    = 24.5E-2     #Distance between actuator 1 and 2
h     = 16.1E-2     #Aileron Height
tsk   = 1.1E-3      #Skin thickness
tsp   = 2.4E-3      #Spar Thickness
tst   = 1.2E-3      #Thickness of stiffener
hst   = 1.3E-2      #Height of stifffener
wst   = 1.7E-2      #Width of stiffeners
nst   = 11          #Number of stiffeners (equalt spaced)
d1    = 0.389E-2    #Vertical displacement of hinge 1
d3    = 1.245E-2    #Vertical displacement of hinge 2
theta = 30          #Maximum upward deflection
P     = 49.2E3      #Load in actuator 2
E     = 73.1E9      # [Pa]
G     = 28E9        # [Pa]


# Calculated parameters
J     = CalcTorsStif(Ca,h,tsk,tsp)
Dx, zbar,ybar,Izz,Iyy, stiff = cross_sect_prop(h,tsk,tsp,tst,hst,wst,Ca,nst)
#zsc = -0.085       # [m]
zsc = CalcShearCenter(h,Ca,tsk,tsp,Iyy,Izz,stiff,Dx,zbar)


# Reading the aerodynamic load
data = np.loadtxt('AERO.dat',delimiter = ',')


# Transforming the coordinates
xl = []
zl = []
Nz = len(data)
Nx = len(data[0])

for i in range(1,Nz+1):
    thetazi  = (i-1)/Nz*np.pi
    thetazi1 = (i)/Nz*np.pi
    z  =  -(1/2)*((Ca/2)*(1-np.cos(thetazi))+(Ca/2)*(1-np.cos(thetazi1)))
    zl.append(z)

for k in range(1,Nx+1):
    thetaxk  = (k-1)/Nx*np.pi
    thetaxk1 = (k)/Nx*np.pi
    x  = (1/2)*((la/2)*(1-np.cos(thetaxk))+(la/2)*(1-np.cos(thetaxk1)))
    xl.append(x)


# Interpolate the aerodynamic load
xlst, zlst, qlst = patchinterpolate(600,600,xl,zl,data)


# Integrate the aerodynamic load
Vlst, Mlst, d1lst, defllst, taulst, Tlst, thetalst = output(xlst, zlst, qlst, zsc)

# Calculate the reaction forces
Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5 = reaction_forces(la,x1,x2,x3,xa,h,d1,d3,theta,P,E,G,zsc,Iyy,
                                                         Izz, J,xlst,Vlst,Mlst,defllst,Tlst,thetalst)

# ------------------------------------------------------------------------------------------------

x = np.linspace(0,la,num = 10)


Mz_vm = np.array([  -357.74357985,  -1594.88443362,  -6651.88885912, -14347.45149789,
       -12912.13233209, -10008.39323038,  -6945.42913784,  -3687.82456212,
         -350.7487045 ,    -73.38058827])
My_vm = np.array([ 1650.44463845, 14665.39197011, 60633.17165021, 80288.98751996,
       61089.45919142, 46494.93546666, 31887.77377816, 16490.97929878,
        1447.92514485,  -414.23138736])
T_vm = np.array([  -99.35419573,   114.26077483,  -477.38919173, -1620.68142007,
        -238.32321721,  -284.39024956,  -210.54230859,  -140.5999087 ,
         -75.23974716,    16.02929465])

Sy_vm = np.array([ 77292.90447293, -31297.19426669, -43902.48206311, -14208.48844372,
        15405.13281629,  16912.12658215,  16944.71270833,  17158.24114619,
        13341.17566935, -13829.13786589])
Sz_vm = np.array([-392038.51188887,  272732.11952853,  260236.00177141,
        -89633.7566658 ,  -68567.88699836,  -91121.42454922,
        -81401.61994014,  -72882.30002097,  -53687.08850463,
        -87871.85096728])

v_vm = np.array([ 4.71348220e-03,  2.71155267e-03,  8.86760041e-04, -2.33103361e-04,
        -4.26927759e-05,  1.40946177e-03,  3.84118281e-03,  6.95039405e-03,
         1.04198310e-02,  1.39618003e-02])
w_vm = np.array([-2.81096280e-03, -1.53506669e-03, -4.34174883e-04,  6.62846336e-05,
        -2.10286374e-04, -1.11228482e-03, -2.48715156e-03, -4.18239589e-03,
        -6.04473649e-03, -7.94045847e-03])
theta_vm = np.array([-0.00255686, -0.00251688, -0.00242618, -0.00366564, -0.00433512,
        -0.00456606, -0.00475794, -0.00491263, -0.00502382, -0.00503579])

Mz_nm = []
My_nm = []
T_nm = []

Sy_nm = []
Sz_nm = []

v_nm = []
w_nm = []
theta_nm = []


for xi in x:
    Myx, Mzx, Tx, Syx, Szx, vsx, vx, wsx, wx, thetax = deflections(xi, x1, x2, x3, xa, h, theta, P, E, G, zsc, Iyy, Izz,
                                                                   J, xlst, Vlst, Mlst, d1lst, defllst, Tlst, thetalst,
                                                                   Ay, Az, By, Bz, Cy, Cz, Fy, Fz, C1, C2, C3, C4, C5)
    Mz_nm.append(Mzx)
    My_nm.append(Myx)
    T_nm.append(Tx)

    Sy_nm.append(Syx)
    Sz_nm.append(Szx)

    v_nm.append(vx)
    w_nm.append(wx)
    theta_nm.append(thetax)

Mz_nm = np.array(Mz_nm)
My_nm = np.array(My_nm)
T_nm = np.array(T_nm)

Sy_nm = np.array(Sy_nm)
Sz_nm = np.array(Sz_nm)

v_nm = np.array(v_nm)
w_nm = np.array(w_nm)
theta_nm = np.array(theta_nm)

RMSE_Mz = np.sqrt(1/len(x)*np.sum((Mz_nm-Mz_vm)**2))
RMSE_My = np.sqrt(1/len(x)*np.sum((My_nm-My_vm)**2))
RMSE_T = np.sqrt(1/len(x)*np.sum((T_nm-T_vm)**2))

RMSE_Sy = np.sqrt(1/len(x)*np.sum((Sy_nm-Sy_vm)**2))
RMSE_Sz = np.sqrt(1/len(x)*np.sum((Sz_nm-Sz_vm)**2))

RMSE_v = np.sqrt(1/len(x)*np.sum((v_nm-v_vm)**2))
RMSE_w = np.sqrt(1/len(x)*np.sum((w_nm-w_vm)**2))
RMSE_theta = np.sqrt(1/len(x)*np.sum((theta_nm-theta_vm)**2))

print('RMSE Mz =',RMSE_Mz)
print('RMSE My =',RMSE_My)
print('RMSE T =',RMSE_T)

print('RMSE Sy =',RMSE_Sy)
print('RMSE Sz =',RMSE_Sz)

print('RMSE v =',RMSE_v)
print('RMSE w =',RMSE_w)
print('RMSE theta =',RMSE_theta)


plt.plot(x,Mz_nm,label='Numerical Model')
plt.plot(x,Mz_vm,label='Verification Model')
plt.title('Moment around the z axis')
plt.xlabel('x [m]')
plt.ylabel('Mz(x) [Nm]')
plt.legend()
plt.grid()
plt.show()

plt.plot(x,My_nm,label='Numerical Model')
plt.plot(x,My_vm,label='Verification Model')
plt.title('Moment around the y axis')
plt.xlabel('x [m]')
plt.ylabel('My(x) [Nm]')
plt.legend()
plt.grid()
plt.show()

plt.plot(x,T_nm,label='Numerical Model')
plt.plot(x,T_vm,label='Verification Model')
plt.title('Torque')
plt.xlabel('x [m]')
plt.ylabel('T(x) [Nm]')
plt.legend()
plt.grid()
plt.show()


plt.plot(x,Sy_nm,label='Numerical Model')
plt.plot(x,Sy_vm,label='Verification Model')
plt.title('Shear force in y direction')
plt.xlabel('x [m]')
plt.ylabel('Sy(x) [N]')
plt.legend()
plt.grid()
plt.show()

plt.plot(x,Sz_nm,label='Numerical Model')
plt.plot(x,Sz_vm,label='Verification Model')
plt.title('Shear force in z direction')
plt.xlabel('x [m]')
plt.ylabel('Sz(x) [N]')
plt.legend()
plt.grid()
plt.show()


plt.plot(x,v_nm,label='Numerical Model')
plt.plot(x,v_vm,label='Verification Model')
plt.title('Deflection in y direction')
plt.xlabel('x [m]')
plt.ylabel('v(x) [m]')
plt.legend()
plt.grid()
plt.show()

plt.plot(x,w_nm,label='Numerical Model')
plt.plot(x,w_vm,label='Verification Model')
plt.title('Deflection in z direction')
plt.xlabel('x [m]')
plt.ylabel('w(x) [m]')
plt.legend()
plt.grid()
plt.show()

plt.plot(x,theta_nm,label='Numerical Model')
plt.plot(x,theta_vm,label='Verification Model')
plt.title('Twist')
plt.xlabel('x [m]')
plt.ylabel('theta(x) [rad]')
plt.legend()
plt.grid()
plt.show()