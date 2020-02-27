import numpy as np
from NEW_Forces_Deflections import locationvalue

def deflections(x,x1,x2,x3,xa,h,theta,P,E,G,zsc,Iyy,Izz,J,xlst,Vlst,Mlst,d1lst,defllst,Tlst,thetalst,Ay,Az,By,Bz,Cy,Cz,Fy,Fz,C1,C2,C3,C4,C5):

    # Transformation of some input parameters

    # Determine the total aerodynamic load
    Vtot = Vlst[-1]

    # Decompose force P along the y and z axis
    Py = P*np.sin(theta*np.pi/180)
    Pz = P*np.cos(theta*np.pi/180)

    # Determine the distance to the actuators
    xf = x2-xa/2
    xp = x2+xa/2

    # Calculate constants used in the displacement functions
    EIzz = E*Izz
    EIyy = E*Iyy
    GJ = G*J


    # Moment around the y axis
    def My(x):

        M = 0

        # Apply the Macaulay step function
        if x-x1>=0:
            M += -Az*(x-x1)
        if x-xf>=0:
            M += -Fz*(x-xf)
        if x-x2>=0:
            M += -Bz*(x-x2)
        if x-x3>=0:
            M += -Cz*(x-x3)

        if x-xp>=0:
            M += Pz*(x-xp)

        return M

    # Moment around the z axis
    def Mz(x):
        M = (-locationvalue(xlst, x, Mlst))  # Mlst is negative, needs to be +
        if x-x1>=0:
            M += -Ay*(x-x1)
        if x-xf>=0:
            M += -Fy*(x-xf)
        if x-x2>=0:
            M += -By*(x-x2)
        if x-x3>=0:
            M += -Cy*(x-x3)
        if x-xp>=0:
            M += Py*(x-xp)
        return M

    # Torque
    def T(x):
        Torque = locationvalue(xlst,x,Tlst)
        if x-x1>=0:
            Torque += Ay*(-(zsc+h/2))
        if x-xf>=0:
            Torque += Fy*(-zsc) - Fz*(h/2)
        if x-x2>=0:
            Torque += By*(-(zsc+h/2))
        if x-x3>=0:
            Torque += Cy*(-(zsc+h/2))
        if x-xp>=0:
            Torque += -Py*(-zsc)+Pz*h/2
        return Torque

    # Slope in y direction
    def vs(x):
        slope = locationvalue(xlst,x,d1lst)/EIzz + C1
        if x-x1>=0:
            slope += Ay/(2*EIzz)*(x-x1)**2
        if x-xf>=0:
            slope += Fy/(2*EIzz)*(x-xf)**2
        if x-x2>=0:
            slope += By/(2*EIzz)*(x-x2)**2
        if x-x3>=0:
            slope += Cy/(2*EIzz)*(x-x3)**2
        if x-xp>=0:
            slope += -Py/(2*EIzz)*(x-xp)**2
        return slope

    # Displacement in y direction
    def v(x):
        defl = locationvalue(xlst,x,defllst)/EIzz + C1*x + C2
        if x-x1>=0:
            defl += Ay/(6*EIzz)*(x-x1)**3
        if x-xf>=0:
            defl += Fy/(6*EIzz)*(x-xf)**3
        if x-x2>=0:
            defl += By/(6*EIzz)*(x-x2)**3
        if x-x3>=0:
            defl += Cy/(6*EIzz)*(x-x3)**3
        if x-xp>=0:
            defl += -Py/(6*EIzz)*(x-xp)**3
        return defl

    # Slope in z direction
    def ws(x):
        slope = C3
        if x-x1>=0:
            slope += Az/(2*EIyy)*(x-x1)**2
        if x-xf>=0:
            slope += Fz/(2*EIyy)*(x-xf)**2
        if x-x2>=0:
            slope += Bz/(2*EIyy)*(x-x2)**2
        if x-x3>=0:
            slope += Cz/(2*EIyy)*(x-x3)**2
        if x-xp>=0:
            slope += -Pz/(2*EIyy)*(x-xp)**2
        return slope

    # Displacement in z direction
    def w(x):
        defl = C3*x + C4
        if x-x1>=0:
            defl += Az/(6*EIyy)*(x-x1)**3
        if x-xf>=0:
            defl += Fz/(6*EIyy)*(x-xf)**3
        if x-x2>=0:
            defl += Bz/(6*EIyy)*(x-x2)**3
        if x-x3>=0:
            defl += Cz/(6*EIyy)*(x-x3)**3
        if x-xp>=0:
            defl += -Py/(6*EIyy)*(x-xp)**3
        return defl

    # Angle
    def ftheta(x):
        twist = locationvalue(xlst,x,thetalst)/GJ + C5
        if x-x1>=0:
            twist += Ay/GJ*(-(zsc+h/2))*(x-x1)
        if x-xf>=0:
            twist += Fy/GJ*(-zsc)*(x-xf) - Fz/GJ*h/2*(x-xf)
        if x-x2>=0:
            twist += By/GJ*(-(zsc+h/2))*(x-x2)
        if x-x3>=0:
            twist += Cy/GJ*(-(zsc+h/2))*(x-x3)
        if x-xp>=0:
            twist += -Py/GJ*(-zsc)*(x-xp)+Pz/GJ*h/2*(x-xp)
        return twist

    # Sum of forces in z direction
    def Sz(x):
        Sumz = 0
        if x-x1>=0:
            Sumz += -Az
        if x-xf>=0:
            Sumz += -Fz
        if x-x2>=0:
            Sumz += -Bz
        if x-x3>=0:
            Sumz += -Cz
        if x-xp>=0:
            Sumz += Pz
        return Sumz

    # Sum of forces in y direction
    def Sy(x):
        Sumy = -Vtot
        if x-x1>=0:
            Sumy += -Ay
        if x-xf>=0:
            Sumy += -Fy
        if x-x2>=0:
            Sumy += -By
        if x-x3>=0:
            Sumy += -Cy
        if x-xp>=0:
            Sumy += Py
        return Sumy


    Myx = My(x)
    Mzx = Mz(x)
    Tx = T(x)
    Syx = Sy(x)
    Szx = Sz(x)
    vsx = vs(x)
    vx = v(x)
    wsx = ws(x)
    wx = w(x)
    thetax = ftheta(x)

    return Myx,Mzx,Tx,Syx,Szx,vsx,vx,wsx,wx,thetax