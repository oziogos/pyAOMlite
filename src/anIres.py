import math
res=0
def anIres(X1,Y1,Z1,a1,lx1,ly1,lz1,X2,Y2,Z2,a2,lx2,ly2,lz2):
    global res
    try:
        if 1 in (lx1,ly1,lz1) and 1 in (lx2,ly2,lz2):
            if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
                res=(4*1.4142135623730951*(a1*a2)**1.25*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
                res=(-8*1.4142135623730951*(a1*a2)**2.25*(Y1-Y2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
                res=(8*1.4142135623730951*a1**2.25*a2**0.75*a2*(Y1-Y2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
                res=(-8*1.4142135623730951*(a1*a2)**2.25*(X1-X2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
                res=(8*1.4142135623730951*a1**2.25*a2**0.75*a2*(X1-X2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
                res=(-16*1.4142135623730951*a1**3.25*a2**1.75*a2*(X1-X2)*(Y1-Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(-8*1.4142135623730951*(a1*a2)**2.25*(Y1-Y2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(4*1.4142135623730951*(a1*a2)**1.25*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(8*1.4142135623730951*a1**2.25*a2**0.75*a2*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(-8*1.4142135623730951*(a1*a2)**2.25*(X1-X2)*(Y1-Y2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(16*1.4142135623730951*a1**3.25*a2**1.75*a2*(X1-X2)*(-Y1+Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(8*1.4142135623730951*a1**2.25*a2**0.75*a2*(X1-X2)*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**2.25*(-Y1+Y2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
                res=(-8*1.4142135623730951*a1**0.75*a1*a2**2.25*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2)))*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
                res=(16*1.4142135623730951*a1**1.75*a1*a2**3.25*(X1-X2)*(Y1-Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
                res=(16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(X1-X2)*(-Y1+Y2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
                res=(-16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(X1-X2)*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(-8*1.4142135623730951*(a1*a2)**2.25*(X1-X2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(8*1.4142135623730951*(a1*a2)**2.25*(-X1+X2)*(Y1-Y2))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(-16*1.4142135623730951*a1**3.25*a2**1.75*a2*(X1-X2)*(Y1-Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(4*1.4142135623730951*(a1*a2)**1.25*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2))))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(8*1.4142135623730951*a1**2.25*a2**0.75*a2*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(-8*1.4142135623730951*a1**2.25*a2**0.75*a2*(-a2+a1*(-1+2*a2*(X1-X2)*(X1-X2)))*(Y1-Y2))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==1:
                res=(-8*1.4142135623730951*a1**0.75*a1*a2**2.25*(X1-X2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==0:
                res=(16*1.4142135623730951*a1**1.75*a1*a2**3.25*(X1-X2)*(Y1-Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==1 and lz2==1:
                res=(-16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(X1-X2)*(Y1-Y2)*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==0:
                res=(-8*1.4142135623730951*a1**0.75*a1*a2**2.25*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==0 and lz2==1:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(a2+a1*(1-2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==1 and ly2==1 and lz2==0:
                res=(-16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(Y1-Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(16*1.4142135623730951*a1**1.75*a1*a2**3.25*(X1-X2)*(Y1-Y2)*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(-8*1.4142135623730951*a1**0.75*a1*a2**2.25*(X1-X2)*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(-16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(X1-X2)*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**2.25*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(-Y1+Y2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(16*1.4142135623730951*a1**1.75*a1*a2**1.75*a2*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(-Y1+Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a2+a1*(1-2*a2*(X1-X2)*(X1-X2)))*(a2+a1*(1-2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
        elif 1 in (lx1,ly1,lz1) and 2 in (lx2,ly2,lz2):
            if lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
                res=(-8*0.816496580927726*a1**1.25*a2**0.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
                res=(-8*0.816496580927726*a1**1.25*a2**1.75*a2*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
                res=(-8*0.816496580927726*a1**1.25*a2**1.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(8*0.816496580927726*a1**1.25*a2**1.75*a2*(-Y1+Y2)*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(-8*0.816496580927726*a1**1.25*a2**0.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(8*0.816496580927726*a1**1.25*a2**1.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2))*(-Y1+Y2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(Y1-Y2)*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**2.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2))*(Y1-Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(8*0.816496580927726*a1**1.25*a2**1.75*a2*(-X1+X2)*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(-8*0.816496580927726*a1**1.25*a2**1.75*a2*(X1-X2)*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2)))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(-8*0.816496580927726*a1**1.25*a2**0.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(X1-X2)*(X1-X2)))*(X1-X2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==2:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(X1-X2)*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==2 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**2.75*a2*(X1-X2)*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==2 and ly2==0 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(X1-X2)*(X1-X2)))*(X1-X2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(16*0.816496580927726*a1**0.75*a1*a2**2.75*a2*(X1-X2)*(Y1-Y2)*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(X1-X2)*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(16*0.816496580927726*a1**0.75*a1*a2**1.75*a2*(-(a1*a2)+a2*a2+2*a1*a1*(-1+a2*(X1-X2)*(X1-X2)))*(X1-X2)*(Y1-Y2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
        elif 2 in (lx1,ly1,lz1) and 1 in (lx2,ly2,lz2):
            if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==1:
                res=(8*0.816496580927726*a1**0.75*a1*a2**1.25*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==1 and lz2==0:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(Y1-Y2)*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==1 and lz2==1:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(Y1-Y2)*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==0 and lz2==0:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(X1-X2)*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==0 and lz2==1:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(X1-X2)*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Z1-Z2)*(Z1-Z2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==1 and ly2==1 and lz2==0:
                res=(16*0.816496580927726*a1**2.75*a1*a2**0.75*a2*(X1-X2)*(Y1-Y2)*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**5.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(8*0.816496580927726*a1**0.75*a1*a2**1.25*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(X1-X2)*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(16*0.816496580927726*a1**2.75*a1*a2**0.75*a2*(X1-X2)*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2)))*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(X1-X2)*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(Y1-Y2)*(Y1-Y2)))*(Y1-Y2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2)))*(Z1-Z2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(8*0.816496580927726*a1**1.75*a1*a2**1.25*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2)))*(Y1-Y2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(16*0.816496580927726*a1**2.75*a1*a2**0.75*a2*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2)))*(Y1-Y2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(8*0.816496580927726*a1**0.75*a1*a2**1.25*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(X1-X2)*(X1-X2)))*(X1-X2))/((a1+a2)**4.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(X1-X2)*(X1-X2)))*(X1-X2)*(Z1-Z2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(16*0.816496580927726*a1**1.75*a1*a2**0.75*a2*(a1*a1-2*a2*a2+a1*a2*(-1+2*a2*(X1-X2)*(X1-X2)))*(X1-X2)*(Y1-Y2))/((a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
        elif 2 in (lx1,ly1,lz1) and 2 in (lx2,ly2,lz2):
            if lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==2:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(-6*a1*a2*(-1+a2*(Z1-Z2)*(Z1-Z2))+a2*a2*(3+2*a2*(Z1-Z2)*(Z1-Z2))+a1*a1*(3-6*a2*(Z1-Z2)*(Z1-Z2)+4*a2*a2*(Z1-Z2)**4)+2*a1*a1*a1*(Z1-Z2)*(Z1-Z2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==2 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2))*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==2 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2))*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2)))*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(-6*a1*a2*(-1+a2*(Y1-Y2)*(Y1-Y2))+a2*a2*(3+2*a2*(Y1-Y2)*(Y1-Y2))+a1*a1*(3-6*a2*(Y1-Y2)*(Y1-Y2)+4*a2*a2*(Y1-Y2)**4)+2*a1*a1*a1*(Y1-Y2)*(Y1-Y2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2))*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2))))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2)))*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2)))*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**0.75*a2*(-6*a1*a2*(-1+a2*(X1-X2)*(X1-X2))+a2*a2*(3+2*a2*(X1-X2)*(X1-X2))+a1*a1*(3-6*a2*(X1-X2)*(X1-X2)+4*a2*a2*(X1-X2)**4)+2*a1*a1*a1*(X1-X2)*(X1-X2)))/(3.*(a1+a2)**5.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
        else:
            if lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(2*1.4142135623730951*a1**0.75*a2**0.75)/((a1+a2)**1.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==1:
                res=(4*1.4142135623730951*a1**1.75*a2**1.25*(Z1-Z2))/((a1+a2)**2.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==2:
                res=(4*0.816496580927726*a1**0.75*a2**0.75*a2*(a1+a2+2*a1*a1*(Z1-Z2)*(Z1-Z2)))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==0:
                res=(4*1.4142135623730951*a1**1.75*a2**1.25*(Y1-Y2))/((a1+a2)**2.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==1 and lz2==1:
                res=(8*1.4142135623730951*a1**2.75*a2**0.75*a2*(Y1-Y2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==0 and ly2==2 and lz2==0:
                res=(4*0.816496580927726*a1**0.75*a2**0.75*a2*(a1+a2+2*a1*a1*(Y1-Y2)*(Y1-Y2)))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==0:
                res=(4*1.4142135623730951*a1**1.75*a2**1.25*(X1-X2))/((a1+a2)**2.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==0 and lz2==1:
                res=(8*1.4142135623730951*a1**2.75*a2**0.75*a2*(X1-X2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==1 and ly2==1 and lz2==0:
                res=(8*1.4142135623730951*a1**2.75*a2**0.75*a2*(X1-X2)*(Y1-Y2))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==0 and lx2==2 and ly2==0 and lz2==0:
                res=(4*0.816496580927726*a1**0.75*a2**0.75*a2*(a1+a2+2*a1*a1*(X1-X2)*(X1-X2)))/((a1+a2)**3.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
                res=(4*1.4142135623730951*a1**1.25*a2**1.75*(-Z1+Z2))/((a1+a2)**2.5*math.exp((a1*a2*((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2)))/(a1+a2)))
            elif lx1==0 and ly1==0 and lz1==2 and lx2==0 and ly2==0 and lz2==0:
                res=(4*0.816496580927726*a1**0.75*a1*a2**0.75*(a1+a2*(1+2*a2*(Z1-Z2)*(Z1-Z2))))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(4*1.4142135623730951*a1**1.25*a2**1.75*(-Y1+Y2))/((a1+a2)**2.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==1 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**2.75*(Y1-Y2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==0 and ly1==2 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(4*0.816496580927726*a1**0.75*a1*a2**0.75*(a1+a2*(1+2*a2*(Y1-Y2)*(Y1-Y2))))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(-4*1.4142135623730951*a1**1.25*a2**1.75*(X1-X2))/((a1+a2)**2.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==0 and lz1==1 and lx2==0 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**2.75*(X1-X2)*(Z1-Z2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==1 and ly1==1 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(8*1.4142135623730951*a1**0.75*a1*a2**2.75*(X1-X2)*(Y1-Y2))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
            elif lx1==2 and ly1==0 and lz1==0 and lx2==0 and ly2==0 and lz2==0:
                res=(4*0.816496580927726*a1**0.75*a1*a2**0.75*(a1+a2*(1+2*a2*(X1-X2)*(X1-X2))))/((a1+a2)**3.5*math.exp((a1*a2*(X1*X1-2*X1*X2+X2*X2+Y1*Y1-2*Y1*Y2+Y2*Y2+Z1*Z1-2*Z1*Z2+Z2*Z2))/(a1+a2)))
    except OverflowError as err:
        pass
    return res
