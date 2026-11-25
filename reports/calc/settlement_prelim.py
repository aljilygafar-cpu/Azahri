import math
# --- Boussinesq: vertical stress under corner of rectangle (Newmark) ---
def Icorner(B,L,z):
    if B<=0 or L<=0: return 0.0
    m,n=B/z,L/z
    V=m*m+n*n+1; V1=(m*n)**2
    a=2*m*n*math.sqrt(V)/(V+V1)*(V+1)/V
    den=V-V1
    t=math.atan(2*m*n*math.sqrt(V)/den) if den!=0 else math.pi/2
    if den<0: t+=math.pi
    return (a+t)/(4*math.pi)
def Ipoint(x,y,x1,x2,y1,y2,z):
    # influence at point (x,y) of uniform load on rectangle; superposition with signs
    def part(dx,dy):  # rectangle from point to (dx,dy) signed
        s=(1 if dx>=0 else -1)*(1 if dy>=0 else -1)
        return s*Icorner(abs(dx),abs(dy),z)
    return part(x2-x,y2-y)-part(x1-x,y2-y)-part(x2-x,y1-y)+part(x1-x,y1-y)
# footings: name -> list of rectangles (x1,x2,y1,y2), net pressure kPa
F={
 'F3 front-L':([(1.80,4.60,-1.41,1.40)],150),
 'F3 front-2':([(6.97,9.77,-1.41,1.40)],150),
 'F5':([(10.25,12.75,2.31,6.71)],150),
 'F3 east':([(14.98,17.78,3.80,6.60)],150),
 'F6':([(1.90,4.51,3.90,6.49),(0.70,3.30,6.49,8.62)],120),
 'F7':([(6.78,9.98,3.60,6.80),(5.10,8.30,6.80,8.92)],120),
 'F1 upper':([(9.65,13.35,8.35,12.05)],155),
 'F2 upper':([(14.85,18.15,8.55,11.85)],150),
 'F8':([(0.70,8.00,9.90,17.10)],130),
 'F1 lower':([(9.77,13.47,13.94,17.64)],155),
 'F2 lower':([(15.35,18.65,14.14,17.44)],150),
 'F4 L':([(0.30,2.70,19.67,22.07)],150),
 'F3 rear':([(5.18,7.97,19.59,22.39)],150),
 'F3 rear-2':([(10.23,13.02,19.59,22.39)],150),
 'F4 R':([(15.79,18.20,19.79,22.19)],150),
}
pts={ # evaluation points (column positions / centers)
 'F8 col SW(2.13,15.79)':(2.13,15.79),'F8 col SE(6.57,15.79)':(6.57,15.79),'F8 col NE(6.57,11.20)':(6.57,11.20),'F8 col NW(2.13,11.20)':(2.13,11.20),'F8 center':(4.35,13.50),
 'F4 L col':(1.50,20.87),'F3 rear col':(6.57,20.99),'F1 lower col':(11.62,15.79),'F1 upper col':(11.50,10.19),
 'F6 lower col':(2.00,7.32),'F7 lower col':(6.70,7.32),'F7 upper col':(8.38,5.20),'F3 front-2 col':(8.38,0.0),
 'F2 upper col':(16.50,10.19),'F2 lower col':(17.00,15.79),'F3 rear-2 col':(11.62,20.99),'F5 C1':(11.50,5.32),'F3 east col':(16.37,5.20),'F4 R col':(17.0,20.99),
}
def settlement(x,y,Eclay,Esand,Efill=50e3,zmax=15.0,dz=0.05):
    s=0.0; z=dz/2
    while z<zmax:
        ds=0.0
        for name,(rects,q) in F.items():
            for (x1,x2,y1,y2) in rects:
                ds+=q*Ipoint(x,y,x1,x2,y1,y2,z)
        E=Efill if z<0.6 else (Eclay if z<2.5 else Esand)
        s+=ds*dz/E
        z+=dz
    return s*1000  # mm
cases=[('soft-ish clay 7MPa, sand 20MPa',7e3,20e3),('firm clay 12MPa, sand 30MPa',12e3,30e3),('stiff clay 20MPa, sand 40MPa',20e3,40e3)]
res={}
for lab,Ec,Es in cases:
    res[lab]={k:settlement(*v,Ec,Es) for k,v in pts.items()}
    print('==',lab)
    for k,v in res[lab].items(): print(f'  {k:24s} {v:6.1f} mm')
pairs=[('F8 col SW(2.13,15.79)','F4 L col'),('F8 col SE(6.57,15.79)','F3 rear col'),('F8 col SE(6.57,15.79)','F1 lower col'),('F8 col NE(6.57,11.20)','F1 upper col'),('F8 col NW(2.13,11.20)','F6 lower col'),('F8 col NE(6.57,11.20)','F7 lower col'),('F7 upper col','F3 front-2 col'),('F1 upper col','F2 upper col'),('F1 lower col','F3 rear-2 col'),('F1 upper col','F5 C1')]
print('== differential (mm) and angular distortion')
for a,b in pairs:
    xa,ya=pts[a]; xb,yb=pts[b]; L=math.hypot(xa-xb,ya-yb)
    out=[]
    for lab,_,_ in cases:
        d=abs(res[lab][a]-res[lab][b]); out.append(f'{d:5.1f}mm 1/{(L*1000/d) if d>0 else 0:5.0f}')
    print(f'{a[:14]:14s}-{b[:14]:14s} L={L:4.2f}m  ', ' | '.join(out))
