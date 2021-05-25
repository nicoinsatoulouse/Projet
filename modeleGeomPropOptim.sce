exec(pwd()+"\PlotSave.sci")
exec(pwd()+"\CartCyl.sci")

function Theta=Modelegeo(alpha,Px,Py,L1,L2,L3)
    l=[L1,L2,L3]
    Wx=Px-l(3)*cos(alpha)
    Wy=Py-l(3)*sin(alpha)
    y=(Wx^2+Wy^2-l(1)^2-l(2)^2)/(2*l(1)*l(2))
    if y<=1 & y>=-1 then
        T2 = atan(sqrt(1-y^2), y)
        T1 = atan(-l(2)*sin(T2)*Wx+(l(1)+l(2)*cos(T2))*Wy, l(2)*sin(T2)*Wy+(l(1)+l(2)*cos(T2))*Wx)
        T3=alpha-T1-T2
        Theta=[T1,T2,T3]
    else
        Theta=%F
    end
endfunction

function bool=membreCasse(Theta, alpha)
    epauleMax = %pi/2 + (%pi/2-alpha)
    epauleMin = -140/180*%pi + (%pi/2-alpha)
    coudeMax = 145/180*%pi
    coudeMin = 0
    poigneMax = 15/180*%pi
    poigneMin = -45/180*%pi
    bool=(Theta(1)>epauleMax | Theta(1)<epauleMin | Theta(2)>coudeMax | Theta(2)<coudeMin | Theta(3)>poigneMax | Theta(3)<poigneMin)
endfunction

function v=V(t)
    v = 20
endfunction

function v=Voppose(t)
    v = -V(t)
endfunction

function Pf=rotation(Pi, a)
    A = [cos(a), -sin(a); sin(a), cos(a)]
    Pf = (A*Pi')'
endfunction

function Theta=ModeleGeometriqueInverseRotationBase(alpha, Px, Py, v, dt, L1, L2, L3)
    P = rotation([Px, Py], %pi/2-alpha)
    Px = P(1)
    Py = P(2)
    i = 0
    Theta = matrix([0, 0, 0], 1, -1)
    vectWP = CylToCart([L3, alpha])
    pasfini = %T
    while pasfini
        Px = Px + dt*v(i*dt)
        i = i+1
        theta = Modelegeo(%pi/2, Px, Py, L1, L2, L3)
        if theta~=%F then
            Theta(i, :) = [theta(1)+alpha-%pi/2, theta(2), theta(3)]
            pasfini = ~membreCasse(theta, alpha)
        else
            pasfini = %F
        end
    end
endfunction

function a=angle(a)
    a = pmodulo(a + %pi, 2*%pi) - %pi
endfunction

function Thetakm=prop(Pxp,Pyp,L1,L2,L3,a,b,Gamma,Beta,Thetam,dt)
    C1 = cos(Thetam(1))
    C12 = cos(Thetam(1) + Thetam(2))
    C123m = cos(Thetam(1) + Thetam(2) + Thetam(3))
    S1 = sin(Thetam(1))
    S12 = sin(Thetam(1) + Thetam(2))
    S123m = sin(Thetam(1) + Thetam(2) + Thetam(3))
    J = [-L1*S1 - L2*S12 - L3*S123m, -L2*S12 - L3*S123m, -L3*S123m ;  L1*C1 + L2*C12 + L3*C123m, L2*C12 + L3*C123m, L3*C123m  ]
    Thetamp= pinv(J)*[Pxp ,Pyp]'
    Thetakm=[Thetam(1)+Thetamp(1), Thetam(2)+Thetamp(2), Thetam(3)+Thetamp(3), 0]
    Thetakm(4) = Beta - Thetakm(1) - Thetakm(2) - Thetakm(3)
    //theta4mSymetrie = 2*Beta - 2*sum(Thetakm(1:3)) - Thetakm(4)
    //Thetakm(4) = theta4mSymetrie
endfunction

function bool=membreCasseProp(Thetam, Gamma)
    epauleMax = %pi/2
    epauleMin = -140/180*%pi
    coudeMax = 145/180*%pi
    coudeMin = 0
    poigneMaxm = 15/180*%pi + %pi/2 - Gamma
    poigneMinm = -45/180*%pi + %pi/2 - Gamma
    bool=(Thetam(1)>epauleMax | Thetam(1)<epauleMin | Thetam(2)>coudeMax | Thetam(2)<coudeMin | Thetam(3)>poigneMaxm | Thetam(3)<poigneMinm)
endfunction

function Thetam=Modele(Beta, Thetam0,dt,L1,L2,L3,a,b,Gamma, A)
    i = 1
    Thetam = matrix(Thetam0, 1, -1)
    pasfini = %T
    while i<1000 && pasfini
        Pxp = A*cos(Beta)*i*dt
        Pyp = A*sin(Beta)*i*dt
        i = i+1
        thetam = prop(Pxp,Pyp,L1,L2,L3,a,b,Gamma,Beta,Thetam(i-1,:),dt)
        alpha = angle(sum(thetam(1:3)) + Gamma)
        /*if  alpha - %pi/2 - Beta <= 0 then //thetam(3) - %pi/2 + Gamma >= 0 then
            pasfini = %F
            Thetam(i, :) = thetam
        else*/
            Thetam(i, :) = thetam
            pasfini = ~membreCasseProp(Thetam(i, :), Gamma)
        //end
    end
    disp(pasfini)
endfunction

function Thetam=Prop(lH, lA, L1, L2, dt, Beta, r)
    /*L3=lH
    Lj = L1+L2
    Gamma = atan(lH, lA)
    
    alpha = Beta + %pi/2
    P = [(L1+L2+L3-r)/2, r]
    P = rotation(P, alpha-%pi/2)
    Px = P(1)
    Py = P(2)

    ThetaMarcheArr = ModeleGeometriqueInverseRotationBase(alpha, Px, Py, Voppose, dt, L1, L2, L3)
    n = max(length(ThetaMarcheArr(:, 1)), 2)

    L3 = sqrt(lH^2 + lA^2)

    Theta0 = ThetaMarcheArr(n - 1, :)
    Theta0m = [Theta0(1), Theta0(2), Theta0(3) + %pi/2 - Gamma, angle(%pi + Gamma)]
    Theta0m(3) = angle(Beta + %pi - Gamma - Theta0m(1) - Theta0m(2))
    
    theta3i = angle(Theta0m(3) - %pi/2 + Gamma)
    theta4i = angle(Theta0m(4) - Gamma)
    theta3f = -0.5212609//-0.7909668//0
    theta4f = -%pi/2
    a = 0.3*(theta4f - theta4i)/(theta3f - theta3i)
    b = theta4i - a*theta3i
    P0=CylToCart([L1,Theta0m(1)])+CylToCart([L2,Theta0m(1)+ Theta0m(2)])+ CylToCart([L3,Theta0m(1)+ Theta0m(2)+Theta0m(3)])
    
    Thetam = Modele(Beta, Theta0m, dt, L1, L2, L3, a, b, Gamma, 1)*/
    
    
    
L3=lH
Lj = L1+L2
Gamma = atan(lH, lA)
    
    alpha = Beta + %pi/2
P = [(L1+L2+L3-r)/2, r]
P = rotation(P, alpha-%pi/2)
Px = P(1)
Py = P(2)

ThetaMarcheArr = ModeleGeometriqueInverseRotationBase(alpha, Px, Py, Voppose, dt, L1, L2, L3)
n = max(length(ThetaMarcheArr(:, 1)), 2)

    
    L3 = sqrt(lH^2 + lA^2)

Theta0 = ThetaMarcheArr(n - 1, :)
Theta0m = [Theta0(1), Theta0(2), Theta0(3) + %pi/2 - Gamma, angle(%pi + Gamma)]
Theta0m(3) = angle(Beta + %pi - Gamma - Theta0m(1) - Theta0m(2))

theta3i = angle(Theta0m(3) - %pi/2 + Gamma)
theta4i = angle(Theta0m(4) - Gamma)
theta3f = -0.5212609//-0.7909668//0
theta4f = -%pi/2
disp(theta3i, theta4i, theta3f, theta4f)
a = 0.3(theta4f - theta4i)/(theta3f - theta3i)
b = theta4i - a*theta3i
P0=CylToCart([L1,Theta0m(1)])+CylToCart([L2,Theta0m(1)+ Theta0m(2)])+ CylToCart([L3,Theta0m(1)+ Theta0m(2)+Theta0m(3)])

Thetam = Modele(Beta, Theta0m, dt, L1, L2, L3, a, b, Gamma, 1)
endfunction

lH = 6.7
lA = 20
L1 = 28.6
L2 = 27.9
dt = 0.005

n=100
m=100
distance=[]
rmax = []
Beta=[]
for i=1:n+1
    distance(i) = 0
    Beta(i)=-%pi/2 + (i-1)/n*%pi
    d = []
    for j=1:m+1
        r = (j-1)/m*(L1 + L2 + lH)
        Thetam = Prop(lH, lA, L1, L2, dt, Beta(i), r)
        
        l = length(Thetam(:, 1))
        p0 = CylToCart([L1,Thetam(1,1)])+CylToCart([L2,Thetam(1,1)+ Thetam(1,2)])+ CylToCart([L3,Thetam(1,1)+Thetam(1,2)+Thetam(1,3)])
        p1 = CylToCart([L1,Thetam(l,1)])+CylToCart([L2,Thetam(l,1)+ Thetam(l,2)])+ CylToCart([L3,Thetam(l,1)+ Thetam(l,2)+Thetam(l,3)])
        d(j)=norm(p0-p1)
        disp(d(j))
    end
    [distance(i), index] = max(d)
    rmax(i) = (index-1)/m*(L1 + L2 + lH)
end

plot(Beta,distance)
show_window(1)
plot(Beta, rmax)
