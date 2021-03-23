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

function bool=membreCasse(Theta)
    epauleMax = %pi/2
    epauleMin = -140/180*%pi
    coudeMax = 145/180*%pi
    coudeMin = 0
    poigneMax = 15/180*%pi
    poigneMin = -45/180*%pi
    bool=(Theta(1)>epauleMax | Theta(1)<epauleMin | Theta(2)>coudeMax | Theta(2)<coudeMin | Theta(3)>poigneMax | Theta(3)<poigneMin)
endfunction

function v=V(t)
    v = 0.2
endfunction

function v=Voppose(t)
    v = -V(t)
endfunction

function Theta=ModeleGeometriqueInverse(alpha, Px, Py, v, dt, L1, L2, L3)
    i = 0
    Theta = matrix([0, 0, 0], 1, -1)
    vectWP = CylToCart([L3, alpha])
    pasfini = %T
    while pasfini
        Px = Px + dt*v(i*dt)
        i = i+1
        theta = Modelegeo(alpha, Px, Py, L1, L2, L3)
        if theta~=%F then
            Theta(i, :) = theta
            pasfini = ~membreCasse(Theta(i, :))
        else
            pasfini = %F
        end
    end
endfunction

function res=optimum(Py, L1, L2, L3, alpha)
    Px0 = (L1+L2+L3-Py)/2
    dt = 0.02
    
    Theta = ModeleGeometriqueInverse(alpha, Px0, Py, Voppose, dt, L1, L2, L3)
    
    n = length(Theta(:, 1))
    if (Theta(1, :)==[0, 0, 0] && n==1) then
        res = 0
    else
        P0 = CylToCart([L1, Theta(n, 1)]) + CylToCart([L2, Theta(n, 1)+Theta(n, 2)]) + CylToCart([L3, Theta(n, 1)+Theta(n, 2)+Theta(n, 3)])
        Px = P0(1)
        Py = P0(2)
        
        Theta = ModeleGeometriqueInverse(alpha, Px, Py, V, dt, L1, L2, L3)
        
        Pn = CylToCart([L1, Theta(n, 1)]) + CylToCart([L2, Theta(n, 1)+Theta(n, 2)]) + CylToCart([L3, Theta(n, 1)+Theta(n, 2)+Theta(n, 3)])
        
        res = Pn(1)-P0(1)
    end
endfunction

L1 = 0.3
L2 = 0.25
L3 = 0.05
alpha = %pi/2
Py = 0:0.01:(L1+L2+L3)
distance = []
argmax = 1
for i=1:length(Py)
    distance(i) = optimum(Py(i), L1, L2, L3, alpha)
    if distance(i)>distance(argmax) then
        argmax = i
    end
end

plot(Py, distance)
disp(Py(argmax))
