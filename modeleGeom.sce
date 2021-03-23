exec(pwd()+"\resolutionEDO.sci")
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

L1 = 0.3
L2 = 0.25
L3 = 0.05
alpha = %pi/2
Py = 0.53
Px0 = (L1+L2+L3-Py)/2
dt = 0.02

Theta = ModeleGeometriqueInverse(alpha, Px0, Py, Voppose, dt, L1, L2, L3)

n = length(Theta(:, 1))
P = CylToCart([L1, Theta(n, 1)]) + CylToCart([L2, Theta(n, 1)+Theta(n, 2)]) + CylToCart([L3, Theta(n, 1)+Theta(n, 2)+Theta(n, 3)])
Px = P(1)
Py = P(2)

Theta = ModeleGeometriqueInverse(alpha, Px, Py, V, dt, L1, L2, L3)

Beta = %pi/2 - alpha
direct = CylToCart([1, Beta])
Lj = L1+L2
n = length(Theta(:, 1))

// Paramètres d'affichage
Lparam = list(init_param('N', 3, 'r1', L1, 'theta1', Theta(1, 1), 'r2', L2, 'theta2', Theta(1, 1)+Theta(1, 2), 'r3', L3, 'theta3', Theta(1, 1)+Theta(1, 2)+Theta(1, 3)), 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', Beta)
for i=2:n
    P = P + dt*V(i*dt)*direct
    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(i, 1), 'r2', L2, 'theta2', Theta(i, 1)+Theta(i, 2), 'r3', L3, 'theta3', Theta(i, 1)+Theta(i, 2)+Theta(i, 3), 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', Beta)
end
for i=n:(n+100)
    P = P + dt*V(i*dt)*direct
    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(n, 1), 'r2', L2, 'theta2', Theta(n, 1)+Theta(n, 2), 'r3', L3, 'theta3', Theta(n, 1)+Theta(n, 2)+Theta(n, 3), 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', Beta)
end

//Plot(NTiges, Lparam, dt, L1+L2+L3, %t)
PlotAndSave(NTiges, Lparam, 'BrasHorizontal.gif', dt, 10, L1+L2, %t)
