exec(pwd()+"\PlotSave.sci")
exec(pwd()+"\CartCyl.sci")

lH = 6.7
lA = 20
L1 = 28.6
L2 = 27.9
L3 = sqrt(lH^2 + lA^2)
Lj = L1+L2
dt=0.01
Gamma = atan(lH, lA)
Beta=%pi/4
x1 = -%pi/4
y1 = -%pi+atan(lH, lA)
x2 = 0
y2 = -%pi/2
a = (y2-y1)/(x2-x1)
b = y1 - a*x1
Theta0 = [0, 145/180*%pi/2, 15/180*%pi/2,%pi]
P0=CylToCart([L1,Theta0(1)])+CylToCart([L2,Theta0(1)+ Theta0(2)])+ CylToCart([L3,Theta0(1)+ Theta0(2)+Theta0(3)+%pi/2-Gamma])

function Thetak=Modelegeo(Pxp,Pyp,L1,L2,L3,a,b,Beta,Theta,dt)
    C1 = cos(Theta(1))
    C12 = cos(Theta(1) + Theta(2))
    C123 = cos(Theta(1) + Theta(2) + Theta(3))
    S1 = sin(Theta(1))
    S12 = sin(Theta(1) + Theta(2))
    S123 = sin(Theta(1) + Theta(2) + Theta(3))
    J=[ - L1*S1 - L2*S12 - L3*S123, - L2*S12 - L3*S123, - L3*S123 ;  L1*C1 + L2*C12 + L3*C123, L2*C12 + L3*C123, L3*S123  ]
    Thetap= pinv(J)*[Pxp ,Pyp]'
    Thetak=[Theta(1)+Thetap(1), Theta(2)+Thetap(2), Theta(3)+Thetap(3), 0]
    Thetak(4) = a*Thetak(3) + b
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

function Theta=ModeleGeometriqueInverse(Beta, Theta0,dt,L1,L2,L3,a,b,Gamma, A)
    i = 1
    Theta = matrix(Theta0, 1, -1)
    pasfini = %T
    while pasfini
        Pxp = A*cos(Beta)*i*dt
        Pyp = A*sin(Beta)*i*dt
        i = i+1
        theta = Modelegeo(Pxp,Pyp,L1,L2,L3,a,b,Beta,Theta(i-1,:),dt)
        alpha=sum(theta)+ %pi/2 -Gamma
        if  alpha < Beta then 
            pasfini = %F
        else 
            Theta(i, :) = theta
            pasfini = ~membreCasse(Theta(i, :))
        end
    end
endfunction

Theta=ModeleGeometriqueInverse( Beta, Theta0,dt,L1,L2,L3,a,b,Gamma, -1)
Theta=ModeleGeometriqueInverse( Beta, Theta(length(Theta(:, 1)),:),dt,L1,L2,L3,a,b,Gamma, 1)


// Paramètres d'affichage

n = length(Theta(:, 1))


Lparam = list(init_param('N', 3, 'r1', L1, 'theta1', Theta(1, 1), 'r2', L2, 'theta2', Theta(1, 1)+Theta(1, 2), 'r3', L3, 'theta3', Theta(1, 1)+Theta(1, 2)+Theta(1, 3) + %pi/2 - Gamma, 'Javelot', %T, 'xJavelot', P0(1), 'yJavelot', P0(2), 'LJavelot', Lj, 'thetaJavelot', sum(Theta(1))))

for i=2:n
    P = P0+[cos(Beta)*(i*dt)^2, sin(Beta)*(i*dt)^2]
    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(i, 1), 'r2', L2, 'theta2', Theta(i, 1)+Theta(i, 2), 'r3', L3, 'theta3', Theta(i, 1)+Theta(i, 2)+Theta(i, 3) + %pi/2 - Gamma, 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', sum(Theta(i)))
end

for i=n:(n+100)
    P = P0+ [cos(Beta)*(i*dt)^2, sin(Beta)*(i*dt)^2]
    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(n, 1), 'r2', L2, 'theta2', Theta(n, 1)+Theta(n, 2), 'r3', L3, 'theta3', Theta(n, 1)+Theta(n, 2)+Theta(n, 3) + %pi/2 - Gamma, 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', Beta)
end

//Plot(NTiges, Lparam, dt, L1+L2+L3, %t)
PlotAndSave(NTiges, Lparam, 'testprop.gif', dt, 10, L1+L3+L2, %t)
