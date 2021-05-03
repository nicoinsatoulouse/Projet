exec(pwd()+"\PlotSave.sci")
exec(pwd()+"\CartCyl.sci")

function a=angle(a)
    a = pmodulo(a + %pi, 2*%pi) - %pi
endfunction

lH = 6.7
lA = 20
L1 = 28.6
L2 = 27.9
L3 = sqrt(lH^2 + lA^2)
Lj = L1+L2
dt=0.01
Gamma = atan(lH, lA)
Beta = 0//-%pi/8
Theta0m = [Beta, 5*%pi/8, 0, 0]
Theta0m(4) = angle(%pi + Gamma)
Theta0m(3) = angle(Beta + %pi - Gamma - Theta0m(1) - Theta0m(2))
theta3i = angle(Theta0m(3) - %pi/2 + Gamma)
theta4i = angle(Theta0m(4) - Gamma)
theta3f = -0.7909668//-0.7909668//0
theta4f = -%pi/2
disp(theta3i, theta4i, theta3f, theta4f)
a = (theta4f - theta4i)/(theta3f - theta3i)
b = theta4i - a*theta3i
P0=CylToCart([L1,Theta0m(1)])+CylToCart([L2,Theta0m(1)+ Theta0m(2)])+ CylToCart([L3,Theta0m(1)+ Theta0m(2)+Theta0m(3)])

function Thetakm=Modelegeo(Pxp,Pyp,L1,L2,L3,a,b,Gamma,Beta,Thetam,dt)
    C1 = cos(Thetam(1))
    C12 = cos(Thetam(1) + Thetam(2))
    C123m = cos(Thetam(1) + Thetam(2) + Thetam(3))
    S1 = sin(Thetam(1))
    S12 = sin(Thetam(1) + Thetam(2))
    S123m = sin(Thetam(1) + Thetam(2) + Thetam(3))
    J = [-L1*S1 - L2*S12 - L3*S123m, -L2*S12 - L3*S123m, -L3*S123m ;  L1*C1 + L2*C12 + L3*C123m, L2*C12 + L3*C123m, L3*C123m  ]
    Thetamp= pinv(J)*[Pxp ,Pyp]'
    Thetakm=[Thetam(1)+Thetamp(1), Thetam(2)+Thetamp(2), Thetam(3)+Thetamp(3), 0]
    theta3 = angle(Thetakm(3) - %pi/2 + Gamma)
    theta4 = a*theta3 + b
    Thetakm(4) = angle(theta4 + Gamma)
    theta4mSymetrie = 2*Beta - 2*sum(Thetakm(1:3)) - Thetakm(4)
    //Thetakm(4) = theta4mSymetrie
endfunction

function bool=membreCasse(Thetam, alpham, Gamma)
    //alpha = alpham - %pi/2 + Gamma
    epauleMax = %pi/2// + (%pi/2-alpha)
    epauleMin = -140/180*%pi// + (%pi/2-alpha)
    coudeMax = 145/180*%pi
    coudeMin = 0
    poigneMaxm = 15/180*%pi + %pi/2 - Gamma
    poigneMinm = -45/180*%pi + %pi/2 - Gamma
    bool=(Thetam(1)>epauleMax | Thetam(1)<epauleMin | Thetam(2)>coudeMax | Thetam(2)<coudeMin | Thetam(3)>poigneMaxm | Thetam(3)<poigneMinm)
endfunction

function Thetam=ModeleGeometriqueInverse(Beta, Thetam0,dt,L1,L2,L3,a,b,Gamma, A)
    i = 1
    Thetam = matrix(Thetam0, 1, -1)
    pasfini = %T
    while i<1000 && pasfini
        Pxp = A*cos(Beta)*i*dt
        Pyp = A*sin(Beta)*i*dt
        i = i+1
        thetam = Modelegeo(Pxp,Pyp,L1,L2,L3,a,b,Gamma,Beta,Thetam(i-1,:),dt)
        alpha = angle(sum(thetam(1:3)) + Gamma)
        if  alpha - %pi/2 - Beta <= 0 then //thetam(3) - %pi/2 + Gamma >= 0 then
            pasfini = %F
            Thetam(i, :) = thetam
            disp("fini")
        else 
            Thetam(i, :) = thetam
            pasfini = ~membreCasse(Thetam(i, :), alpha, Gamma)
        end
    end
endfunction

Thetam=ModeleGeometriqueInverse(Beta, Theta0m, dt, L1, L2, L3, a, b, Gamma, 1)
//Theta=ModeleGeometriqueInverse( Beta, Theta(length(Theta(:, 1)),:),dt,L1,L2,L3,a,b,Gamma, 1)


// Paramètres d'affichage

n = length(Thetam(:, 1))
disp(Thetam(n, 3) - %pi/2 + Gamma)

Lparamm = list(init_param('N', 4, 'r1', L1, 'theta1', Thetam(1, 1), 'r2', L2, 'theta2', Thetam(1, 2), 'r3', L3, 'theta3', Thetam(1, 3), 'r4', Lj, 'theta4', Thetam(1, 4)))

for i=2:n
    Lparamm(i) = init_param('N', 4, 'r1', L1, 'theta1', Thetam(i, 1), 'r2', L2, 'theta2', Thetam(i, 2), 'r3', L3, 'theta3', Thetam(i, 3), 'r4', Lj, 'theta4', Thetam(i, 4))
end

Lparam = list(init_param('N', 5, 'r1', L1, 'theta1', Thetam(1, 1), 'r2', L2, 'theta2', Thetam(1, 2), 'r3', lH, 'theta3', Thetam(1, 3) - %pi/2 + Gamma, 'r4', lA,'theta4', %pi/2, 'r5', Lj, 'theta5', Thetam(1, 4) - Gamma))

for i=2:n
    Lparam(i) = init_param('N', 5, 'r1', L1, 'theta1', Thetam(i, 1), 'r2', L2, 'theta2', Thetam(i, 2), 'r3', lH, 'theta3', Thetam(i, 3) - %pi/2 + Gamma, 'r4', lA,'theta4', %pi/2, 'r5', Lj, 'theta5', Thetam(i, 4) - Gamma)
end

//for i=n:(n+100)
//    P = P0+ [cos(Beta)*(i*dt)^2, sin(Beta)*(i*dt)^2]
//    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(n, 1), 'r2', L2, 'theta2', Theta(n, 1)+Theta(n, 2), 'r3', L3, 'theta3', Theta(n, 1)+Theta(n, 2)+Theta(n, 3) + %pi/2 - Gamma, 'Javelot', %T, 'xJavelot', P(1), 'yJavelot', P(2), 'LJavelot', Lj, 'thetaJavelot', Beta)
//end

//Plot(Propulseur, Lparam, dt, L1+L2+L3, %f)
PlotAndSave(Propulseur, Lparam, 'testprop.gif', dt*4, 10, L1+L3+L2, %f)
//Propulseur(Lparam(2))
