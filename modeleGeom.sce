exec(pwd()+"\resolutionEDO.sci")
exec(pwd()+"\PlotSave.sci")
exec(pwd()+"\CartCyl.sci")

function Theta=Modelegeo(alpha,Px,Py,L1,L2,L3)
    l=[L1,L2,L3]
    Wx=Px-l(3)*cos(alpha)
    Wy=Py-l(3)*sin(alpha)
    y=(Wx^2+Wy^2-l(1)^2-l(2)^2)/(2*l(1)*l(2))
    T2 = atan(sqrt(1-y^2), y)
    T1 = atan(-l(2)*sin(T2)*Wx+(l(1)+l(2)*cos(T2))*Wy, l(2)*sin(T2)*Wy+(l(1)+l(2)*cos(T2))*Wx)
    T3=alpha-T1-T2
    Theta=[T1,T2,T3]
endfunction

function v=V(t)
    v = 0.2
endfunction

function Theta=ModeleGeometriqueInverse(alpha, Px0, Py, v, dt, L1, L2, L3)
    t = 0
    i = 1
    Px = Px0
    theta = Modelegeo(alpha, Px0, Py, L1, L2, L3)
    Theta = matrix(theta, 1, -1)
    vectWP = CylToCart([L3, alpha])
    while norm([Px, Py]-vectWP)+v(t)*dt<=L1+L2
        t = t + dt
        Px = Px0 + t*v(t)
        i = i+1
        Theta(i, :) = Modelegeo(alpha, Px, Py, L1, L2, L3)
    end
endfunction

L1 = 0.3
L2 = 0.3
L3 = 0.1
alpha = %pi/2
Px0 = -0.57
Py = L3 + 0.1
dt = 0.02

Theta = ModeleGeometriqueInverse(alpha, Px0, Py, V, dt, L1, L2, L3)

// Paramètres d'affichage
Lparam = list(init_param('N', 3, 'r1', L1, 'theta1', Theta(1, 1), 'r2', L2, 'theta2', Theta(1, 1)+Theta(1, 2), 'r3', L3, 'theta3', Theta(1, 1)+Theta(1, 2)+Theta(1, 3)))
for i=2:length(Theta(:, 1))
    Lparam(i) = init_param('N', 3, 'r1', L1, 'theta1', Theta(i, 1), 'r2', L2, 'theta2', Theta(i, 1)+Theta(i, 2), 'r3', L3, 'theta3', Theta(i, 1)+Theta(i, 2)+Theta(i, 3))
end

//Plot(NTiges, Lparam, dt, L1+L2+L3)
PlotAndSave(NTiges, Lparam, 'BrasHorizontal.gif', dt, 10, L1+L2, %t)
