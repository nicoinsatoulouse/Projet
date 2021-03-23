function bool=membreCasse(Theta, alpha)
    epauleMax = %pi/2 + (%pi/2-alpha)
    epauleMin = -140/180*%pi + (%pi/2-alpha)
    coudeMax = 145/180*%pi
    coudeMin = 0
    poigneMax = 15/180*%pi
    poigneMin = -45/180*%pi
    bool=(Theta(1)>epauleMax | Theta(1)<epauleMin | Theta(2)>coudeMax | Theta(2)<coudeMin | Theta(3)>poigneMax | Theta(3)<poigneMin)
endfunction

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

function Pf=rotation(Pi, a)
    A = [cos(a), -sin(a); sin(a), cos(a)]
    Pf = (A*Pi')'
endfunction

function x=ptInit(a, b, alpha, L1, L2, L3)
    rand("uniform")
    casse = %T
    i=0
    while i<10000 && casse
        lamb = rand()
        x = lamb*a + (1-lamb)*b
        theta = Modelegeo(alpha, x(1), x(2), L1, L2, L3)
        if theta==%F then
            casse = %T
        else
            casse = membreCasse(theta, alpha)
        end
        i = i+1
    end
endfunction

function c=butee(ptOK, ptPasOK, alpha, L1, L2, L3, prec)
    a = ptOK
    b = ptPasOK
    while norm(b-a)>prec
        c = (a+b)/2
        theta = Modelegeo(%pi/2, c(1), c(2), L1, L2, L3)
        if theta==%F then
            casse = %T
        else
            casse = membreCasse(theta, alpha)
        end
        if casse then
            b = c
        else
            a = c
        end
    end
    c = (a+b)/2
endfunction

function l=longueur(r, alpha, L1, L2, L3, prec)
    tmp = sqrt((L1+L2+L3)^2-r^2)
    maxd = [tmp, r]
    maxg = [-tmp, r]
    P = ptInit(maxg, maxd, alpha, L1, L2, L3)
    P = rotation(P, %pi/2-alpha)
    a = butee(P, maxg, alpha, L1, L2, L3, prec)
    b = butee(P, maxd, alpha, L1, L2, L3, prec)
    l = norm(b-a)
endfunction

L1 = 0.286
L2 = 0.279
L3 = 0.067
alpha = 0.1:0.1:%pi/2-0.1
r = 0:0.1:(L1+L2+L3)
distance = []
for i=1:length(alpha)
    d = []
    for j=1:length(r)
        d(j) = longueur(r(j), alpha(i), L1, L2, L3, 0.1)
    end
    distance(i) = max(d)
end
plot(r, d)
//plot(distance, alpha)

disp(butee(ptOK, [], %pi/4, L1, L2, L3, 0.1))
