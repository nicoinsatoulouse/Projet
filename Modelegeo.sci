function Theta=Modelegeo(alpha,Px,Py)
    l=[0.3,0.3,0.1]
    Wx=Px-l(3)*cos(alpha)
    Wy=Py-l(3)*sin(alpha)
    y=(Wx^2+Wy^2-l(1)^2-l(2)^2)/2*l(1)*l(2)
    T2 = atan(sqrt(1-y^2),y)
    T1 = atan(-l(2)*sin(T2)*Wx+(l(1)+l(2)*cos(T2))*Wy,l(2)*sin(T2)*Wy+(l(1)+l(2)*cos(T2))*Wx)
    T3=alpha-T1-T2
    Theta=[T1,T2,T3]
endfunction
