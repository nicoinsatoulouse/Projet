exec(pwd()+"\CartCyl.sci")

function res=PlotAndSave(obj, Lparamobj, outgif, dt, fps, tailleMax, persistent)
    if ~exists("persistent","local") then
        persistent = %f
    end
    fig = figure
    fig.background = 8
    h = obj(Lparamobj(1))
    isoview on
    h_axes = gca()
    h_axes.data_bounds = [-1.1*tailleMax, -1.1*tailleMax; 1.1*tailleMax, 1.1*tailleMax]
    h_axes.tight_limits = "on"
    
    delay = round(1000/fps - modulo(1000/fps, 10)) //Pour que l'échelle temporelle soit OK dans le fichier enregistré
    mdelete(outgif);
    idGif = animaGIF(gcf(), outgif, delay, 2)
    
    pas = delay/(1000*dt) //Same shit d'échelle temporelle
    for i=2:pas:length(Lparamobj)
        drawlater()
        if ~persistent then
            delete(h)
        else
            if is_param(Lparamobj(i), 'Javelot') && get_param(Lparamobj(i), 'Javelot')(1)==%T then
                delete(h(length(h)))
            end
        end
        h = obj(Lparamobj(i))
        drawnow()
        idGif = animaGIF(gcf(), idGif)
    end
    
    animaGIF(idGif)
    res = %T
endfunction

function res=Plot(obj, Lparamobj, dt, tailleMax, persistent)
    if ~exists("persistent","local") then
        persistent = %f
    end
    fig = figure
    fig.background = 8
    tempsAct = getdate()(9:10)
    h = obj(Lparamobj(1))
    isoview on
    h_axes = gca()
    h_axes.data_bounds = [-1.1*tailleMax, -1.1*tailleMax; 1.1*tailleMax, 1.1*tailleMax]
    h_axes.tight_limits = "on"
    
    for i=2:length(Lparamobj)
        drawlater()
        if ~persistent then
            delete(h)
        else
            if is_param(Lparamobj(i), 'Javelot') && get_param(Lparamobj(i), 'Javelot')(1)==%T then
                delete(h(length(h)))
            end
        end
        tempsAct = Pause(tempsAct, dt) //Pour un resultat en temps réel (hors lags)
        h = obj(Lparamobj(i))
        drawnow()
    end
    
    res = %T
endfunction

function tempsAct=Pause(tempsPrec, dt)
    tempsAct = getdate()(9:10)
    temps = max(((dt-tempsAct(1)+tempsPrec(1))*1000-tempsAct(2)+tempsPrec(2)), 0)
    sleep(temps)
endfunction

function h=NTiges(params)
    // Affiche N tiges bout à bout
    N = get_param(params, 'N')(1)
    if is_param(params, 'x0') then
        x0 = get_param(params, 'x0')
    else
        x0 = [0, 0]
    end
    for i=1:N
        r = get_param(params, 'r'+string(i))(1)
        theta = get_param(params, 'theta'+string(i))(1)
        xy = CylToCart([r, theta])
        h(i) = plot([x0(1), x0(1)+xy(1)], [x0(2), x0(2)+xy(2)])
        x0 = x0+xy
    end
    if is_param(params, 'Javelot') && get_param(params, 'Javelot')(1)==%T then
        x = get_param(params, 'xJavelot')(1)
        y = get_param(params, 'yJavelot')(1)
        xy = [x, y]
        L = get_param(params, 'LJavelot')(1)
        theta = get_param(params, 'thetaJavelot')(1)
        dir = CylToCart([1, theta])
        depart = xy - L/2*dir
        fin = xy + L/2*dir
        h(N+1) = plot([depart(1), fin(1)], [depart(2), fin(2)], 'r')
    end
endfunction

function h=Propulseur(params)
    // Affiche N tiges bout à bout dont la dernière est rouge
    N = get_param(params, 'N')(1)
    if is_param(params, 'x0') then
        x0 = get_param(params, 'x0')
    else
        x0 = [0, 0]
    end
    theta = 0
    for i=1:N
        r = get_param(params, 'r'+string(i))(1)
        theta = theta + get_param(params, 'theta'+string(i))(1)
        xy = CylToCart([r, theta])
        if i==N then
            coul = "red"
        elseif i==N-1 then
            coul = "green"
        else
            coul = "blue"
        end
        h(i) = plot([x0(1), x0(1)+xy(1)], [x0(2), x0(2)+xy(2)], coul)
        x0 = x0+xy
    end
endfunction
