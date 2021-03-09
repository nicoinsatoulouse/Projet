function rtheta=CartToCyl(xy)
    // Cartésiennes à polaires
    theta = atan(xy(2), xy(1))
    r = norm(xy)
    rtheta = [r, theta]
endfunction

function xy=CylToCart(rtheta)
    // Polaires à cartésiennes
    x = rtheta(1)*cos(rtheta(2))
    y = rtheta(1)*sin(rtheta(2))
    xy = [x, y]
endfunction
