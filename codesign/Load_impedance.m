function ZL = Load_impedance(w, x, C)


% can call the controller from here. at the momem we assume it is just a PI

ZL = x(1) - 1i*x(2)./w(:);


