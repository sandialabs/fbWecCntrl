function Pmech = oneDof_mech_power(Zi, Zin, Fe)
    Pmech = 1/2 * abs( Fe ./ (Zi + Zin) ).^2 .* real(Zin);
end