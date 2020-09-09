function C = fbc(x,cinfo)
    
    switch cinfo.type
        case 'PI'
            C = x(1) - 1i*x(2)./cinfo.w(:);
        case 'P'
            C = x(1);
        case 'CC'
            C = conj(cinfo.Zi);
        otherwise
            error('Invalid value for ''cinfo.type''')
    end
end