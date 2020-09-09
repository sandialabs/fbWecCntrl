function [ZL] = Zi2ZL(Zpto, C)

ZL = ((squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ./ ...
    ( squeeze(Zpto(1,1,:)) - C )) - squeeze(Zpto(2,2,:)));
