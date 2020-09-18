function Z_mat = PTO_Impedance(Gamma, w)

Nfreq = length(w);

Z_mat = zeros(2,2,Nfreq);

N = Gamma(1);
Id = Gamma(2);
Bd = Gamma(3);
Kd = Gamma(4);
Kt = Gamma(5);
Rw = Gamma(6);
Lw = Gamma(7);

Zd = Bd + 1i*(w(:)*Id - Kd./w(:));
Zw = Rw + 1i*w(:)*Lw;

Z_mat(1,1,:) = Zd * N^2;
Z_mat(1,2,:) = -Kt * N * sqrt(3/2);
Z_mat(2,1,:) = Kt * N * sqrt(3/2);
Z_mat(2,2,:) = Zw;

% Nfreq = length(w);
% 
% Z_mat = zeros(2,2,Nfreq);
% 
% N = PTO_param(1);
% Id = PTO_param(2);
% Bd = PTO_param(3);
% Kd = PTO_param(4);
% Kt = PTO_param(5);
% Rw = PTO_param(6);
% Lw = PTO_param(7);
% 
% 
% Zd = Bd + 1i*(w(:)*Id - Kd./w(:));
% Zw = Rw + 1i*w(:)*Lw;
% 
% Z_mat(1,1,:) = Zd * N^2;
% Z_mat(1,2,:) = -Kt * N * sqrt(3/2);
% Z_mat(2,1,:) = Kt * N * sqrt(3/2);
% Z_mat(2,2,:) = Zw;
% 