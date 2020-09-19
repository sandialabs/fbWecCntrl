clear

ex1 = load('out_var_ex1.mat')
ex2 = load('out_var_ex2.mat')


ex1.out_var_coopt.PTO_param_names

Kp1 = ex1.out_var_coopt.cntrl(1);
Ki1 = ex1.out_var_coopt.cntrl(2);

Kp2 = ex2.out_var_coopt.cntrl(1);
Ki2 = ex2.out_var_coopt.cntrl(2);

w = ex1.out_var_coopt.w;
% Kp1(1) - 1i*x(2)./w(:);


C1 = frd(Kp1(1) - 1i*Ki1./w(:), w);
C2 = frd(Kp1(1) - 1i*Ki1./w(:), w);

G1 = frd(1 ./ ex1.out_var_coopt.Z_th, w);


Nfreq = length(w);

Z_mat = zeros(2,2,Nfreq);

N2 = ex1.out_var_coopt.PTO_param_out(1);
Id2 = ex1.out_var_coopt.PTO_param_out(2);
Bd2 = ex1.out_var_coopt.PTO_param_out(3);
Kd2 = ex1.out_var_coopt.PTO_param_out(4);
Kt2 = ex1.out_var_coopt.PTO_param_out(5);
Zi = ex1.out_var_coopt.Zi;

Zd2 = Bd2 + 1i*(w(:)*Id2 - Kd2./w(:));
G2 = frd(sqrt(3/2) * Kt2 * N2^2 ./ ( N2^2 * Zd2 + Zi ), w);

L1 = G1*C1;
L2 = G2*C2;
S1 = -L1 / (1 + L1);
S2 = -L2 / (1 + L2);


figure(1)
bode(L1, L2)
grid on
legend('colocated', 'non-colocated')

figure(2)
bode(S1, S2)
grid on
legend('colocated', 'non-colocated')




