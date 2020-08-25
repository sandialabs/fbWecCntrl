clear
close all

cf = 60;

mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Fe = sqrt(8*real(Zi))*1; % const. power excitation

%%

Zpto = PTO_Impedance(w,[1, 0, 0, 0, 1, 1e-3, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]

% CC on hydro
legCel{1} = 'CC on hydro';
C{1} = conj(Zi);
ZL{1} = Load_impedance(Zpto,C{1});

% CC on full sys. (eq. 22), output matching condition only
legCel{2} = 'CC on full sys.';
C{2} = conj( squeeze(Zpto(2,2,:)) ...
    - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
    ./ (squeeze(Zpto(1,1,:)) + Zi) );
ZL{2} = C{2};


Pmax = abs(Fe).^2 ./ (8*real(Zi));


for ii = 1:length(C)
    Zin{ii} = input_impedance(Zpto,ZL{ii});
    Pmech(:,ii) = oneDof_mech_power(Zi, Zin{ii}, Fe);
    [~,Pelec(:,ii)] = oneDof_PTO_power(ZL{ii},Zpto,Zi,Fe);
end



fig = figure;
fig.Position = fig.Position .* [1 1 1 0.5];
% set(gca,'yscale','log')
hold on
grid on

for ii = 1:length(C)
    plt(1,ii) = plot(f,Pmech(:,ii),'--','LineWidth',1.5);
end
ax = gca; ax.ColorOrderIndex = 1;
for ii = 1:length(C)
    plt(2,ii) = plot(f,Pelec(:,ii),'-','LineWidth',1.5);
end

l1 = legend([plt(2,:)],[legCel(:)']);
set(l1,'location','southwest')
ylim([-5,1])
xlim([0.2, 1])

ylabel('Efficiency [ ]')
xlabel('Frequency [Hz]')
exportgraphics(gcf,'codesign_freqPowerComp.pdf','ContentType','vector')
