% KerFt = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
% KerFt = exp(-KerFt.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);

PrefStim = linspace(-Width, Width, NetPars.N+1);
TunWidth  = NetPars.TunWidth;
Width = NetPars.Width;
W = angle(exp(1i*(PrefStim - PrefStim(1)) *pi/Width))* Width/pi;
W = exp(-W.^2/(2*TunWidth^2))/(sqrt(2*pi)*TunWidth);

W = gallery('circul', W);

%%

subplot(1,2,1)
imagesc(PrefStim, PrefStim,W)
colorbar
axis xy; axis square
% cMap = hot(64);
cMap = gray(64);
colormap(flipud(cMap));
set(gca, 'xtick', -180:90:180, 'ytick', -180:90:180)


subplot(1,2,2)
Wrand = randn(20);
imagesc(Wrand);
axis xy; axis square