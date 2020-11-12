f_samp = 330e3;

%Band Edge speifications
fs1 = 57.6e3;
fp1 = 61.6e3;
fp2 = 81.8e3;
fs2 = 85.6e3;
ws1 = fs1*2*pi/f_samp;
wp1 = fp1*2*pi/f_samp;
wp2 = fp2*2*pi/f_samp;
ws2 = fs2*2*pi/f_samp;

Wc1 = (ws1+wp1)/2;
Wc2 = (ws2+wp2)/2;
w_t = (fp1-fs1)*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*w_t));           %empirical formula for N_min

%Window length for Kaiser Window
N=N_min + 15;
N
%Ideal bandpass impulse response of length "n"
% bp_ideal = ideal_lp(0.494*pi,n) - ideal_lp(0.373*pi,n);
a = (N-1)/2;
n = [0:1:(N-1)];
m = n - a + eps;

bp_ideal = 2 .* cos((Wc2+Wc1).*m./2) .* sin((Wc2-Wc1).*m./2) ./ (pi*m);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(N,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
% fvtool(FIR_BandPass);         %frequency response


%mplotting response
f = 0:3.14159e-04:pi;

h = freqz(FIR_BandPass, 1, length(f));
plot(f,abs(h))
grid
xticks([ws1 wp1 wp2 ws2])
yticks([0.15 0.85 1.15])
title("DTFT of BandPass Filter")
ylabel("Magnitude Response")
xticklabels({'ws1', 'wp1', 'wp2', 'ws2'})

figure,plot(f,angle(h))
grid
title("DTFT of BandPass Filter")
ylabel("Phase Response")
xticks([ws1 wp1 wp2 ws2])
xticklabels({'ws1', 'wp1', 'wp2', 'ws2'})
xlim([0 pi])