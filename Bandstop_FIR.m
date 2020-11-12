f_samp = 260e3;

%Band Edge speifications
fp1 = 45.6e3;
fs1 = 49.6e3;
fs2 = 69.6e3;
fp2 = 73.6e3;

wp1 = fp1*2*pi/f_samp;
wp2  = fp2*2*pi/f_samp;
ws1 = fs1*2*pi/f_samp;
ws2  = fs2*2*pi/f_samp;

wc1 = (ws1+wp1)/2;
wc2 = (ws2+wp2)/2;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*0.04*pi));       %empirical formula for N_min

%Window length for Kaiser Window
N=N_min + 13;
N
%Ideal bandstop impulse response of length "n"
a = (N-1)/2;
n = [0:1:(N-1)];
m = n - a + eps;

bs_ideal = (sin(wc1.*m) ./ (pi*m))  + (2 .* cos((pi+wc2).*m./2) .* (sin((pi-wc2).*m ./2) ./ (pi*m))) ;

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(N,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
% fvtool(FIR_BandStop);         %frequency response
% FIR_BandStop
%magnitude response
f = 0:3.1416e-04:pi;

h = freqz(FIR_BandStop, 1,length(f));
figure,plot(f,abs(h))
grid
xticks([wp1 ws1 ws2 wp2])
yticks([0.15 0.85 1.15])
xticklabels({'wp1', 'ws1', 'ws2', 'wp2'})
title("DTFT of BandStop Filter")
ylabel("Magnitude Response")

figure,plot(f,angle(h))
grid
xticks([ws1 wp1 wp2 ws2])
title("DTFT of BandPass Filter")
ylabel("Phase Response")
xticklabels({'ws1', 'wp1', 'wp2', 'ws2'})