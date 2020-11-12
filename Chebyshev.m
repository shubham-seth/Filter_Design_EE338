%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))]; % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Band Edge speifications
fs1 = 57.6;
fp1 = 61.6;
fp2 = 81.6;
fs2 = 85.6;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%analog lpf
[num, dem] = numden(analog_lpf(s));
num = sym2poly(expand(num));                          
dem = sym2poly(expand(dem));                          %collect coeffs into matrix form
k = dem(1);    
dem = dem/k;
num = num/k;
fvtool(num,dem)

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
% fvtool(nz,dz)                                       %frequency response in dB


%magnitude plot (not in log scale)
[H,f] = freqs(ns, ds,10e3);
figure,plot(f,abs(H))
grid
title("FT of Band-Pass Filter in Analog doomain")
ylabel("Magnitude response")
xlim([0 2])

figure,plot(f,angle(H))
grid
title("FT of Band-Pass Filter in Analog doomain")
ylabel("phase response")
xlim([0 2])

% DTFT
f = -pi:3.1416e-04:pi;
[H] = freqz(ns,ds,f, 2*pi);
figure,plot(f,abs(H))
grid
title("DTFT of Band-Pass Filter")
ylabel("Magnitude response")
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2', '3\pi/4','\pi'})
xlim([-pi pi])

figure,plot(f,angle(H))
grid
title("DTFT of Band-Pass Filter")
ylabel("Phase response")
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2', '3\pi/4','\pi'})
xlim([-pi pi])
