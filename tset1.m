clc;
clear all;
x1 = [1,2,3,7,9,8,3,7]';
x2 = [4,5,6,5,4,3,8,2]';
N = length(x1)+length(x2)-1;
NFFT = 64;
range = NFFT/2+1 - (N-1)/2:NFFT/2+1 + (N-1)/2;
r1 = xcorr(x1,x2);
[k1,ind1] = max(r1);
rtmp = fft(x1,NFFT).*conj(fft(x2,NFFT));
% r2 = ifft(fft(x1,NFFT).*conj(fft(x2,NFFT)));
r2 = ifft(rtmp./abs(rtmp));
r3 = fftshift(r2);
r4 = r3(range);
[k2,ind2] = max(r4);


M = length(x2);
xc = xcorr(x1,x2,'biased');
[k,ind] = max(xc);

fs = 8000;
d = 0.04;
an = acos((ind-M)/fs*340*d)*180/pi;
