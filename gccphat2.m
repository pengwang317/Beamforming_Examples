function [r,tau] = gccphat2( x, y, fs,d,c )
M = max(numel(x),numel(y));
if sum(abs(x(:))) < 5 && sum(abs(y(:)))< 5
    tau = 90;
    r=0;
else
%%Transform both vectors
% X = fft(x,2^nextpow2(2*M-1));
% Y = fft(y,2^nextpow2(2*M-1));
% 
% % Compute cross-correlation
% 
% R = X.*conj(Y);
% c = ifft(R./abs(R));

%%
N = 2*M-1; 
Nfft = 2^nextpow2(N);

R = bsxfun(@times, ...
        fft(y,Nfft), ...
        conj(fft(x,Nfft)));
% rtmp = fftshift( ...
%         ifft(exp(1i*angle(R))) ,1);

tmp = fft(y,Nfft).*conj(fft(x,Nfft));
% r2 = ifft(tmp./(abs(tmp)+0.00001));
rtmp = fftshift(ifft(tmp));

r = rtmp(Nfft/2+1-(N-1)/2:Nfft/2+1+(N-1)/2,:);

lags = (-(N-1)/2:(N-1)/2).';
lags = lags/fs;
[k,idx] = max(abs(r));
tau = lags(idx);
if k < 0.1 || tau <= -8 || tau >= 8
    tau = 90;
else
% tau = N/(2*fs)+lags(idx);

tmp_cos = tau*c/d;
if tmp_cos > 1
    tmp_cos = 1;
elseif tmp_cos < -1
    tmp_cos = -1;
end
tau = acos(tmp_cos) * 180 / pi;

end
end
end