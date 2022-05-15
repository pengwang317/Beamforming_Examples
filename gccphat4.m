function [r,tau,xcorr_tmp_out] = gccphat4( x, y, fs,d,c,tmppha,xcorr_tmp )
M = max(numel(x),numel(y));
a = sum(abs(x(:)));
b = sum(abs(y(:)));
if sum(abs(x(:))) < 200/32767/M && sum(abs(y(:)))< 200/32767/M
    tau = tmppha;
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
% N = length(x1)+length(x2)-1;
% NFFT = 64;
% range = NFFT/2+1 - (N-1)/2:NFFT/2+1 + (N-1)/2;
% r1 = xcorr(x1,x2);
% r2 = ifft(fft(x1,NFFT).*conj(fft(x2,NFFT)));
% r3 = fftshift(r2);
% r4 = r3(range);

N = 2*M-1; 
Nfft = 2^nextpow2(N);

% R = bsxfun(@times, ...
%         fft(x,Nfft), ...
%         conj(fft(y,Nfft)));
% rtmp = fftshift( ...
%         ifft(exp(1i*angle(R))) ,1);
tmp = fft(y,Nfft).*conj(fft(x,Nfft));
tmp_buffer = 0.9*xcorr_tmp + 0.1*tmp;
xcorr_tmp_out = tmp_buffer;
% r2 = ifft(tmp./(abs(tmp)+0.00001));
tmp2 = angle(tmp_buffer);
tmp3 = tmp./(abs(tmp));
% r2 = ifft(exp(1i*angle(tmp)));
r2 = ifft(tmp_buffer);
rtmp = fftshift(r2);
range = Nfft/2+1-(N-1)/2:Nfft/2+1+(N-1)/2;
r = rtmp(range);

lags = (-(N-1)/2:(N-1)/2).';
% lags = lags/fs;
[k,idx] = max(abs(r));
tau0 = lags(idx);
threshold = floor(d/c*fs);
if k < 0.15 || tau0 < -threshold || tau0 > threshold
    tau = tmppha;
else
    % tau = N/(2*fs)+lags(idx);

    tmp_cos = tau0/fs*c/d;
    if tmp_cos > 1
        tmp_cos = 1;
    elseif tmp_cos < -1
        tmp_cos = -1;
    end

    tau = 0.5*tmppha + 0.5*acos(tmp_cos) * 180 / pi;
end

if tau >170 || tau<10
    tau = tau;
end

end
end