
function [phi] = gcc_phat(arraySpectrum, fs, d, c)

X1(:) = arraySpectrum(1,1,:);
X2(:) = arraySpectrum(2,1,:);
a2 = sum(abs(X1(:).*X1(:)));
b2 = sum(abs(X2(:)));
% if sum(abs(X1(:))) < 10 && sum(abs(X2(:)))< 10
%     phi = 90;
% else
%range = 1024/2+1 -
G = X2 .* conj(X1);
corr = ifft(G ./ abs(G));

corr2 = fftshift(corr);
[~,I2] = max(abs(corr2));

[~,I] = max(abs(corr));
% tdoa = lag(I);  
tdoa = I2-256; 

phi = acos(tdoa/fs*c/d) * 180 / pi;
phi = a2;
% end
end
