function [phi] = tdoa(arraySignals, fs, d, c, tmppha,xcorr_tmp)

x1 = arraySignals(:,1);
x2 = arraySignals(:,2);
a = sum(abs(x1(:)));
b = sum(abs(x2(:)));
if sum(abs(x1(:))) < 100/32767/length(x1) && sum(abs(x2(:)))< 100/32767/length(x1)
    phi = tmppha; % 90
else
[acor,lag] = xcorr(x2,x1);
[k,I] = max(abs(acor));
t = lag(I);   
tmp_cos = t/fs*c/d;
if tmp_cos > 1
    tmp_cos = 1;
elseif tmp_cos < -1
    tmp_cos = -1;
end

threshold = floor(d/c*fs);
if k < 0.15 || t < -threshold || t> threshold
    phi = tmppha; %90
else
    phi = 0.5*tmppha + 0.5*acos(tmp_cos) * 180 / pi;
end

if phi >170 || phi<10
    phi = phi;
end

end
end

function [phi] = gcc_phat(arraySpectrum, fs, d, c)

X1 = arraySpectrum(:,1);
X2 = arraySpectrum(:,2);

G = X1 .* conj(X2);
corr = ifft(G ./ abs(G));

[~,I] = max(abs(corr));
tdoa = lag(I);   

phi = acos(tdoa/fs*c/d) * 180 / pi;

end