%DAS_BF
clear;
close all;
load('Computed_RIRs.mat');

%==================Generate array signals==================================%
speechfilename = {'wav/6319-275224-0008.flac', 'wav/6319-275224-0011.flac'};
noisefilename = {'wav/noise1.wav', 'wav/noise2.wav'};
arrayname = {'wav/array4-60.wav'};
% [source1, fs] = audioread(speechfilename{1});
% [source2, ~] = audioread(noisefilename{1});
% [noise, fs_n] = audioread(noisefilename{1});

[array, fs] = audioread(arrayname{1});

% noise = resample(noise,fs,fs_n);
% 
% n_f = fs * 10; %2 seoconds
% source1 = source1(1:n_f);
% source2 = source2(1:n_f) ./ 3;
% noise = noise(1:n_f);
% 
% rir = RIR_sources(:,:,1);
% speech1 = fftfilt(rir, source1).*30;
% 
% rir = RIR_sources(:,:,2);
% speech2 = fftfilt(rir,source2).*30;
% 
% array = speech1 + repmat(noise, 1, 5);
% audiowrite('output/output.wav', array, fs);
%==========================================================================%

%=============================STFT==============================%

frame_length = 480;
frame_shift = 480;
fft_len = 1024;
d = 0.048;%0.048
c = 340;
[nFrames, ffts] = arrayStft(array, frame_length, frame_shift, fft_len);
%[frames, ffts] = multi_fft(arraySignal, frame_length, frame_shift, fft_len);
%nFrames = size(frames, 2)
figure(1);
subplot(8,1,1);
plotSpectrogram((1:nFrames)*frame_shift/fs, (1:fft_len/2)*fs/fft_len, squeeze(ffts(1,:,:)));
phi = zeros(nFrames,1);
phi2 = zeros(nFrames,1);
phi3 = zeros(nFrames,1);
phi4 = zeros(nFrames,1);
phi5 = zeros(nFrames,1);
phi6 = zeros(nFrames,1);
tmp1 = 90;
tmp2 = 90;
tmp3 = 90;
tmp4 = 90;
tmp5 = 90;
xcorr_tmp1 = zeros(fft_len,1);
xcorr_tmp2 = zeros(fft_len,1);
xcorr_tmp3 = zeros(fft_len,1);
xcorr_tmp4 = zeros(fft_len,1);
xcorr_tmp5 = zeros(fft_len,1);

for i=1:1:nFrames
    % abs
    [phi(i)] = gcc_phat(ffts(1:1:2,i,:), fs, d, c);
    if real(phi(i)) <0 
       phi(i) = 0;
    end
    if i >= 415
        i = i;
    end
    % 1-3 tdoa
    [phi2(i)] = tdoa(array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),:), fs, d, c, tmp1,xcorr_tmp1);
    tmp1 = phi2(i);
    if real(phi2(i)) <0 
       phi2(i) = 0;
    end
    % 1-4 gcc
    [r,phi3(i),xcorr_tmp2] = gccphat3(array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),2),array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),3), fs, d*1, c, tmp2,xcorr_tmp2);
    tmp2 = phi3(i);
    if real(phi3(i)) <0 
       phi3(i) = 0;
    end
    % 1-3 freq
    [r2,phi4(i),xcorr_tmp3] = gccphat4(array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),1),array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),3), fs, d*2, c,tmp3,xcorr_tmp3);
    tmp3 = phi4(i);
    if real(phi4(i)) <0 
       phi4(i) = 0;
    end
    % 1-4 freq
    [r3,phi5(i),xcorr_tmp4] = gccphat4(array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),1),array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),4), fs, d*3, c,tmp4,xcorr_tmp4);
    tmp4 = phi5(i);
    if real(phi5(i)) <0 
       phi5(i) = 0;
    end
    % 2-4 freq
    [r4,phi6(i),xcorr_tmp5] = gccphat4(array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),2),array(frame_shift*(i-1)+1:1:frame_length+frame_shift*(i-1),4), fs, d*2, c,tmp5,xcorr_tmp5);
    tmp5 = phi6(i);
    if real(phi6(i)) <0 
       phi6(i) = 0;
    end
end
subplot(8,1,2);
angle1 = mean(real(phi(:)));
plot((1:nFrames)*frame_shift/fs, real(phi(:)));
subplot(8,1,3)
angle2 = mean(real(phi2(:)));
plot((1:nFrames)*frame_shift/fs, real(phi2(:)));
subplot(8,1,4)
angle3 = mean(real(phi3(:)));
plot((1:nFrames)*frame_shift/fs, real(phi3(:)));
subplot(8,1,5)
angle4 = mean(real(phi4(:)));
plot((1:nFrames)*frame_shift/fs, real(phi4(:)));
subplot(8,1,6)
angle5 = mean(real(phi4(:)));
plot((1:nFrames)*frame_shift/fs, real(phi5(:)));
subplot(8,1,7)
angle6 = mean(real(phi4(:)));
plot((1:nFrames)*frame_shift/fs, real(phi6(:)));
subplot(8,1,8)
plot((1:nFrames)*frame_shift/fs, 0.5*real(phi6(:))+ 0.5*real(phi4(:)));
% [nFrames, spec_s] = arrayStft(speech1(:, 1), frame_length, frame_shift, fft_len);
% subplot(3,1,2);
% plotSpectrogram((1:nFrames)*frame_length/fs, (1:fft_len/2)*fs/fft_len, squeeze(spec_s(1,:,:)));
%==========================================================================%

%=============================apply CGMM=============================%
% lambda_y = zeros(size(ffts, 2), size(ffts, 3));
% lambda_n = zeros(size(ffts, 2), size(ffts, 3));
% 
% mini_batch = size(ffts,2); %estimate cgmm per  mini_batch 
% iters = 10;
% for t = 1:mini_batch:size(ffts, 2)
%     idx =  t:min(t+mini_batch-1, size(ffts,2));
%     [lambda_y_one, lambda_n_one,  Ry, Rn, Q] = cgmm_em(ffts(:, idx, :), iters);
%     lambda_y(idx, :) = lambda_y_one;
%     lambda_n(idx, :) = lambda_n_one;
% end
%==========================================================================%
% figure(2);
% subplot(2,1,1);
% plotMask(lambda_y.');
% subplot(2,1,2);
% plotMask(lambda_n.');
% Rx = Ry -Rn;      

% [M, T, F]  = size(ffts); %fft bins number
% d = zeros(M, F);         %steering vectors
% hw = d;                   %mvdr beamforming weight 
% output = zeros(T, F);    %beamforming outputs

%========== Steering vectors esimation  ===========%
% for f= 1:F
%     d(:, f) = rtf(squeeze(Ry(:, :, f)), squeeze(Rn(:, :, f)), 'EVD');
%     hw(:, f) = mvdr(d, Rn(:, :, f));
%     output(:, f) =   hw(:, f) .'* squeeze(ffts(:, :, f));
% end
% figure(1);
% subplot(3,1,3)
% plotSpectrogram((1:nFrames)*frame_length/fs, (1:fft_len/2)*fs/fft_len, squeeze(output(1,:,:)));

%====================Quality evaluation====================================%
%output = [output, fliplr(conj(output(:, 2:end-1)))];
%rec_frames = real(ifft(output, fft_len, 2));
%rec_frames = rec_frames(:,1:frame_length);
%sig = overlapadd(rec_frames, hamming(frame_length, 'periodic'), frame_shift);

% sig = invStft(output, frame_length, frame_shift);
% audiowrite('output/output_cgmm.wav', sig./max(abs(sig)), fs);
% save('output/testCGMM.mat')

%==========================================================================%