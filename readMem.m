fileID = fopen('InpOut.bin', 'r');
mat = fread(fileID);
fclose(fileID);

oplen = 500;
signalin = zeros(1, oplen);

for i=1:oplen
    if mat(2*i) > 127
        signalin(1, i) = mat(2*i)-256;
        signalin(1, i) = signalin(1, i) + mat(2*i-1)/256;
    else
        signalin(1, i) = mat(2*i);
        signalin(1, i) = signalin(1, i) + mat(2*i-1)/256;
    end
end

figure(1)
hold on;
freqz(signalin);
title('Input frequency response');
hold off;

fileID = fopen('FilOut.bin', 'r');
mat = fread(fileID);
fclose(fileID);

oplen = 33;
signalfil = zeros(1, oplen);

for i=1:oplen
    if mat(2*i) > 127
        signalfil(1, i) = mat(2*i)-256;
        signalfil(1, i) = signalfil(1, i) + mat(2*i-1)/256;
    else
        signalfil(1, i) = mat(2*i);
        signalfil(1, i) = signalfil(1, i) + mat(2*i-1)/256;
    end
end

figure(2)
hold on;
freqz(signalfil);
title('Filter response');
hold off;

fileID = fopen('MemOut.bin', 'r');
mat = fread(fileID);
fclose(fileID);

oplen = 532;
signal = zeros(1, oplen);

for i=1:oplen
    if mat(2*i) > 127
        signal(1, i) = mat(2*i)-256;
        signal(1, i) = signal(1, i) + mat(2*i-1)/256;
    else
        signal(1, i) = mat(2*i);
        signal(1, i) = signal(1, i) + mat(2*i-1)/256;
    end
end

figure(3)
hold on;
freqz(signal);
title('Output frequency response');
hold off;

figure(4)
hold on;
periodogram(signalin);
title('Input PSD');
hold off;


figure(5)
hold on;
periodogram(signalfil);
title('Filter PSD');
hold off;

figure(6)
hold on;
periodogram(signal);
title('Output PSD');
hold off;
