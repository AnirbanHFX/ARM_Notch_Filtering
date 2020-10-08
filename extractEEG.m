idx = 1;
channel = 1;
len = 500;
wordlen = 16;
fraclen = 8;

inp = load('eegdata.mat');
mat = inp.data;
signal = mat{idx}{4}(channel,1:len);

signalfixed = fi(signal, 1, wordlen, fraclen);

fileID = fopen('extractedEEG.txt', 'w');

for i=1:length(signalfixed)
    fp = signalfixed(1, i);
    if i < length(signalfixed)
        fprintf(fileID, '0x%s, ', fp.hex);
    else
        fprintf(fileID, '0x%s', fp.hex);
    end
end

fclose(fileID);
