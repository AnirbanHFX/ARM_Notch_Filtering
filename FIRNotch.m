% Owner
% Anirban Majumder
% Git : https://github.com/AnirbanHFX
% Provided as is
%%
clear all;
%%
% Parameters of the filter
Fs = 250;
f0 = 50;
alpha = 1;
M = 32;
N = floor(M/2+1);
inplen = 500;

Wt = 200;
epsilon = 0;

w0 = 2*pi*f0/Fs;
a = 2*pi*alpha/Fs;

%%
% Calculation of Q and P matrices
% Q = integral(UU^t)|w=(0,w0-a) + integral(Wt*UU^t)|w=(w0-a,w0+a) + integral(UU^t)|w=(w0+a,pi)
% P = integral(-2U)|w=(0,w0-a) + integral(Wt*epsilon*(-2U))|w=(w0-a,w0+a) + integral(-2U)|w=(w0+a,pi)

Q = zeros(N, N);
for i=1:N           % Calculation of symmetric Q matrix
    for j=1:N
        Q(i, j) = integral(@(x) cos((i-1)*x).*cos((j-1)*x), 0, w0-a) + integral(@(x) Wt*cos((i-1)*x).*cos((j-1)*x), w0-a, w0+a) + integral(@(x) cos((i-1)*x).*cos((j-1)*x), w0+a, pi);
    end
end

P = zeros(N, 1);    
for i=1:N           % Calculation of P matrix
    P(i, 1) = -2*integral(@(x) cos((i-1)*x), 0, w0-a) - 2*epsilon*Wt*integral(@(x) cos((i-1)*x), w0-a, w0+a) - 2*integral(@(x) cos((i-1)*x), w0+a, pi);
end
%%
% Fixed point parameters

wordlen = 16; 
fraclen = 8;
signed = 1;

%%
% Determine filter coefficients using Particle Swarm Optimization

% PSO parameters
args.wl = wordlen;
args.fl = fraclen;
args.Q = Q;
args.P = P;

args.iMax = 400;
args.PopSize = 200;
args.c1 = 2;
args.c2 = 2;
args.b = 0.99;
args.w = 1;

args.dim = 2;                   % Particles modelled as a column vector [a1 b1; a2 b2; a3 b3; ...] where a = mantissa (signed 8 bit integer), b = fractional part (unsigned 8 bit integer)
args.order = M/2+1;             % Number of rows for each particle
range.maxP = [];
range.minP = [];
args.range = repmat(range, args.dim, 1);
args.range(1).maxP = 2^(wordlen-fraclen-1)-1;
args.range(1).minP = -2^(wordlen-fraclen-1);
args.range(2).maxP = 2^(fraclen)-1;
args.range(2).minP = 0;

% Execute PSO
PSO_Out = PSO_fp(args);    % X calculated using PSO
while(PSO_Out.cost > -2.6)      % Retry if PSO fails to converge
    PSO_Out.cost
    disp("failed to converge, retrying")
    PSO_Out = PSO_fp(args);    % X calculated using PSO
end
X_PSO = PSO_Out.X;

%%
% Calculate X using Gaussian inversion of Q

X = -0.5*(Q\P);          % X calculated using Gaussian elimination

%%
% Calculate filter coefficients from X and X_PSO

h = zeros(1, M+1);  
h(1, N) = X(1, 1);
for k=1:M/2
    h(1, N-k) = 0.5*X(k+1, 1);
    h(1, N+k) = 0.5*X(k+1, 1);              % Filter coefficients from X
end

figure(1)       % Response of Gaussian Inversion filter
hold on;
freqz(h);       
title('Fitler using Gaussian Inversion');
hold off;

h_PSO = zeros(1, M+1);
h_PSO(1, N) = X_PSO(1, 1);
for k=1:M/2
    h_PSO(1, N-k) = 0.5*X_PSO(k+1, 1);
    h_PSO(1, N+k) = 0.5*X_PSO(k+1, 1);      % Filter coefficients from PSO
end

h_d = double(h_PSO);

figure(2)       % Response of PSO filter
hold on
freqz(h_d);
title('Filter using PSO');
hold off

h_PSO = fi(h_PSO, signed, wordlen, fraclen);

h_fx = fi(h, signed, wordlen, fraclen);

%%
% Load data from EEG signal

load('eegdata.mat')

IP = data{1}{4}(1,1:inplen);    % EEG Signal

IP_fx = fi(IP, signed, wordlen, fraclen);

%%
% Convolve input with h and h_PSO

OP_fx = fi(zeros(1, inplen+M), signed, wordlen, fraclen);

for i=1:inplen+M
    for j=1:M+1
        if i-j+1 >= 1 & i-j+1 <= inplen
            OP_fx(1, i) = OP_fx(1, i) + IP_fx(1, i-j+1)*h_fx(j);    % Convolution
        end
    end
end

OP_fx_PSO = fi(zeros(1, inplen+M), signed, wordlen, fraclen);

for i=1:inplen+M
    for j=1:M+1
        if i-j+1 >= 1 & i-j+1 <= inplen
            OP_fx_PSO(1, i) = OP_fx_PSO(1, i) + IP_fx(1, i-j+1)*h_PSO(j);
        end
    end
end

%%
% Write fixed point filter coefficients to external text file for both filters (Gaussian elimination and PSO)

file1 = fopen('filter_gauss.txt', 'w');
for i=1:length(h_fx)
    fp = h_fx(1, i);
    if i < length(h_fx)
        fprintf(file1, '0x%s, ', fp.hex);
    else
        fprintf(file1, '0x%s', fp.hex);
    end
end
fclose(file1);

file2 = fopen('filter_pso.txt', 'w');
for i=1:length(h_PSO)
    fp = fi(h_PSO(1, i), signed, wordlen, fraclen);
    if i < length(h_PSO)
        fprintf(file2, '0x%s, ', fp.hex);
    else
        fprintf(file2, '0x%s', fp.hex);
    end
end
fclose(file2);

%%
% Plot frequency response of input and output with both filters

figure(3)
hold on;
freqz(IP);
title('Input signal');  
hold off;

temp = double(OP_fx);
figure(4)
hold on;
freqz(temp);
title('Output signal using Gaussian filter');
hold off;

temp = double(OP_fx_PSO);
figure(5)
hold on;
freqz(temp);
title('Output signal using PSO filter');
hold off;

%%



