clear;

%% Three frequency model

% parameters

N = 64; % Number of observation
K = 3;  % Number of frequencies contained in the signal (considered known in the algorithm)

f = [0.2 0.2+1/N 0.2+2/N]; % [Values of the frequencies]

Ei = [20 6.3246 20]; % [Values of the amplitudes of the harmonics]

argf = [0 pi/4 pi/3]; % [Phase of the harmonics]

ac = Ei.*cos(argf);  % [Cosine components of the harmonics]

as = Ei.*sin(argf);  % [Sin ecomponents of the harmomnics]

SNR_db = 5; % SNR as described in the paper

sn = sqrt( Ei(1)./(2.*exp(SNR_db/10)) ); % Noise std
t = 0:1:N-1; % Time axis
y = zeros(1, N); % signal initialization 

for i = 1:K 

    y = y + ac(i) .* cos(2.*pi.*f(i).*t) + as(i) .* sin(2.*pi.*f(i).*t); % signal model

end

n = normrnd(0, sn, [1 N]);
y = y + n; % Added with standard normal noise 

% -atan(as./ac)/pi

% lam = 0.2;
% sig_RW = 1/(5*N);
% nu_0 = 0; gamma_0 = 0;

%% MCMC 
Iter = 10000;   % Number of iterations

IN = eye(N, N); % Indentitiy matrix 

fk0 = [3 4 15]; % Initial guess of the frequencies

fk(1, :) = fk0;
gamma_0 = 0;    % parameter from paper
lam = 0.2;      % parameter from paper
delta = 0.1;    % parameter from paper - related to SNR
nu_0 = 0;       % parameter from the paper

fki(1).f = [];  % to record the frequencies with iteration for 1st frequency
fki(2).f = [];  % for the second one
fki(3).f = [];  % for the third one

for o = 1:Iter
    
    for m = 1:K
        u = rand();
            if u < lam
                [fk_, t_] = meshgrid(fk, t);
                D = [cos(2 .* pi .* fk_ .* t_) sin(2 .* pi .* fk_ .* t_)];
                Sigk = delta^(-2) * (D')*D;
                M = D'*D+inv(Sigk);
                Pk = IN - D*M*(D');
                pidfkm = (gamma_0 + y*Pk*y').^(-(N+nu_0)/2); % p(w_{j}|z, w_{-j})
                Np = 1024;
                pl = abs((fft(y, Np))).^2;
                pl = pl(1:Np/2);
                freq = linspace(0, 1/2, Np/2);
%                 figure; plot(freq, pl);
                [plmaxf, plmaxf_ind] = max(pl); % If u < lam, Find the maximum of the DFT 
                fk(m) = freq(plmaxf_ind);
                fki(m).f = [fki(m).f fk(m)];
%                 fk(m) = 
            else
                [fk_, t_] = meshgrid(fk, t);
                D = [cos(2 .* pi .* fk_ .* t_) sin(2 .* pi .* fk_ .* t_)];
                Sigk = delta^(-2) * (D')*D;
                M = D'*D+inv(Sigk);
                Pk = IN - D*M*(D');
                pidfkm = (gamma_0 + y*Pk*y').^(-(N+nu_0)/2); % p(w_{j}|z, w_{-j}) 
                fk(m) = normrnd(fk(m), 1/(5*N), 1); % If u > lam, choose a new frequency from a normal dist
                fki(m).f = [fki(m).f fk(m)];
            end

    end
end


for i = 1:K
    figure(i); plot(fki(i).f); hold on; plot(ones(1, Iter).*f(i), '-.', 'Color', 'k');
    legend({'MCMC', 'Ground Truth'})
end
