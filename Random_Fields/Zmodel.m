function [z_model] = Zmodel(N, r0, dr, ph0, dph, Nt, SNR_db, lambda, dt, u_mean, u_sigma, v_mean, v_sigma) 

r = r0 - dr/2 + dr .* rand(1, N);
ph = ph0 - dph/2 + dph .* rand(1, N);

% x0 = 10 + 5 .* rand(1, N);
% % y0 = 10 + 5 .* rand(1, N); 
% y0 = zeros(1, N);

x0 = r .* cos(ph); y0 = r .* sin(ph);


%% Plot positions of scatterers
% txt = ['Resolution cell'];
% dtext = ['Scaterer position at time t = 0'];
% xl = 'x[m]';
% 
% f = figure(107); hold on; f.Position = [10 10 1000 1000];
% color = 'k';
% yl =  ['y[m]'];
% 
% marker = markers(2);
% 
% plott2(x0, y0, xl, yl, txt, 2, dtext, color, marker)
%%
% figure(101); plot(x0, y0, '*');

x(1, :) = x0;
y(1, :) = y0;

% v = linspace(1, 5, N);
u = normrnd(u_mean, u_sigma, [1 N]);
v = normrnd(v_mean, v_sigma, [1 N]);

A0 = ones(1, N);
D(1, :) = sqrt(x(1, :).^2 + y(1, :).^2);

z(1) = sum(A0 .* exp(1j .* 4 * pi / lambda * x0));

SNR = 10^(SNR_db/10);

for i = 2:Nt
    x(i, :) = x(i - 1, :) + u .* dt;
    y(i, :) = y(i - 1, :) + v .* dt;
    
%     figure(101); hold on; plot(x(i, :), y(i, :), '+');
    
    D(i, :) = sqrt(x(i, :).^2 + y(i, :).^2);
    z(i) = sum(A0 .* exp(1j .* 4 * pi / lambda .* D(i, :)));
end

Noise = sum(abs(z).^2)./(Nt .* SNR);
sigma_n = sqrt(Noise);

z_model = z + sigma_n .* (randn(1, Nt));

end