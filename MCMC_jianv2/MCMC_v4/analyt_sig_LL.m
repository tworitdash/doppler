clear;
markers = load('markers.mat');
markers = markers.markers;
colors = load('colors.mat');
colors = colors.color;


dt = 1e-3;

lambda = 0.03;

L = 4 * pi / lambda * dt;

Nt = 10; 

k = 0:1:Nt-1;

Nu = 128;

u = linspace(0, lambda/(2 * dt), Nu);

Mu = 7.5;

Sigmavec = linspace(eps, 8, 100);

for s = 1:length(Sigmavec)

Sigma = Sigmavec(s);




for i = 1:Nu
    
%     if s == 1
%         LLt(i) = -(sin((1 - 2 .* Nt) .* L/2 .* (Mu - u(i)))./(sin(L./2 .* (Mu - u(i)))) + 2 .* Nt - 1);
%     end

    LL(i) = -sum(1 + exp(-k.^2 .* L.^2 .* Sigma.^2) - ...
    2 .* exp(-k.^2 .* L.^2 .* Sigma.^2/2) .* cos(L .* k .* (Mu - u(i))));

  

end



LLo = LL - min(LL);
% 
  if s == 1
        figure; plot(u, LLo); % hold on; plot(u, LLt);
  end

dv = u(2) - u(1);

dvres = lambda/(2 * dt)/(Nt - 1);

PT(s) = sum(abs(LLo).^2  .* dv);
mu(s) = 1./PT(s) .* sum(u .* abs(LLo).^2 .* dv);
sigma(s) = 2 .* sqrt(sum(1./PT(s) .* (u - mu(s)).^2 .* abs(LLo).^2 .* dv));


end

txt = ['\sigma_{True} vs \sigma_{re}'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', SNR = ', num2str(SNR_db), ' [dB]'];
% dtext = ['Nt = ', num2str(Nt), ', dv = ', num2str(dv), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]'];

dtext = ['Nt = ', num2str(Nt)];
% ', dv = ', num2str(dv_obs), ' [m/s]', ', \mu_{re} = ', num2str(mu), ' [m/s]', ', \sigma_{re} = ', num2str(sigma), ' [m/s]', ...
%    ', \mu_{True} = ', num2str(mu_true), ' [m/s]', ', \sigma_{True} = ', num2str(sigma_true), ' [m/s]'];

xl = '\sigma_{True} [m/s]';

f = figure(1); hold on; f.Position = [10 10 1000 1000];
color = colors(4).c;
yl =  ['\sigma_{re}'];

marker = markers(1);
plott2(Sigmavec, sigma, xl, yl, txt, 2, dtext, color, marker);

marker = markers(2);

plott2(3.2./(L .* Nt) .* ones(1, length(Sigmavec)), ...
    sigma, xl, yl, txt, 2, dtext, color, marker);

% marker = markers(3);
% 
% plott2(Sigmavec, ...
%     Sigmavec/sqrt(Nt), xl, yl, txt, 2, dtext, color, marker);
