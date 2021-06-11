dth = pi/180;  % 1 deg
dph = pi/180;

th = -pi/2:dth:pi/2;
ph = eps:dph:2*pi;

[th, ph] = meshgrid(th, ph);

U = cos(th).^2; % .* cos(ph).^2;


P_rad_i = U(:, 91:end) .* sin(th(:, 91:end)) .* dth .* dph; 

P_rad = sum(sum(P_rad_i));

D_ = U/P_rad .* (4*pi);



figure; plot(th(1, :).*180/pi, db(D_(1, :))./2); ylim([-50 50]); grid on;