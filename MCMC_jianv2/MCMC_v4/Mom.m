function [PT, mu, sigma] = Mom(FFT, u)

dv = u(2) - u(1);

PT = sum(abs(FFT).^2  .* dv);
mu = 1./PT .* sum(u .* abs(FFT).^2 .* dv);
sigma = sqrt(sum(1./PT .* (u - mu).^2 .* abs(FFT).^2 .* dv));

% pd = fit(u.', abs(FFT).', 'gauss1');

end