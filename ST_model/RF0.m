function [R1] = RF0(x_, y_, fx, fy)

R0 = rand(length(x_), (length(y_))) + 1j .* rand(length(x_), (length(y_))); 


% figure; surface(x_, y_, abs(R0).'); shading flat;
F = fft2(R0);

nx1 = 1 + floor(0.9999.*fx*length(x_));
nx2 = length(x_) - nx1 + 1; 

ny1 = 1 + floor(0.9999.*fy*length(y_));
ny2 = length(y_) - ny1 + 1;


F(nx1:nx2, :) = 0;
F(:, ny1:ny2) = 0;

R1 = real(ifft2(F));

end