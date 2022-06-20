function [E0] = Eta_0(x_, y_, x, y, dx, dy, lambdax, lambday, fx, fy)

    R0 = (rand(size(x)) + 1j .* rand(size(x)) - 0.5);
    kx0 = 2*pi/lambdax;
    ky0 = 2*pi/lambday;

    Kx_ = linspace(-kx0, kx0, 100);
    dKx = Kx_(2) - Kx_(1);

    Ky_ = linspace(-ky0, ky0, 100);
    dKy = Ky_(2) - Ky_(1);

    [Kx, Ky] = meshgrid(Kx_, Ky_);


    for i = 1:length(Kx_)
        for l = 1:length(Ky_)

            Spec(l, i) = sum(sum(R0 .* exp(-1j .* (Kx_(i) .* x + Ky_(l) .* y)) .* dKx .* dKy));

        end
    end


    % figure; surface(Kx, Ky, abs(Spec)); shading flat; colormap('jet'); colorbar;
    Spec_modified = Spec;

    nx1 = 1 + floor(0.9999.*fx*length(Kx_));
%     nx2 = length(Kx_) - nx1 + 1; 

    ny1 = 1 + floor(0.9999.*fy*length(Ky_));
%     ny2 = length(Ky_) - ny1 + 1;


    Spec_modified(nx1:end, :) = 0;
    Spec_modified(:, ny1:end) = 0;



    % figure; surface(Kx,Ky, abs(Spec_modified)); shading flat; colormap('jet'); colorbar;

    for i = 1:length(x_)
        for l = 1:length(y_)
            E0(l, i) = sum(sum(Spec_modified .* exp(1j .* (Kx .* x_(i) + Ky .* y_(l))) .* dx .* dy));
        end
    end


%     figure; surface(x, y, (abs(E0))); shading flat; colormap('jet'); colorbar;

end