clear; 
close all;

L = 20;
N = 10240;

dx = 0.01;
dt = 0.001;

x_ = 0:dx:L;

t_ = 0:dt:N*dt;

[x, t] = meshgrid(x_, t_);

eta = sin(2 .* acot(exp((2 .* pi * t + L .* log(cot(pi .* x/L + eps)))./L)));

txt = ['Intial field defined at time =  ', num2str(t(1)), ' [s]'];
figure; plot(x_, (eta(1, :)), 'o', 'DisplayName', txt); 


txt = ['Analytical at time =  ', num2str(t(end)), ' [s]'];
hold on ; plot(x_, (eta(end, :)), '*-', 'DisplayName', txt);grid on;
 
eta_n = zeros(size(x));

eta_n(1, :) = sin(2 .* pi .* (x_)/L);

for k = 1:length(t_)
    for i = 1:length(x_)
        if i == 1
            eta_n(k + 1, i) = eta_n(k, i) - sin(2 * pi * x_(i)/L)./(2 * dx) .* dt .* (eta_n(k, i+1) - eta_n(k, i));
        elseif i == length(x_)
            eta_n(k + 1, i) = eta_n(k, i) - sin(2 * pi * x_(i)/L)./(2 * dx) .* dt .* (eta_n(k, i) - eta_n(k, i - 1));
        else
            eta_n(k + 1, i) = eta_n(k, i) - sin(2 * pi * x_(i)/L)./(2 * dx) .* dt .* (eta_n(k, i+1) - eta_n(k, i - 1));
        end
    end 
end
txt = ['Numerical at time =  ', num2str(t(end)), ' [s]'];

hold on; plot(x_, (eta_n(end, :)), '*-', 'DisplayName', txt);grid on;

legend;

% plot(x_, (eta_n(1, :)), 'o'); hold on ; plot(x_, (eta_n(end, :)), '*-');grid on;