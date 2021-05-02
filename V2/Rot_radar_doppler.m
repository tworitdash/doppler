clear;
% close all;

%% Configuration 0 is off, 1 in on

% Want plot for Error vs SNR? 
SNR_enable = 0;

% Want to plot for Error vs BW?
BW_enable = 0;

% Want plot for Erorr vs Omega and Phi?
OP_enable_error = 0;

% Want Doppler Azimuth plot?

DS_Azimuth_plots = 1; % Doppler Azimuth Plot

% Want to run Monte Carlo simulation?
MC_enable = 1;

% Want to plot mean Doppler and Doppler spectrum width?
OP_enable = 0;

%%

if MC_enable == 1
    n_MC = 32; 
else
    n_MC = 1;
end


if SNR_enable == 1

    SNR_i_db = -30; 
    SNR_f_db = 30; 
    n_SNR = 10; % Number of points wanted in the SNR axis
    SNR_db = linspace(SNR_i_db, SNR_f_db, n_SNR);
    SNR = 10.^(SNR_db/10);

else
    n_SNR = 1;
    SNR_db = 30;
    SNR = 10^(SNR_db/10);
end

if BW_enable == 1
    n_BW = 10;
    BW_deg = linspace(1, 10, n_BW);
else
    BW_deg = 1;
    n_BW = 1;
end


%% RADAR constants and wind direction

beta_wind = eps; % wind direction
mu = 5;
sigma = 0.2;



PRT = 1e-3;
lambda = 3e-2;
n = 2^10;

v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity

% for l = 1:n_MC
%     for s = 1:n_SNR
%         [data_(l, s, :), data_f_(l, s, :)] = DS_simulator(SNR(s), mu, sigma, n, v_amb); % One time simulation of the spectrum and the signal with fine sampling
%     end
% end

% data = squeeze(mean(data_, 1));
% data_f = squeeze(mean(data_f_, 1));

% n_ = 4;
% vel_ = linspace(-v_amb, v_amb, n);
% vel__ = linspace(-v_amb, v_amb, n_);
% 
% figure; plot(vel_, abs(data_f)); hold on; plot(vel__, 1/sqrt(n) .* abs(fftshift(fft(data, n_))));

%% 

% Omega_rpm = linspace(1, 60, 4); % RPM axis for the rotation of the radar
% Omega_rpm = 160;
Omega_rpm = linspace(1, 60, 4);
% Omega_rpm = 40.33;

%% 

for l = 1:n_MC % Monte Carlo iterations v
    

for m = 1:n_BW

Phi_0_deg = 0;
Phi_end_deg = 360;

BW = BW_deg(m)*pi/180; % Beam width in radians
phi_0 = Phi_0_deg * pi/180; % Start of azimuth angle
phi_end = Phi_end_deg * pi/180; % End of azimuth angle

Phi(m).Phi = phi_0:BW:phi_end;

% mean_Phi = mean([Phi(1:end-1); Phi(2:end)]); % This is done to take the mid angles of all possible angular resolution cell

%% resampling based on the rotation of the radar

    for s = 1:n_SNR
                            
%         [data(l, s, :), data_f(l, s, :)] = DS_simulator(SNR(s), 1e3, mu, sigma, n, v_amb);
        
%         [data, data_f] = DS_simulatorV2(SNR, m0, mu, sigma, n, v_amb, N)
        
        for i = 1:length(Omega_rpm)

            Omega = Omega_rpm(i) .* 2 * pi ./ 60;

            T(m ,i) = BW/Omega;

            time_axis(m, i).axis = eps:PRT:T(m, i);
            

            hits_scan(m, i) = length(time_axis(m, i).axis); % length of time axis
%             if hits_scan(i) == 4
%                 hits_scan(i) = 5;
%             end

           

%             hits_scan_(m, i) = 2^(nextpow2(hits_scan(m, i)) - 1) + 1; % hits scan for Doppler processing
             hits_scan_(m, i) = 2^(nextpow2(hits_scan(m, i)) - 1);
             
             if hits_scan_(m, i) < 32
                hits_scan_1(m, i) = 64;
             else
                 hits_scan_1(m, i) = hits_scan_(m, i); % +1 is used for odd number of samples
             end
%             time_axis_full = linspace(eps, PRT.*n, n); % High resolution time axis
%             
%             
            time_axis_(m, i).axis = linspace(eps, PRT.*n, hits_scan_(m, i));
%           Axis(m, i).axis =((0:hits_scan_1(m, i)-1) -ceil((hits_scan_1(m, i)-1)/2))/hits_scan_1(m, i);
            Axis(m, i).axis = linspace(-hits_scan_1(m, i)/2, hits_scan_1(m, i)/2 - 1, hits_scan_1(m, i))/hits_scan_1(m, i);
          
            vel_axis(m, i).axis =2*v_amb*Axis(m, i).axis;
            
%             vel_axis(m, i).axis = linspace(-v_amb, v_amb, hits_scan_1(m, i)); % To show the Doppler spectrum
            
            
            vel_axis_(m, i).axis = linspace(-v_amb, v_amb, hits_scan_(m, i)); % For mean and spectrum width calculation 
            vel_axis_full = linspace(-v_amb, v_amb, n);
            
            
            delta_v(m, i) = lambda/(2*hits_scan_(m, i)*PRT);
            
            

            for k = 1:length(Phi(m).Phi)-1
                
                %disp(m);
                beta_scan = beta_wind - linspace(Phi(m).Phi(k), Phi(m).Phi(k + 1), hits_scan_(m, i));
                beta_scan_hd = beta_wind - linspace(Phi(m).Phi(k), Phi(m).Phi(k + 1), n);
                [data(i, m).data(l, s, k, :), data(i, m).data_f(l, s, k, :)] = DS_simulatorV2(SNR(s), 1, mu, sigma, n, v_amb, hits_scan_(m, i));
             
%                [Signal(i, m).sig(l, s, k, :)] = DS_simulatorV3_with_az(SNR(s), 1, mu, sigma, n, v_amb, hits_scan_(m, i), beta_scan_hd, 0);
               
                
             Signal(i, m).sig(l, s, k, :) = ((abs(squeeze(data(i, m).data(l, s, k, :)))...
                 .* exp(1j .* unwrap(angle(squeeze(data(i, m).data(l, s, k, :)))) .* cos(beta_scan).')));
             
             
%                 Signal(i, m).doppler_whole(l, s, k, :) = 1./sqrt(n) .* fftshift(fft(Signal(i, m).sig(l, s, k, :), n));
                
%                 idx = 1:n;
%                 idxq = linspace(min(idx), max(idx), hits_scan_(m, i));
% 
%                 Signal(i, m).sig_resampled(l, s, k, :) = interp1(idx, squeeze(Signal(i, m).sig(l, s, k, :)), idxq, 'spline'); % resample the HR signal


%                 
%                 figure; plot(time_axis_(m, i).axis, real(squeeze(Signal(i, m).sig_resampled(l, s, k, :))), '*');
%                 hold on; plot(time_axis_full, real(squeeze(Signal(i, m).sig(l, s, k, :))), 'o');
                
%                 Signal(i, m).doppler_resampled(l, s, k, :) = 1./sqrt(hits_scan_(m, i)) .* fftshift(fft(Signal(i, m).sig_resampled(l, s, k, :), hits_scan_(m, i)));
                
                Signal(i, m).doppler(l, s, k, :) = 1./sqrt(hits_scan_1(m, i)) .* fftshift(fft(Signal(i, m).sig(l, s, k, :), hits_scan_1(m, i)));
                
%                 figure; plot(vel_axis(m, i).axis, db(abs(squeeze(Signal(i, m).doppler(l, s, k, :)).^2))/2, '-o'); 
%                 
%                 hold on;  plot(vel_axis(m, i).axis, db(abs(squeeze(data(i, m).data_f(l, s, k, :))).^2)/2, '*'); grid on;
                
%                 figure; 
%                 plot(vel_axis_full, abs(squeeze(Signal(i, m).doppler_whole(l, s, k, :))));
%                 hold on; 
%                 plot(vel_axis_(m, i).axis, abs(squeeze(Signal(i, m).doppler_resampled(l, s, k, :))), 'x');
%                 hold on; 
%                 plot(vel_axis_(m, i).axis, abs(squeeze(Signal(i, m).doppler(l, s, k, :))), 'o', 'color', 'k');
                
                %% Signal Doppler for calculation of the mean Doppler velocity and Doppler spectrum width
                
%                 Signal(i, m).doppler_(l, s, k, :) = 1./sqrt(hits_scan_(m, i)) .* (fft(Signal(i, m).sig(l, s, k, :), hits_scan_(m, i)));
                
                
%                 Signal(i, m).dop_for_process(l, s, k, :) = 1./sqrt(hits_scan_(m, i)) .* (fft(Signal(i, m).sig(l, s, k, :), hits_scan_(m, i)));
                
                
%                 if cos(beta_scan(end)) < 0
%                     vel_axis_proc(m, i).axis = (vel_axis(m, i).axis(1:hits_scan_(m, i)/2+1));
%                     
%                     Signal_dop_required = squeeze(Signal(i, m).dop_for_process(l, s, k, hits_scan_(m, i)/2:hits_scan_(m, i)));
%                 else
%                     vel_axis_proc(m, i).axis = vel_axis(m, i).axis(hits_scan_(m, i)/2:hits_scan_(m, i));
%                     Signal_dop_required = squeeze(Signal(i, m).dop_for_process(l, s, k, 1:hits_scan_(m, i)/2+1));
%                 end
%                 
% %                 plot(vel_axis_proc(m, i).axis, abs(squeeze(Signal_dop_required)), 'o', 'color', 'k');
%                 
%                 PT_integrand = abs(Signal_dop_required).^2 .* delta_v(m, i);
%                 PT = sum(PT_integrand); % Total power of the Doppler Spectrum
% 
%                 v_mean_integrand = vel_axis_proc(m, i).axis.' .* abs(Signal_dop_required).^2 .* delta_v(m, i);
%                 v_mean_l(l, m, s, i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 
% 
%                 v_spread_integrand = (vel_axis_proc(m, i).axis.' - v_mean_l(l, m, s, i, k)).^2 .* abs(Signal_dop_required).^2 .* delta_v(m, i);
%                 v_spread_l(l, m, s, i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
                
                
                PT_integrand = abs(squeeze(Signal(i, m).doppler(l, s, k, :))).^2 .* delta_v(m, i);
                PT = sum(PT_integrand); % Total power of the Doppler Spectrum

                v_mean_integrand = vel_axis(m, i).axis.' .* abs(squeeze(Signal(i, m).doppler(l, s, k, :))).^2 .* delta_v(m, i);
                v_mean_l(l, m, s, i, k) = sum(v_mean_integrand) ./ PT; % Mean Doppler velocity 

                v_spread_integrand = (vel_axis(m, i).axis.' - v_mean_l(l, m, s, i, k)).^2 .* abs(squeeze(Signal(i, m).doppler(l, s, k, :))).^2 .* delta_v(m, i);
                v_spread_l(l, m, s, i, k) = sqrt(sum(v_spread_integrand)./ PT); % Doppler spectrum width
            end
            
    %         figure; surface(vel_axis(i).axis, mean_Phi * 180/pi, abs(squeeze(Signal(i).doppler(s, :, :)))); shading flat; colorbar; %colormap('jet'); 
    %         ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    %         xlabel('Velocity axis [m/s]', 'FontSize', 16)
    %         zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    %         title(['Normalized Doppler spectrum at ', num2str(Omega_rpm(i)), 'RPM', ', SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 16);
        end

    end
end

end

% v_mean = squeeze(mean(v_mean_l, 1));
% v_spread = squeeze(mean(v_spread_l, 1));

%% Error calculations

for m = 1:n_BW
    for s = 1:n_SNR
       for i = 1:length(Omega_rpm)
          for k = 1:length(Phi(m).Phi)-1
               v_mean(m, s, i, k) = sum(v_mean_l(:, m, s, i, k))/n_MC;
               v_spread(m, s, i, k) = sum(v_spread_l(:, m, s, i, k))/n_MC;
               v_mean_e(m, s, i, k) = sqrt(sum((v_mean_l(:, m, s, i, k) - mu).^2)/n_MC);
               v_spread_e(m, s, i, k) = sqrt(sum((v_spread_l(:, m, s, i, k) - sigma).^2)/n_MC);
          end
       end
    end
end


%% Plot Doppler spectrum (Azimuth vs Doppler Velocity)


if DS_Azimuth_plots == 1
    SI = 1; % Index for SNR
    OI = 2; % Index for Omega
    BI = 1; % Index of Beamwidth
    Plot2DDoppler(vel_axis(BI, OI).axis, Phi, Signal, BI, SI, OI, BW_deg, SNR_db, Omega_rpm);
end
%%
% for k = 1:length(Phi(BI).Phi)-1
%     s = squeeze(mean(abs(squeeze(Signal(OI, BI).doppler(:, SI, k, :))), 1));
%     figure(101); plot(vel_axis(BI, OI).axis, db(abs(s)));hold on;
% end
%% Plot Erros with BW

if BW_enable == 1
    OI = 1;
    PI = 1; 
    SI = 1;
    PlotDopplerBW(BW_deg, v_mean_e, v_spread_e, SI, OI, PI, SNR_db, Omega_rpm, Phi);
end


%% Error with With SNR

if SNR_enable == 1
    OI = 4; % Index for Omega
    PI = 1; % Index for Phi
    BI = 1;
    PlotDopplerSNR(SNR_db, v_mean_e, v_spread_e, BI, OI, PI, BW_deg, Omega_rpm, Phi);
end
% %%
% legend show;

%% 2D plots of mean Doppler and Doppler spread and Erros

if OP_enable == 1 || OP_enable_error == 1
    SI = 10; % Index of the SNR axis
    BI = 1; % Index of the Beamwidth axis
    PlotDopplerOP(OP_enable, OP_enable_error, BW_deg, BI, SNR_db, SI, Phi, Omega_rpm, v_mean, v_spread, v_mean_e, v_spread_e);
end

%% 1D plots of velocity with azimuth at different rotation speeds but a given Beamwidth


figure;
SI = 1; % Index of the SNR axis
BI = 1; % Index of BeamWidth
Length_Phi_axis = length(Phi(BI).Phi) - 1;


for i = 1:length(Omega_rpm)
   txt = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
   hold on; plot(Phi(BI).Phi(1:Length_Phi_axis) * 180/pi, squeeze(v_mean(BI, SI, i, 1:Length_Phi_axis)), 'DisplayName', txt);
end
h = legend;
grid on;
ylabel('Mean Doppler velocity [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Mean Doppler velocity when ', ' SNR = ', num2str(SNR_db(SI)), ' dB', ', BW = ', num2str(BW_deg(BI)), ' deg'], 'FontSize', 12);
figure;
for i = 1:length(Omega_rpm)
   txt = ['\Omega = ', num2str(Omega_rpm(i)), ' [rpm]'];
   hold on; plot(Phi(BI).Phi(1:Length_Phi_axis) * 180/pi, squeeze(v_spread(BI, SI, i, 1:Length_Phi_axis)), 'DisplayName', txt);
end
h = legend;
grid on;
ylabel('Doppler spectrum width [m/s]', 'FontSize', 12);
xlabel('Azimuth \Phi [deg]', 'FontSize', 12)
set(h,'FontSize',12, 'FontWeight', 'bold', 'Location','north');
title(['Doppler spectrum width when ', ' SNR = ', num2str(SNR_db(SI)), ' dB',', BW = ', num2str(BW_deg(BI)), ' deg'], 'FontSize', 12);

%% Cosine compensation basics


% figure;
% SI = 1; % Index of the SNR axis
% BI = 1; % Index of BeamWidth
% OI = 1; % Index of the rotation speed vector
% Length_Phi_axis = length(Phi(BI).Phi) - 1; % length of the azimuth axis based on the beamwidth of interest
% 
% Phi_axis = Phi(BI).Phi(1:Length_Phi_axis);
% 
% v_mean_s = squeeze(v_mean(BI, SI, OI, 1:Length_Phi_axis)); % mean Doppler for specific rotation speed, SNR and Beamwidth with respect to the azimuth angles
% 
% figure; plot(Phi_axis*180/pi, v_mean_s, 'LineWidth', 2); grid on;
% 
% v_mean_s_fft = abs(fftshift(fft(v_mean_s))).^2/Length_Phi_axis;
% 
% figure; plot(v_mean_s_fft, 'LineWidth', 2); grid on;
% 
%% To do list
%% Make cosine compensation and use for multiple rotations

