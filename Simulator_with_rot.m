function [sig, sig_doppler, sig_with_Omega, hits_scan, delta_v, vel_axis, time_axis] = Simulator_with_rot(Omega_rpm, BW, SNR_db, mean_Phi, beta_wind, PRT, lambda, mu, sigma)

SNR = 10^(SNR_db./20); % SNR in linear scale


N_burst = 256; % Only used when the rotation speed is 0 [RPM]

%% Signal generator and Doppler processing
 
       Omega = Omega_rpm .* 2.*pi./60; % Angular velocity of the radar beam in rad/sec
     
                if Omega ~= 0
                    T = BW/Omega;
                else
                    T = N_burst .* PRT;
                end

                time_axis = eps:PRT:T; % Time axis in terms of multiples of PRT for one resolution step in angle
                hits_scan = length(time_axis); % length of time axis

                delta_v = lambda/(2*hits_scan*PRT);% velocity resolution
                v_amb = lambda/(4*PRT);% Doppler ambiguity limits in velocity
                vel_axis = linspace(-v_amb,v_amb,hits_scan);% velocity axis

                Nifft = hits_scan; % Number of DTFT points for ifft

                beta = beta_wind - mean_Phi; % Angle exerted by the radar beam with the wind direction at each time step

             
                X = rand(1, hits_scan);  % Random number generator for the simulator amplitude for all the velocities
                Theta = rand(1, hits_scan) .* 2 .* pi; % Random phase generator for all the velocities (uncorelated with X in previous line)
                

             

                [sig, sig_doppler] = Doppler_spectrum_TD(vel_axis, mu, sigma, Nifft, SNR, lambda, X, Theta); 

               
                sig_with_Omega = abs((sig)) .* exp(1j * angle(sig).* cos(beta));% + noise(i).n;

        

end
