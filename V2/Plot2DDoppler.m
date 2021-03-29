function [] = Plot2DDoppler(vel_axis, mean_Phi, Signal, SI, OI) 
    SNR_interest = SNR(SI);
    Omega_interest = Omega_rpm(OI);
    
    figure; surface(vel_axis, mean_Phi * 180/pi, abs(squeeze(Signal(OI).doppler(SI, :, :)))); shading flat; colorbar; %colormap('jet'); 
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Velocity axis [m/s]', 'FontSize', 16)
    zlabel('Normalized Doppler spectrum', 'FontSize', 16);
    title(['Normalized Doppler spectrum at ', num2str(Omega_interest), 'RPM', ', SNR = ', num2str(SNR_interest), ' dB'], 'FontSize', 16);

end