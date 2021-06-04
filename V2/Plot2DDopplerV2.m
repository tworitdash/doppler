function [] = Plot2DDopplerV2(N, mean_Phi, v_amb, Signal, OI, SNR_db, Omega_rpm) 
   
    Omega_interest = Omega_rpm(OI);
    
    s = Signal(OI).doppler;
    vel_axis = linspace(-v_amb, v_amb, N);
   
    figure; imagesc(vel_axis, mean_Phi * 180/pi, db(abs(s))); shading interp; colorbar; colormap('jet'); 
    xlim([vel_axis(1) vel_axis(end)]);
    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 14);
    xlabel('Velocity axis [m/s]', 'FontSize', 14)
    zlabel('Doppler spectrum [dB]', 'FontSize', 14);
    title(['Doppler spectrum [dB] at ', num2str(Omega_interest), 'RPM', ', SNR = ', num2str(SNR_db), ' dB'], 'FontSize', 12);

end