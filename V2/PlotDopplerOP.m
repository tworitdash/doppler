function [] = PlotDopplerOP(OP_enable, OP_enable_error, BW_deg, BI, SNR_db, SI, Phi, Omega_rpm, v_mean, v_spread, v_mean_e, v_spread_e)

SNR_interest = SNR_db(SI);
BW_interest = BW_deg(BI);
Length_Phi_axis = length(Phi(BI).Phi) - 1;

if OP_enable == 1
    figure; surface(Omega_rpm, Phi(BI).Phi(1:Length_Phi_axis) .* 180/pi, squeeze(v_mean(BI, SI, :, 1:Length_Phi_axis)).'); shading flat; colorbar;colormap('jet') % Plot the error in mean Doppler

    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
    zlabel('Mean Doppler velocity [m/s]', 'FontSize', 16);
    title(['Mean Doppler velocity [m/s] at ', ', SNR = ', num2str(SNR_interest), ' dB', ', BW = ', num2str(BW_interest), ' deg'], 'FontSize', 12);
    
    figure; surface(Omega_rpm, Phi(BI).Phi(1:Length_Phi_axis) .* 180/pi, squeeze(v_spread(BI, SI, :, 1:Length_Phi_axis)).'); shading flat; colorbar; colormap('jet'); % Plot the mean Doppler

    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
    zlabel('Doppler Spectrum Width [m/s]', 'FontSize', 16);
    title(['Doppler spectrum width [m/s] at ', ', SNR = ', num2str(SNR_interest), ' dB' ', BW = ', num2str(BW_interest), ' deg'], 'FontSize', 12);
end
if OP_enable_error == 1
    
    figure; imagesc(Omega_rpm, Phi(BI).Phi(1:Length_Phi_axis) .* 180/pi, squeeze(v_mean_e(BI, SI, :,  1:Length_Phi_axis)).'); shading flat; colorbar;colormap('jet') % Plot the error in mean Doppler

    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
    zlabel('Error in Mean Doppler velocity [m/s]', 'FontSize', 16);
    title(['Error in Mean Doppler velocity [m/s] at ', ', SNR = ', num2str(SNR_interest), ' dB' ', BW = ', num2str(BW_interest), ' deg'], 'FontSize', 12);
    
    figure; imagesc(Omega_rpm, Phi(BI).Phi(1:Length_Phi_axis) .* 180/pi, squeeze(v_spread_e(BI, SI, :, 1:Length_Phi_axis)).'); shading flat; colorbar; colormap('jet'); % Plot the mean Doppler

    ylabel('Azimuth Angle \Phi [deg]', 'FontSize', 16);
    xlabel('Angular speed of radar beam in azimuth [rpm]', 'FontSize', 16)
    zlabel('Error in Doppler Spectrum width [m/s]', 'FontSize', 16);
    title(['Error in Doppler spectrum width[m/s] at ', ', SNR = ', num2str(SNR_interest), ' dB'', BW = ', num2str(BW_interest), ' deg' ], 'FontSize', 12);
end

end