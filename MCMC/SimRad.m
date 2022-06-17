%% Radar data simulator for one specific resolution cell

% INPUT:

%       1. Number of Scatters                                               N_{scatters}
%       2. Mean Range                                                       r_0
%       3. Range resolution                                                 dr
%       4. Azimuth resolution                                               dphi
%       3. Mean Azimuth                                                     phi_0
%       4. Standard Deviaton ratio                                          p
%       5. Spatial Distribution Dist_space                                 ´Gaussian´, ´Poisson´, etc...
%       6. Gap samples                                                      N_gap
%       7. Number of rotations                                              N_rot
%       8. Number of pulses per sector                                      N_pulses
%      10. plot_Geometry [0 or 1]                                           Plot_Geo
%      11. plot_Spectrum [0 or 1]                                           Plot_spectrum
%      12. SNR (dB)


%% OUTPUT:
%       1. signal available
%       2. signal in HD

% function [output_info.sig, output_info.sig_HD] = SimRad(input_info)
function [] = SimRad(input_info)


r0 = input_info.r0;
phi0 = input_info.r0;

dr = input_info.dr;
dph = input_info.dph;

%% Geometry of the scatterers

if input_info.spatial_dist.type == 1
    
    drho = makedist('poisson', 'lambda', inputinfo.spatial_dist.lamr);
    tdrho = drho.truncate(r0-dr/2, r0+dr/2);
    r = tdrho.random(1, input_info.NScatters);
    
    dphi = makedist('poisson', 'lambda', inputinfo.spatial_dist.lamp);
    tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
    phi = tdphi.random(1, input_info.NScatters);
    
elseif input_info.spatial_dist.type == 2
    
    drho = makedist('Gaussian', 'mu', inputinfo.spatial_dist.mu, 'sigma', inputinfo.spatial_dist.sigma);
    tdrho = drho.truncate(r0-dr/2, r0+dr/2);
    r = tdrho.random(1, input_info.NScatters);
    
    dphi = makedist('Gaussian', 'mu', inputinfo.spatial_dist.mu, 'sigma', inputinfo.spatial_dist.sigma);
    tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
    phi = tdphi.random(1, input_info.NScatters);
    
else
    drho = makedist('poisson', 'lambda', inputinfo.spatial_dist.lamr);
    tdrho = drho.truncate(r0-dr/2, r0+dr/2);
    r = tdrho.random(1, input_info.NScatters);
    
    dphi = makedist('poisson', 'lambda', inputinfo.spatial_dist.lamp);
    tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
    phi = tdphi.random(1, input_info.NScatters);
end

x0 = r .* cos(phi);
y0 = r .* sin(phi);
r0 = sqrt(x0.^2 + y0.^2); 


%% Plot the scatterers

if input_info.plot_geo == 1
    figure; scatter(x0, y0, 50, r0); colormap('copper'); colorbar; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 
end
