function [geo] = Geo(input_info)

r0 = input_info.r0;
phi0 = input_info.phi0;

dr = input_info.dr;
dph = input_info.dph;

if input_info.spatial_dist.type == 1
    
    drho = makedist('poisson', 'lambda', input_info.spatial_dist.lamr);
%     tdrho = drho.truncate(r0-dr/2, r0+dr/2);
    r = drho.random(1, input_info.NScatters);
    r = r0 - dr/2 + dr .* r/max(r);
  
    dphi = makedist('poisson', 'lambda', input_info.spatial_dist.lamp);
%     tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
    phi = dphi.random(1, input_info.NScatters);
    phi = phi0 - dph/2 + dph .* phi/max(phi);
    
elseif input_info.spatial_dist.type == 2
    
    drho = makedist('Gaussian', 'mu', input_info.spatial_dist.mu, 'sigma', inputinfo.spatial_dist.sigma);
    tdrho = drho.truncate(r0-dr/2, r0+dr/2);
    r = tdrho.random(1, input_info.NScatters);
    
    dphi = makedist('Gaussian', 'mu', input_info.spatial_dist.mu, 'sigma', inputinfo.spatial_dist.sigma);
    tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
    phi = tdphi.random(1, input_info.NScatters);
    
else
%     drho = makedist('poisson', 'lambda', input_info.spatial_dist.lamr);
% %     tdrho = drho.truncate(r0-dr/2, r0+dr/2);
%     r = drho.random(1, input_info.NScatters);
%     r = r0 - dr/2 + dr .* r/max(r);
%   
%     dphi = makedist('poisson', 'lambda', input_info.spatial_dist.lamp);
% %     tdphi = dphi.truncate(phi0-dph/2, phi0+dph/2);
%     phi = dphi.random(1, input_info.NScatters);
%     phi = phi0 - dph/2 + dph .* phi/max(phi);
    
    r = r0 - dr/2 + dr .* rand([1 input_info.NScatters]);
    phi = phi0 - dph/2 + dph .* rand([1 input_info.NScatters]);
end

xc0 = r .* cos(phi);
yc0 = r .* sin(phi);
rc0 = sqrt(xc0.^2 + yc0.^2); 


%% Plot the scatterers

if input_info.plot_geo == 1
    figure; scatter(xc0, yc0, 50, rc0); 
    colormap('copper'); colorbar; 
    xlabel('x [m]'); 
    ylabel('y [m]'); 
    zlabel('z [m]'); 
    
    
   if input_info.spatial_dist.type == 1
      title('Distributed according to Poisson distribution');
   elseif input_info.spatial_dist.type == 2
       title('Distributed according to Gaussian distribution');
   else
       title('Distributed according to Uniform distribution');
   end
   
   
end

geo.r = r;
geo.phi = phi;
geo.xc0 = xc0;
geo.yc0 = yc0;
geo.rc0 = rc0;

end