function [c, w] = change_sigma(u_axis, cells, epsilon, u_obs, sigma_obs, w, c)
    var = abs(u_axis(end-cells:end) - mean(u_axis(end-cells:end)))./(norm(u_axis(end-cells:end)));
    corrs = xcorr(var, var); 
    corrs = corrs(2:end-1);
    corrs = corrs(floor(length(corrs)/2):end);
    corrs = mean(corrs);
   
    
%     if (corrs < epsilon) && (  (mean(u_axis(end-cells:end)) > u_obs + sigma_obs ) &&  (mean(u_axis(end-cells:end)) < u_obs - sigma_obs ) )
%         c = 1;
%         u_mean = mean(u_axis(end-cells:end));
%         w = [w sqrt(abs(u_mean - u_obs).^2./abs(u_mean - muu_obs).^2)];
%     end
    if (corrs < epsilon)
        u_mean = mean(u_axis(end-cells:end));
        w = [w sqrt(abs(u_mean - u_obs).^2./abs(u_obs).^2)];
        if w(end) < 0.01
            c = [c 1];
        else
            c = [c 0];
        end
    end
        
    
end