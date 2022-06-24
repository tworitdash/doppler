
%% Weather simulator

clear; 
close all;

input_info.r0 = 100;
input_info.phi0 = 0*pi/180;
input_info.dr = 0;
input_info.dph = 0; % 0.001*pi/180;

input_info.spatial_dist.type = 1;
input_info.plot_geo = 1;
input_info.NScatters = 1;

input_info.spatial_dist.lamr = 10;
input_info.spatial_dist.lamp = 20;

input_info.RADAR.dT = 1e-3;

input_info.RADAR.dTjittervec = 0; % linspace(0, 4*input_info.RADAR.dT, 5);

input_info.N_pulse = 5;

input_info.Ngap_avg = 20; 

input_info.sig_gap = 2 * input_info.N_pulse * 2; % linspace(0, 2, 3);

% Temp = linspace(1, 10, 10);

% for k = 1:length( input_info.sig_gap )

input_info.RADAR.dTjitter = 0; % input_info.RADAR.dTjittervec(k); % input_info.RADAR.dT/2;
input_info.RADAR.lambda = 3e-2;


input_info.N_rot         = 2;

input_info.N_gap         = [0 randi([input_info.Ngap_avg-input_info.sig_gap/2 input_info.Ngap_avg+input_info.sig_gap/2], 1, input_info.N_rot)];
input_info.velocity.u.mu = 3;
input_info.velocity.v.mu = 0;
input_info.velocity.u.sigma = 0;
input_info.velocity.v.sigma = 0;

input_info.SNR = 20;
input_info.velocity.type = 2;
input_info.Doppler_plot = 0;
% input_info.Ngt = 128;
input_info.Ngt = [];
[out] = SimRad(input_info);

%% Likelihood check

llinfo = input_info;
llinfo.v_amb = input_info.RADAR.lambda/(2*out.dTavg);
llinfo.Nutest = 10000;

llinfo.velocity.u.muvec = linspace(eps, llinfo.v_amb, llinfo.Nutest);



llinfo.plot_geo = 0;

[geo] = Geo(llinfo);

llinfo.xc0 = geo.xc0;
llinfo.yc0 = geo.yc0;
llinfo.rc0 = geo.rc0;


llinfo.t_avail = out.t_avail;

% [u_axis, t_axis] = meshgrid(llinfo.velocity.u.muvec, llinfo.t_avail);

%% Deterministic Annealing

% Temp = linspace(1e4, 0.1, 10);
Temp = [1e4 1];

muvec = zeros(llinfo.Nutest, length(Temp));
muvec(:, 1) = llinfo.velocity.u.muvec;

for l = 2:length(Temp)
    
%     disp(l);
    
    for i = 1:llinfo.Nutest

        llinfo.velocity.u.mu = muvec(i, l - 1); % Overwrite the mean velocity of u to compute the likelihood

        [LLout(i)] = LL(llinfo, out);
        
        Pr(i, l) = (eps + exp(LLout(i)/Temp(l)));

%         E(i, l) = (-LLout(i)) * Pr(i, l);

        Num1(i, l) = ( imag(out.Z_avail_vec) * out.t_avail.' .* Pr(i, l) );
        Num2(i, l) = ( real(out.Z_avail_vec) * out.t_avail.' .* Pr(i, l) );
       
    end
    
    figure(101); hold on; plot(llinfo.velocity.u.muvec, db(Pr(:, l))); 
    
%     Norm(l) = sum(Pr(:, l));
    
%     muvec(:, l) = atan2( Num1(:, l), Num2(:, l) ) .* input_info.RADAR.lambda / (4 * pi);
    muvec(:, l) = muvec(:, l -1);
    
%     figure(102); hold on; histogram(muvec(:, l));
    
end



