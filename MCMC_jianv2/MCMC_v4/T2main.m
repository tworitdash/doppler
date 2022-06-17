%% Weather simulator

clear; 
close all;

input_info.r0 = 100;
input_info.phi0 = 0*pi/180;
input_info.dr = 0;
input_info.dph = 0; % 0.001*pi/180;

input_info.spatial_dist.type = 1;
input_info.plot_geo = 1;
input_info.NScatters = 2;

input_info.spatial_dist.lamr = 10;
input_info.spatial_dist.lamp = 20;

input_info.RADAR.dT = 1e-3;

input_info.RADAR.dTjittervec = 0; %input_info.RADAR.dT/3; % linspace(0, 4*input_info.RADAR.dT, 5);
input_info.N_pulse = 5;

input_info.Ngap_avg = 0;

input_info.sig_gap = 0; %2 * input_info.N_pulse .* linspace(0, 2, 3);

for k = 1:length( input_info.sig_gap )
    
input_info.RADAR.dTjitter = 0; % input_info.RADAR.dTjittervec(k); % input_info.RADAR.dT/2;
input_info.RADAR.lambda = 3e-2;


input_info.N_rot         = 32;

input_info.N_gap         = [0 randi([input_info.Ngap_avg-input_info.sig_gap(k)/2 input_info.Ngap_avg+input_info.sig_gap(k)/2], 1, input_info.N_rot)];
input_info.velocity.u.mu = 3;
input_info.velocity.v.mu = 0;
input_info.velocity.u.sigma = 1;
input_info.velocity.v.sigma = 0;

input_info.SNR = 20;
input_info.velocity.type = 2;
input_info.Doppler_plot = 0;
% input_info.Ngt = 128;
input_info.Ngt = [];


[out] = SimRad(input_info); % Signal generator



%% Likelihood check



llinfo = input_info;
llinfo.v_amb = input_info.RADAR.lambda/(2*out.dTavg);
llinfo.v_ambscan = input_info.RADAR.lambda/(2*max(input_info.RADAR.dT*input_info.N_gap(2:end)));
llinfo.Nutest = 1e6;

llinfo.velocity.u.muvec = linspace(eps, llinfo.v_amb, llinfo.Nutest); %0:llinfo.v_ambscan/2000:llinfo.v_amb; %linspace(eps, llinfo.v_amb, llinfo.Nutest);
llinfo.plot_geo = 0;

llinfo.t_known = 1; % 1 if time instances are known: 0 if not known

if llinfo.t_known == 0

    for l = 1:input_info.N_rot
        t1(l, :) = (l - 1)*(input_info.Ngap_avg + input_info.N_pulse)+[0:input_info.N_pulse-1];
    end

    t2 = reshape(t1.', [size(t1, 1) * size(t1, 2) 1]);


    llinfo.t_avail = t2 .* input_info.RADAR.dT;
else
    llinfo.t_avail = out.t_avail;
end

llinfo.scatters = 1;
[geo] = Geo(llinfo);

llinfo.xc0 = geo.xc0;
llinfo.yc0 = geo.yc0;
llinfo.rc0 = geo.rc0;

for i = 1:length(llinfo.velocity.u.muvec)
    
    llinfo.velocity.u.mu = llinfo.velocity.u.muvec(i); % Overwrite the mean velocity of u to compute the likelihood

    [LLout2(i)] = LL(llinfo, out);

end



%% plot likelihood

% dtext = ['Sample gap: ', num2str(input_info.N_gap), ', Number of rotations: ', num2str(input_info.N_rot), ...
%     ', Jitter: ', num2str(input_info.RADAR.dTjitter/input_info.RADAR.dT),' x dT ' ];
dtext = ['N_{sec, std}: ', num2str(input_info.sig_gap(k)/input_info.N_pulse), ', Number of rotations: ', num2str(input_info.N_rot), ...
    ', Jitter: ', num2str(input_info.RADAR.dTjitter/input_info.RADAR.dT),' x dT ', ' Time instances: unknown ' ];
% figure(101); hold on; plot(llinfo.velocity.u.muvec, db(-(LLout).*sinc(llinfo.velocity.u.muvec - input_info.velocity.u.mu)), 'DisplayName', dtext); grid on;

figure(101); hold on; plot(llinfo.velocity.u.muvec, (LLout2), 'DisplayName', dtext); grid on;

% figure(101); hold on; plot(llinfo.velocity.u.muvec, ((-LLout1+LLout2)), 'DisplayName', dtext); grid on;


end

xl = 'Velocity [m/s]';
yl = 'Log Likelihood [dB]';
% tl = 'Likelihood';
% title(tl,  'FontSize', 30);
xlabel(xl, 'FontSize', 30);
ylabel(yl, 'FontSize', 30);
xticks(linspace(0, llinfo.v_amb, 16));
lgd = legend;
lgd.FontSize = 10;
lgd.FontWeight = 'bold';

