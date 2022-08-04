% clear;
% close all;
% 
% 
% Mu = 7.5/4 .* 0.4189;
% Sigma = 0.3 .*  0.4189;
% 
% % rng(1, 'twister');
% Nu = 100000;
% 
% beta = 2 * pi .* rand([1 Nu]); 
% 
%     
% U = normrnd(Mu, Sigma, [1 Nu]);
% 
% 
% Nt = 1024;
% K = 0:1:Nt-1;
% % K = Nt;
% 
% x = linspace(-2*pi*Nt, 2*pi*Nt, 10000);
% dx = x(2) - x(1);
% 
% for k = 1:length(K)
%      Sigp = exp(1j .* K(k) .* U + 1j .* beta) ;
%      
%      
% %      if mod(k, 200) == 0 
%      if k == 2
%       figure(1000); hold on; histogram(real(Sigp));  figure(2000); hold on; histogram(imag(Sigp)); 
% %       figure(3000); hold on; plot(real(Sigp));  figure(4000); hold on; plot(imag(Sigp)); 
%         disp('\mu'); mean(real(Sigp))
%         disp('\sigma'); std(real(Sigp))
%         disp('sum of sigp'); sum(Sigp)
%      end
%      
%      sig(k) = sum( exp(1j .* K(k) .* U + 1j .* beta) );
%     
% %     figure(100); hold on; histogram(K(k) .* U + beta);
% %     Y = 1/2 .* sign(K(k)) .* sign(Sigma) .* cos(x)...
% %         .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) ));
% 
%      Y =  1/2 .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) ));
% %     figure(200); hold on; plot(x, Y);
%     
% %      Ycos(k) = sum(1/2 .* sign(K(k)) .* sign(Sigma) .* cos(x)...
% %          .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )));
% %      
% %      Ysin(k) = sum(1/2 .* sign(K(k)) .* sign(Sigma) .* sin(x)...
% %          .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )));
%     
%      Ycos(k) = Nu .* sum(1./2 .* (cos(x)) .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )) .* dx);
%      
%      Ysin(k) = Nu .* sum(1/2 .* (sin(x)) .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )) .* dx);
%     
%      Ycossq(k) = Nu.^2 .* sum(1./2 .* (cos(x)).^2 .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )) .* dx);
%      
%      Ysinsq(k) = Nu.^2 .* sum(1/2 .* (sin(x)).^2 .* (erf( (K(k).*Mu-x+2*pi)./(sqrt(2) .* K(k) .* Sigma) ) - erf( (K(k) .* Mu - x)./(sqrt(2) .* K(k) .* Sigma) )) .* dx);
%     
%      Ycosvar(k) = Ycossq(k) - (Ycos(k)).^2;
%      
%      Ysinvar(k) = Ysinsq(k) - (Ysin(k)).^2;
%     
% end
% 
% figure(300); hold on; plot(K, Ycos);
% figure(400); hold on; plot(K, Ysin);
% 
% figure(500); hold on; plot(K, Ycosvar);
% figure(600); hold on; plot(K, Ysinvar);
% 
% Mur = mean(real(sig));
% Sr = std(real(sig));
% 
% Mui = mean(imag(sig));
% Si = std(imag(sig));
% 
% 
% % X_cos = @(x) 1./(4 ) .* exp(-1./2 .* K .* (K .* Sigma.^2 + 2 .* 1j .* Mu)) .* ...
% %     ( 1j .* exp(1j .* 2 .* K .* Mu) .* ...
% %     erfz( (1j .* K.^2 .* Sigma.^2 + K .* (Mu - x) + 2 * pi)./(sqrt(2) .* K .* Sigma) ) - ...
% %     1j .* exp(1j .* 2 .* K .* Mu) .* erfz( (1j .* K .* Sigma.^2 + Mu - x)./(sqrt(2).*Sigma) ) - ...
% %     2 .* exp(1./2 .* K .* (K .* Sigma.^2 + 2 .* 1j .* Mu)) .* sin(K .* x) .* erf( (Mu - x)./(sqrt(2) .* Sigma) ) + ...
% %     2 .* exp(1/2 .* K .* (K .* Sigma.^2 + 2 .* 1j .* Mu)) .* sin(K .* x) .* erf( (K .* Mu - K .* x + 2 * pi)./(sqrt(2) .* K .* Sigma) ) - ...
% %     erfi( (K.^2 .* Sigma.^2 + 1j .* K .* (Mu - x) + 2 .* 1j.* pi)./(sqrt(2) .* K .* Sigma) ) + ...
% % %     erfi( (K.*Sigma.^2 + 1j .* Mu - 1j .* x)./(sqrt(2) .* Sigma) ));
% % 
% % X_cosm = @(x)  Nu * sign(Sigma) .* cos(K .* x) .* ( erf( (2*pi+Mu-x)./(sqrt(2) .* Sigma) ) - erf( (Mu - x)./(sqrt(2) .* Sigma) ) );
% % X_sinm =  @(x) Nu * sign(Sigma) .* sin(K .* x) .* ( erf( (2*pi+Mu-x)./(sqrt(2) .* Sigma) ) - erf( (Mu - x)./(sqrt(2) .* Sigma) ) );
% % 
% % X_cossq = @(x) Nu^2 *  sign(Sigma) .* (cos(K .* x)).^2 .* ( erf( (2*pi+Mu-x)./(sqrt(2) .* Sigma) ) - erf( (Mu - x)./(sqrt(2) .* Sigma) ) );
% % X_sinsq =  @(x) Nu^2 * sign(Sigma) .* (sin(K .* x)).^2 .* ( erf( (2*pi+Mu-x)./(sqrt(2) .* Sigma) ) - erf( (Mu - x)./(sqrt(2) .* Sigma) ) );
% % 
% % X_cosvar = @(x) X_cossq(x) - (X_cosm(x)).^2;
% % X_sinvar = @(x) X_sinsq(x) - (X_sinm(x)).^2;
% % 
% % figure; plot( X_cosvar(2*pi) - X_cosvar(0) );
% % figure; plot( X_sinvar(2*pi) - X_sinvar(0) );
% % 
% % figure; plot( X_cosm(2*pi) - X_cosvar(0) );
% % figure; plot( X_sinm(2*pi) - X_sinvar(0) );
% 

%%  
clear;
% close all;

lambda = 0.03;
dt = 1e-3;

Mu = 7.5 .* 0.4189;
Sigma = 0.2 .*  0.4189;

% rng(1, 'twister');
Nu = 100000;

beta = 2 * pi .* rand([1 Nu]); 

    
U = normrnd(Mu, Sigma, [1 Nu]);


Nt = 64;
K = 0:1:Nt-1;
% K = Nt

x = linspace(-2*pi*Nt, 2*pi*Nt, 10000);
dx = x(2) - x(1);

for k = 1:length(K)
    
     if mod(K(k), 10) == 0
      figure(1000); hold on; histogram(K(k) .* U + beta);
     end
     
     sig(k) = sum( exp(1j .* K(k) .* U + 1j .* beta) );
%      sig(k) = sum( exp(1j .* K(k) .* U ) );
%     
    
end

figure; plot(real(sig)); hold on; plot(imag(sig));

figure(100); hold on; histogram(real(sig)); figure(200); hold on; histogram(imag(sig));
% hold on; plot(Nu .* Mu .* exp(-Sigma.^2 .* K.^2/2) .* cos(Mu)); 
% hold on; plot(Nu .* Mu .* exp(-Sigma.^2 .* K.^2/2) .* sin(Mu)); 

Nfft = 1024;

v = linspace(0, lambda/(2 * dt), Nfft);
figure; plot(v, db(abs(fft(sig, Nfft))) )

