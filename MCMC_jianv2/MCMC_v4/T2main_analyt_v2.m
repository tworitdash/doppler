clear;

Nu = 10000;
Mu = pi; 
Sigma = 0.5 * 4 * pi / 0.03 * 1e-3;
U_ = normrnd(Mu, Sigma, [1 Nu]);

% Psi_ = -pi + 2 * pi * rand([1 Nu]);
Psi_ = 0 .* rand([1 Nu]);

[Psim, Psin] = meshgrid(Psi_, Psi_);

[U1, U2] = meshgrid(U_, U_);

Np = 5;
Ng = 95;
Ns = 2;
K_ = [];

for i = 1:Ns
    K_ = [K_ (i-1)*(Np + Ng):(i-1)*(Np+Ng)+Np-1];
end

NToT = Ns * Np + Ng * (Ns - 1);
Const = -(1 + Nu) * (Ns-1) * Np;

FirstTint = - csc((U1 - U2)./2) .* sin((U1 - U2)/2 .* Np) .* csc((U1 - U2)./2 .* (Np + Ng))...
    .* sin((U1 - U2)./2 .* (Np + Ng) .* Ns) .* cos((U1 - U2)./2 .* NToT + (Psim - Psin));

FirstT = nansum(nansum(FirstTint)) + Nu;
Nup = 10000;
u_ = linspace(0, 2 * pi, Nup);
[Uw, w] = meshgrid(U_, u_);
[Psiw, ~] = meshgrid(Psi_, u_);

SecondTint = 2 .* csc((Uw - w)./2) .* sin((Uw - w)/2 .* Np) .* csc((Uw - w)./2 .* (Np + Ng))...
    .* sin((Uw - w)./2 .* (Np + Ng) .* Ns) .* cos((Uw - w)./2 .* NToT + Psiw);

SecondT = sum(SecondTint, 2);
LL = Const + FirstT + SecondT;

figure(1);
hold on; plot(u_, LL - max(LL));

% figure; histogram(U_);