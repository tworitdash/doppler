function [S] = Gauss_PDF(m0, mu, sigma, x)



    S = m0 .* exp(-(x - mu).^2./(2 * sigma.^2));



end