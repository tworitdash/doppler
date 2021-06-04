C = rand(5, 4) + 1j .* rand(5, 4);

k = linspace(eps, 2*pi, 360);

dk = k(2) - k(1);
count = 0;

for o = 1:length(k)
    for i = 1:size(C, 2)
        for l = 1:size(C, 2)
            if l ~= i
%                 count = count + 1
                si(o) = sum((C(:, i) - C(:, l) .* exp(1j.*k(o)))./((C(:, i) - C(:, l) .* exp(1j.*k(o)))).^2);
            end
        end
        F(i) = sum(si .* dk);
    end
end



