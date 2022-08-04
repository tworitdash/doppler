function [qwjk, ind] = q1(yfft, wjk, freq)

[indxlow] = (freq(2:end) >= wjk);
[indxup] = (freq(1:end-1) <= wjk);

indx = (indxlow & indxup);

Ones = linspace(1, length(yfft), length(yfft));

ind = Ones(indx);

qwjk = yfft(indx);

end