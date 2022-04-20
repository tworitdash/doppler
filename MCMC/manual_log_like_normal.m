function ret = manual_log_like_normal(x, data)
    ret = sum(-log(x(2) .* sqrt(2 * pi)) - ((data-x(1)).^2)./(2 * x(2).^2));
end