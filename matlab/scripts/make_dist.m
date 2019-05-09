function out = make_dist(x)

out = maxv(x, 10^-300);
out = filt(out);
% out = maxv(out, 10^-300);
out = out ./ sum(out(:));

end

