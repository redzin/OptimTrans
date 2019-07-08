function [wasserstein_dist,v,w] = Sinkhorn(im1,im2)

global sigma;
global sinkhorn_iterations;

iter = sinkhorn_iterations;

filter = @(x) filt(x);
if ndims(im1) == 3
filter = @(x) filt3(x);
end

v = ones(size(im1));
w = ones(size(im1));

for i = 1:iter
    w = maxv(filter(v), 10^-100);
    w = im2 ./ w;
    
    v = maxv(filter(w), 10^-100);
    v = im1 ./ v;
end

wasserstein_dist = im1 .* log(maxv(v, 10^-100)) + im2 .* log(maxv(w, 10^-100));
wasserstein_dist = sum(wasserstein_dist(:)) * sigma;

end







