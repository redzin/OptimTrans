function [wasserstein_dist,v,w] = Sinkhorn(im1,im2)

global sigma;

iter = 250;

v = ones(size(im1,1), size(im1,2));
w = ones(size(im1,1), size(im1,2));

for i = 1:iter
    v = maxv(filt(w), 10^-300);
    v = im1 ./ v;

    w = maxv(filt(v), 10^-300);
    w = im2 ./ w;
end

wasserstein_dist = im1 .* log(maxv(v, 10^-300)) + im2 .* log(maxv(w, 10^-300));
wasserstein_dist = sum(wasserstein_dist(:)) * sigma;

end