function [wasserstein_dist,v,w] = Sinkhorn(im1,im2)

global sigma;
global sinkhorn_iterations;
global sinkhorn_eps;

iter = sinkhorn_iterations;

filter = @(x) filt(x);
if ndims(im1) == 3
filter = @(x) filt3(x);
end

v = ones(size(im1));
w = ones(size(im1));

prev_wasserstein_dist = Inf;
change = Inf;
i=0;
stop = false;
while(i < iter && ~stop)
    i = i+1;
    w = maxv(filter(v), 10^-100);
    w = im2 ./ w;
    
    v = maxv(filter(w), 10^-100);
    v = im1 ./ v;
    
    wasserstein_dist = im1 .* log(maxv(v, 10^-100)) + im2 .* log(maxv(w, 10^-100));
    wasserstein_dist = sum(wasserstein_dist(:)) * sigma;
    
    if isfinite(prev_wasserstein_dist)
        change = abs(prev_wasserstein_dist - wasserstein_dist);
    end
    
    if change < sinkhorn_eps
        stop = true;
        disp("Convergence reached, no further iterations.");
    end
    
    prev_wasserstein_dist = wasserstein_dist;
    
    disp("("+i+"/"+iter+"), change = "+change)
end


end







