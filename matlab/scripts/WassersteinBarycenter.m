function im_out = WassersteinBarycenter(im1, im2, alpha)

    global sigma;
    
    filter = @(x) filt(x);
    if ndims(im1) == 3
    filter = @(x) filt3(x);
    end

    iter = 50;

    v1 = ones(size(im1));
    w1 = ones(size(im1));
    d1 = ones(size(im1));
    v2 = ones(size(im1));
    w2 = ones(size(im1));
    d2 = ones(size(im1));
    im_out = ones(size(im1));

    for j = 1:iter
        im_out = ones(size(im1));
        
        w1 = im1 ./ maxv(filter(v1), 10^-100);
        d1 = v1 .* filter(w1);
        im_out = im_out .* d1.^alpha;
        
        w2 = im2 ./ maxv(filter(v2), 10^-100);
        d2 = v2 .* filter(w2);
        im_out = im_out .* d2.^(1.0-alpha);
        
        v1 = v1 .* im_out ./ maxv(d1, 10^-100);
        v2 = v2 .* im_out ./ maxv(d2, 10^-100);
        
    end
    
end