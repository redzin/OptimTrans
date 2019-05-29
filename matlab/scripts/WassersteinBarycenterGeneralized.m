function im_out = WassersteinBarycenterGeneralized(imgs, alpha, entropic_factor)

    global sigma;
    global enable_entropic_sharpening;
    
    filter = @(x) filt(x);
    if ndims(imgs{1}) == 3
    filter = @(x) filt3(x);
    end

    iter = 50;
    lower_limit = 10^-300;
%     epsilon = 10^-15;
    
    alpha = alpha ./ sum(alpha(:));
    
    for i = 1:length(imgs)
        v{i} = ones(size(imgs{1}));
        w{i} = ones(size(imgs{1}));
        d{i} = ones(size(imgs{1}));
    end
    
    im_out = ones(size(imgs{1}));
    im_out_prev = ones(size(imgs{1}));
    
    if (enable_entropic_sharpening)
        H0 = Entropy(imgs{1}); % Entropic Sharpening
    end
    
    
    for j = 1:iter
        
        im_out_prev = im_out;
        im_out = ones(size(imgs{1}));
        
        for i=1:length(imgs)
            w{i} = imgs{i} ./maxv(filter(v{i}), lower_limit);
            d{i} = v{i} .* filter(w{i});
            im_out = im_out .* d{i}.^alpha(i);
            
            
            if (enable_entropic_sharpening)
                H = Entropy(imgs{i}); % Entropic Sharpening
                if (H > H0) % Entropic Sharpening
                    H0 = H*entropic_factor; % downscale slightly to impose entropic sharpening on the input for which H = H0
                end
            end
            
        end
        
        if (enable_entropic_sharpening)
            im_out = EntropicSharpening(im_out, H0); % Entropic Sharpening
        end
        
        for i=1:length(imgs)
            v{i} = v{i} .* im_out ./ maxv(d{i}, lower_limit);
        end
        
%         if(sum((im_out-im_out_prev).^2,'all') < epsilon)
%             j
%             break;
%         end
        if (isnan(sum(im_out(:))))
            j
            im_out = im_out_prev;
            break;
        end
    end
end