function out = SinkhornEvalR(v,w,a)
    
    filter = @(x) filt(x);
    if ndims(v) == 3
        filter = @(x) filt3(x);
    end

    out = v .* filter( w .* a);
end

