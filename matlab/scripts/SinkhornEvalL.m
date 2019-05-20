function out = SinkhornEvalL(v,w,a)
    
    filter = @(x) filt(x);
    if ndims(v) == 3
        filter = @(x) filt3(x);
    end
    
    out = w .* filter( v .* a);
end

