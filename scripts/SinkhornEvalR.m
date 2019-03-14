function out = SinkhornEvalR(v,w,a)
    out = v .* filt( w .* a);
end

