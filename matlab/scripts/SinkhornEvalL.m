function out = SinkhornEvalL(v,w,a)
    out = w .* filt( v .* a);
end

