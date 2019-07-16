function out = filt3(a)
    global filter_size;
    global filter_padding_value;
    global preprocess_sigma;
    out = imgaussfilt3(a, preprocess_sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
end
