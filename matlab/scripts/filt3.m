function out = filt3(a)
    global filter_size;
    global filter_padding_value;
    global sigma;
    out = imgaussfilt3(a, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
end
