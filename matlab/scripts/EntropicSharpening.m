function mu = EntropicSharpening(mu, H0)

H_mu = Entropy(mu);

beta = 1;
if(H_mu + sum(mu(:)) > H0 + 1)
    options = optimset('TolX',1e-7);
    fun = @(x) sum(mu.^x,'all') + Entropy(mu.^x) - (1+H0);
%     x = fminbnd(fun,0,10,options);
    x = fminsearch(fun,0,options); % technically fminbnd is correct, but fminsearch gives more consistently good results with fewer parameters
    beta = x;
end
mu = mu.^beta;
end