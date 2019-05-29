%% Compute the entropy the input
function H = Entropy(mu)

H = mu .* log(maxv(mu, 10^-100));
H = -sum(H(:));

end

