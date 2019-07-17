function kernel = shiftedGaussian(p, sigma, kernel_size, resolution)

kernel = zeros(resolution, resolution, resolution);
kernel_size_half = (kernel_size-1)/2;
pr = round(p);

for x = pr(1)-kernel_size_half:pr(1)+kernel_size_half
    for y = pr(2)-kernel_size_half:pr(2)+kernel_size_half
        for z = pr(3)-kernel_size_half:pr(3)+kernel_size_half
            if x > 0 && y > 0 && z > 0 && x <= resolution && y <= resolution && z <= resolution
                kernel(x,y,z) = 1 / (sigma^3 * (2*pi)^(3/2)) *  exp(-((x-p(1))^2 + (y-p(2))^2 + (z-p(3))^2)/(2*sigma));
            end
        end
    end
end

end