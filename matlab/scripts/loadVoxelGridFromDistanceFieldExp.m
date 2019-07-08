function voxels = loadVoxelGridFromDistanceField(filename, resolution, shift, alpha, beta)

data = importdata(filename);
voxels = zeros([resolution resolution resolution]);
voxels(:,:,:) = exp(-alpha-beta);


for idx=1:size(data,1)
    x = data(idx,:);
    
    d = x(4);
    if d > alpha
        d = alpha;
    elseif d < -alpha
        d = -alpha;
    end
    P = exp(-d-beta);
    
    voxels(shift+x(1), shift+x(2), shift+x(3)) = P;
    
end

