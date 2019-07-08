function voxels = loadVoxelGridFromDistanceField(filename, resolution, shift, alpha)

data = importdata(filename);
voxels = zeros([resolution resolution resolution]);
voxels(:,:,:) = exp(-alpha);


for idx=1:size(data,1)
    x = data(idx,:);
    
    d = x(4);
    if d > alpha
        P = exp(-alpha);
    elseif d < -1
        P = exp(1);
    else
        P = exp(-d);
    end
    
    voxels(shift+x(1), shift+x(2), shift+x(3)) = P;
    
end

