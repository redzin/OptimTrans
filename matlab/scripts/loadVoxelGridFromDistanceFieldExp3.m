function voxels = loadVoxelGridFromDistanceField(filename, resolution, shift)

data = importdata(filename);
voxels = zeros([resolution resolution resolution]);
voxels(:,:,:) = exp(-alpha);


for idx=1:size(data,1)
    x = data(idx,:);
    voxels(shift+x(1), shift+x(2), shift+x(3)) = exp(-x(4));
end

