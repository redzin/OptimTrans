function voxels = loadVoxelGridFromDistanceField(filename, resolution, shift)

data = importdata(filename);

voxels = zeros([resolution resolution resolution]);

for idx=1:size(data,1)
    x = data(idx,:)+shift;
    voxels(x(1), x(2), x(3)) = 1;
end