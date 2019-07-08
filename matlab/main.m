%% Global Settings

close all
clear all
clc
format long

addpath('scripts')
addpath('libraries/WOBJ_toolbox_Version2b')
addpath('libraries/export_fig')
addpath('../data/objects')
addpath('../data/images')
addpath('../data/cochleas')

global sigma;
global filter_size;
global filter_padding_value;
global enable_entropic_sharpening;
global enable_prints;
global sinkhorn_iterations;

sinkhorn_iterations = 50;
sigma = 1.05;
filter_size = 129;
filter_padding_value =  0.0;

sigmas = [1.0]; %[1.0 2.5 4.1 5.0];
filter_sizes = [129]; %[5 41 45 101];
% print_types = {{"-dpng", ".png"}};
print_types = {{"-depsc", ".eps"}, {"-dpng", ".png"}};
enable_prints = true;

% For Wasserstein Barycenter computaions
enable_entropic_sharpening = true;

% Cochlea files


cochlea_folder = "../data/cochleas/";
cochlea_files = [
    "shape05_pca",
    "shape06_pca",
    "shape08_pca",
    "shape09_pca",
    "shape10_pca",
    "shape11_pca",
    "shape12_pca",
    "shape15_pca",
    "shape16_pca",
    "shape18_pca",
    "shape19_pca",
    "shape20_pca",
    "shape21_pca",
    "shape22_pca",
    "shape23_pca",
    "shape24_pca",
    "shape5876_pca",
    "shape6317_pca"
];

%% Wasserstein distance test
format long g


imgs = {
         im2double(rgb2gray(imread('../images/one-blob.png')))
         im2double(rgb2gray(imread('../images/one-blob-moved.png')))
         im2double(rgb2gray(imread('../images/one-blob-moved-even-more.png')))
         im2double(rgb2gray(imread('../images/one-blob-moved-even-more-again.png')))
};


for filter_size=filter_sizes
    for sigma=sigmas
        
        disp('----------------')
        disp("Meta: sigma = "+sigma+", filter_size = "+filter_size)
        
        dists = {
                make_dist(imgs{1}),
                make_dist(imgs{2}),
                make_dist(imgs{3}),
                make_dist(imgs{4})
            };
            
        FigH = figure('Position', [0, 0, 350, 120]);
        colormap('gray')
        
        subplot(1,4,1)
        hold on
        axis equal
        axis([0, 64, 0 ,64])
        imagesc(dists{1})
        title('Original')
        
        
        for i = 2:4
            [wd,v,w] = Sinkhorn(dists{1},dists{i});
            ed = sqrt(sum(  (dists{1}(:) - dists{i}(:)).^2  ));
            marg = SinkhornEvalR(v,w,ones(size(v,1),size(v,2)));
        %     wd
        
            disp("#"+(i-1)+": "+"W2 = " + wd + ", D2 = " + ed + ", sum(\pi1) = " + sum(marg(:)));
            subplot(1,4,i)
            hold on
            axis equal
            axis([0, 64, 0, 64])
            imagesc(dists{i})
            title("#"+(i-1))
        
        %     title("\pi1 = p. sum(\pi1) = " + sum(img(:)))
        %     title("W2 = " + wd + ", D2 = " + ed)
        
        end
        
        if enable_prints
            for type = print_types
                print("prints/delta-comparison-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2},type{1}{1})
            end
        end

    end
end


disp('------------------------------')
disp('------------------------------')
disp('------------------------------')

%% Coupling Test #1

format long g

imgs = {
         im2double(rgb2gray(imread('images/one-blob.png')))
         im2double(rgb2gray(imread('images/one-blob-moved.png')))
         im2double(rgb2gray(imread('images/one-blob-moved-even-more.png')))
         im2double(rgb2gray(imread('images/one-blob-moved-even-more-again.png')))
};



for filter_size=filter_sizes
    for sigma=sigmas
        disp('-----------------------------')
        disp("Meta: sigma = "+sigma+", filter_size = "+filter_size)
        
        dists = {
                make_dist(imgs{1}),
                make_dist(imgs{2}),
                make_dist(imgs{3}),
                make_dist(imgs{4})
            };
        
        for i = 2:4
            figure('Position', [0, 0, 200, 200]);
            colormap('gray')
            
            subplot(2,2,1)
            hold on
            axis equal
            axis([0, 64, 0, 64])
            imagesc(dists{1})
            title('p')
        
            [wd,v,w] = Sinkhorn(dists{1},dists{i});
            ed = sqrt(sum(  (dists{1}(:) - dists{i}(:)).^2  ));
            subplot(2,2,3)
            hold on
            axis equal
            axis([0, 64, 0, 64])
            imagesc(dists{i})
            title('q')
        %     title("W2 = " + wd + ", D2 = " + ed)
            
            img = SinkhornEvalR(v,w,ones(size(v,1),size(v,2)));
            subplot(2,2,2)
            hold on
            axis equal
            axis([0, 64, 0, 64])
            imagesc(img)
        %     title('\pi1')
            title("\pi1, sum = " + sum(img(:)))
            
            img = SinkhornEvalL(v,w,ones(size(v,1),size(v,2)));
            subplot(2,2,4)
            hold on
            axis equal
            axis([0, 64, 0, 64])
            imagesc(img)
        %     title('\pi^T1')
            title("\pi^T1, sum = " + sum(img(:)))
            
            if enable_prints
                for type = print_types
                    print("prints/delta-comparison-with-marginals-"+sigma+"-sigma-"+filter_size+"-filter-size-"+i+type{1}{2},type{1}{1})
                end
            end
        end
    end
end

%% Coupling Test #2

format long g

im1 = im2double(rgb2gray(imread('images/smile_s.png')));
im2 = im2double(rgb2gray(imread('images/dots_s.png')));

for filter_size=filter_sizes
    for sigma=sigmas
        disp('-----------------------------')
        disp("Meta: sigma = "+sigma+", filter_size = "+filter_size)
        
        dist1 = make_dist(im1);
        dist2 = make_dist(im2);
        
        [wd,v,w] = Sinkhorn(dist1, dist2);
        disp("Wasserstein distance = " + wd)
        
        figure()
        figure('Position', [0, 0, 200, 200]);
        colormap('gray')
            
        hold on
        subplot(2,2,1);
        axis equal
        axis([0 64 0 64])
        imagesc(dist1)
        title("p")
        
        subplot(2,2,3);
        axis equal
        axis([0 64 0 64])
        imagesc(dist2)
        title("q")
        
        img = SinkhornEvalR(v,w,ones(size(v,1),size(v,2)));
        subplot(2,2,2);
        axis equal
        
        axis([0 64 0 64])
        imagesc(img)
        title("\pi1, sum = " + sum(img(:)))
        
        img = SinkhornEvalL(v,w,ones(size(v,1),size(v,2)));
        subplot(2,2,4)
        axis equal
        axis([0 64 0 64])
        imagesc(img)
        title("\pi^T1, sum = " + sum(img(:)))
        
        if enable_prints
            for type = print_types
                print("prints/smiley-comparison-"+sigma+"-sigma-"+filter_size+"-filter-size-with-marginals"+type{1}{2},type{1}{1})
            end
        end
    end
end

%% Coupling Test #3

format long g

dist1 = make_dist(im2double(rgb2gray(imread('images/dots_s.png'))));
dist2 = make_dist(im2double(rgb2gray(imread('images/smile_s.png'))));

disp('-----------------------------')
[wd,v,w] = Sinkhorn(dist1, dist2);
disp("Wasserstein distance = " + wd)

p = [
    linspace(35,35,64)
    1:64
];
mask = ones(size(v,1), size(v,2));
% b = boxes(size(v,1), size(v,2));
% mask = b{1};
% mask = zeros(size(v,1), size(v,2));
% mask(sub2ind(size(mask), p(1,:), p(2,:))) = 1;
img = SinkhornEvalL(v,w,mask);
composition = [mask img];
composition = composition > 0;


dists = [dist1 dist2];
figure()
hold on
axis equal
imagesc([dist1.*mask img])

figure()
hold on
axis equal
imagesc([dist1 dist2])


%% Quadrant test

format long g

dist1 = make_dist(im2double(rgb2gray(imread('images/smile_s.png'))));
dist2 = make_dist(im2double(rgb2gray(imread('images/dots_s.png'))));


% figure()
% hold on
% axis equal
% imshowpair(imcomplement(imresize(imread('images/dots_s.png'),2)),imcomplement(imresize(imread('images/smile_s.png'),2)), 'montage')

[wd,v,w] = Sinkhorn(dist1, dist2);
disp('---------------------------------------')
disp("Wasserstein distance = " + wd)

img_out = zeros(size(dist1,1),size(dist1,2),4);
cmp_out = zeros(size(dist1,1),size(dist1,2),4);
s = 0;
for b = boxes(size(dist1,1), size(dist1,2));
    
    b = b{1};
    
    % Get random color
    color = randn(1,3);
    color = color .* color;
    color = color / max(color);
    
    img = SinkhornEvalR(v,w,b);
    s = s + sum(img(:));
    col_img = zeros(size(img,1), size(img,2), 3);
    
    alpha = minv(img(1:size(img,1),1:size(img,2)) / 0.002, 1.0);
    col_img(1:size(img,1),1:size(img,2),1) = color(1) * alpha;
    col_img(1:size(img,1),1:size(img,2),2) = color(2) * alpha;
    col_img(1:size(img,1),1:size(img,2),3) = color(3) * alpha;
    col_img(1:size(img,1),1:size(img,2),4) = alpha;
    
    mask = img_out(:,:,4) < alpha;
    img_out(mask(:,:,[1,1,1,1])) = col_img(mask(:,:,[1,1,1,1]));
    
    cmp_tmp = b .* dist2;
    cmp_img = zeros(size(cmp_tmp,1), size(cmp_tmp,2), 3);
    alpha = minv(cmp_tmp(1:size(cmp_tmp,1),1:size(cmp_tmp,2)) / 0.002, 1.0);
    cmp_img(1:size(cmp_tmp,1),1:size(cmp_tmp,2),1) = color(1) * alpha;
    cmp_img(1:size(cmp_tmp,1),1:size(cmp_tmp,2),2) = color(2) * alpha;
    cmp_img(1:size(cmp_tmp,1),1:size(cmp_tmp,2),3) = color(3) * alpha;
    cmp_img(1:size(cmp_tmp,1),1:size(cmp_tmp,2),4) = alpha;
    
    mask = cmp_out(:,:,4) < alpha;
    cmp_out(mask(:,:,[1,1,1,1])) = cmp_img(mask(:,:,[1,1,1,1]));
    
    
end

disp("Marginal sum = " + s)

figure()
hold on
axis equal
imshowpair(imcomplement(imresize(cmp_out(:,:,1:3),2)),imcomplement(imresize(img_out(:,:,1:3),2)), 'montage')

if enable_prints
    for type = print_types
        print("prints/quadrant-test-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1})
    end
end


%% 3d Basic Setup Test
format long g

close all
clc

showE = 1;
renewVars = 0;


if(renewVars | exist('orig_A','var') ~= 1)
orig_A = rand([64,64,64]);
orig_A(orig_A > 0.99) = 1;
orig_A(orig_A < 0.99) = 0;
orig_B = zeros([64,64,64]);
orig_B(16,16,48) = 1;

orig_C = zeros([64,64,64]);
orig_C(32,32,32) = 1;

orig_D = zeros([64,64,64]);
orig_D(48,48,16) = 1;

orig_E = zeros([64,64,64]); % Only used for display purposes
orig_E(orig_B > 0 | orig_C > 0 | orig_D > 0) = 1;
end
A = orig_A;
B = orig_B;
C = orig_C;
D = orig_D;
E = orig_E;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(orig_A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test1-A-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

if(showE)
    figure
    volshow(orig_E,...
            'Renderer', 'MaximumintensityProjection',...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
    if enable_prints
        set(gcf,'PaperPositionMode','auto')
        for type = print_types
            print("prints/3d-test1-E-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
        end
    end
end

% Gaussian filter pre-processing
A = imgaussfilt3(A, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
A = A ./ sum(A(:));
B = imgaussfilt3(B, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
B = B ./ sum(B(:));
C = imgaussfilt3(C, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
C = C ./ sum(C(:));
D = imgaussfilt3(D, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
D = D ./ sum(D(:));
E = imgaussfilt3(E, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
E = E ./ sum(E(:));

% Post-Gaussian-filter display
figure
volshow(A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);

if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test1-A-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

if(showE)
    figure
    volshow(E,...
            'Renderer', 'MaximumintensityProjection',...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
    if enable_prints
        set(gcf,'PaperPositionMode','auto')
        for type = print_types
            print("prints/3d-test1-E-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
        end
    end
end

% Print meta info
disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% A-B Comparison
[wd,v,w] = Sinkhorn(A,B);
ed = sqrt(sum(  (orig_A(:) - orig_B(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-B: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

% A-C Comparison
[wd,v,w] = Sinkhorn(A,C);
ed = sqrt(sum(  (orig_A(:) - orig_C(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-C: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

% A-D Comparison
[wd,v,w] = Sinkhorn(A,D);
ed = sqrt(sum(  (orig_A(:) - orig_D(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-D: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))


%% 3d Basic Setup Test 2
format long g

close all
clc

showE = 1;
A = orig_A;
B = orig_B;
C = orig_C;
D = orig_D;
E = orig_E;

A(1:64,1:64,1:32) = 0;
A(32:64,1:64,1:64) = 0;
A(1:64,32:64,1:64) = 0;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test2-A-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
if(showE)
    figure
    volshow(orig_E,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
    if enable_prints
        set(gcf,'PaperPositionMode','auto')
        for type = print_types
            print("prints/3d-test2-E-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
        end
    end
end    

% Gaussian filter pre-processing
A = imgaussfilt3(A, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
A = A ./ sum(A(:));
B = imgaussfilt3(B, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
B = B ./ sum(B(:));
C = imgaussfilt3(C, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
C = C ./ sum(C(:));
D = imgaussfilt3(D, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
D = D ./ sum(D(:));
E = imgaussfilt3(E, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
E = E ./ sum(E(:));

% Post-Gaussian-filter display
figure
volshow(A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test2-A-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
if(showE)
    figure
    volshow(E,...
            'Renderer', 'MaximumintensityProjection',...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
    if enable_prints
        set(gcf,'PaperPositionMode','auto')
        for type = print_types
            print("prints/3d-test2-E-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
        end
    end
end

% Print meta info
disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% A-B Comparison
[wd,v,w] = Sinkhorn(A,B);
ed = sqrt(sum(  (orig_A(:) - orig_B(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-B: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

% A-C Comparison
[wd,v,w] = Sinkhorn(A,C);
ed = sqrt(sum(  (orig_A(:) - orig_C(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-C: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

% A-D Comparison
[wd,v,w] = Sinkhorn(A,D);
ed = sqrt(sum(  (orig_A(:) - orig_D(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("A-D: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))


%% 3d Basic Setup Test, Display Marginals and/or Parameter Sensitivity

format long g

close all
clc

renewVars = 0;

if (renewVars)
X = rand([64,64,64]);
X(X > 0.9999) = 1;
X(X < 0.9999) = 0;
orig_X = X;
Y = rand([64,64,64]);
Y(Y > 0.9999) = 1;
Y(Y < 0.9999) = 0;
orig_Y = Y;
end
X = orig_X;
Y = orig_Y;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(orig_X,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-pre-gaussian-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
figure
volshow(orig_Y,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-pre-gaussian-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

for sigma=sigmas
X = orig_X;
Y = orig_Y;


% Gaussian filter pre-processing
X = imgaussfilt3(X, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
X = X ./ sum(X(:));
Y = imgaussfilt3(Y, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
Y = Y ./ sum(Y(:));

% Post-Gaussian-filter display
figure
volshow(X,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
figure
volshow(Y,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

% Print meta info
disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% A-B Comparison
[wd,v,w] = Sinkhorn(X,Y);
ed = sqrt(sum(  (orig_X(:) - orig_Y(:)).^2  ));
margR = SinkhornEvalR(v,w,ones(size(v)));
margL = SinkhornEvalL(v,w,ones(size(v)));
disp("A-B: WD = " + wd + ", ED = " + ed + ", sum of right marginal = " + sum(margR(:)) + ", sum of left marginal = " + sum(margL(:)) )

figure
volshow(margR,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

figure
volshow(margL,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

end

%% 3d Pulse Test

format long g

close all
clc

showE = 1;
renewVars = 1;


if(renewVars | exist('orig_A','var') ~= 1)
orig_B = zeros([64,64,64]);
orig_B(16,16,48) = 1;

orig_C = zeros([64,64,64]);
orig_C(32,32,32) = 1;

orig_D = zeros([64,64,64]);
orig_D(48,48,16) = 1;

orig_E = zeros([64,64,64]); % Only used for display purposes
orig_E(orig_B > 0 | orig_C > 0 | orig_D > 0) = 1;
end
B = orig_B;
C = orig_C;
D = orig_D;
E = orig_E;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
% figure
% volshow(orig_A,...
%         'Renderer', 'MaximumintensityProjection',...
%         'Colormap', color,...
%             'CameraTarget',[0 0 0],...
%             'CameraViewAngle',30,...
%             'CameraUpVector',[0 1 0],...
%             'CameraPosition',[2 1.5 2]);
% if enable_prints
%     set(gcf,'PaperPositionMode','auto')
%     for type = print_types
%         print("prints/3d-test1-A-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
%     end
% end

if(showE)
    figure
    volshow(orig_E,...
            'Renderer', 'MaximumintensityProjection',...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
%     if enable_prints
%         set(gcf,'PaperPositionMode','auto')
%         for type = print_types
%             print("prints/3d-test1-E-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
%         end
%     end
end

% Gaussian filter pre-processing
% A = imgaussfilt3(A, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
% A = A ./ sum(A(:));
B = imgaussfilt3(B, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
B = B ./ sum(B(:));
C = imgaussfilt3(C, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
C = C ./ sum(C(:));
D = imgaussfilt3(D, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
D = D ./ sum(D(:));
E = imgaussfilt3(E, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
E = E ./ sum(E(:));

% Post-Gaussian-filter display
% figure
% volshow(A,...
%         'Renderer', 'MaximumintensityProjection',...
%         'Colormap', color,...
%             'CameraTarget',[0 0 0],...
%             'CameraViewAngle',30,...
%             'CameraUpVector',[0 1 0],...
%             'CameraPosition',[2 1.5 2]);
% 
% if enable_prints
%     set(gcf,'PaperPositionMode','auto')
%     for type = print_types
%         print("prints/3d-test1-A-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
%     end
% end

if(showE)
    figure
    volshow(E,...
            'Renderer', 'MaximumintensityProjection',...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
%     if enable_prints
%         set(gcf,'PaperPositionMode','auto')
%         for type = print_types
%             print("prints/3d-test1-E-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
%         end
%     end
end

% Print meta info
disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% B-C Comparison
[wd,v,w] = Sinkhorn(B,C);
ed = sqrt(sum(  (orig_B(:) - orig_C(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("B-C: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

% B-D Comparison
[wd,v,w] = Sinkhorn(B,D);
ed = sqrt(sum(  (orig_B(:) - orig_D(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("B-D: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))

%% Finer grid time test

% Test notes:
% 240x240x240 seems to be max resolution before crash on desktop with 8 gigs of ram and ubuntu 19.04
% 200x200x200 takes about 30 seconds
% 64x64x64 is <10 seconds.


format long g

close all
clc

renewVars = 1;
arraySize = 64;
enable_prints = 0;
tic
if (renewVars)
X = rand([arraySize,arraySize,arraySize]);
X(X > 0.9999) = 1;
X(X < 0.9999) = 0;
orig_X = X;
Y = rand([arraySize,arraySize,arraySize]);
Y(Y > 0.9999) = 1;
Y(Y < 0.9999) = 0;
orig_Y = Y;
end
X = orig_X;
Y = orig_Y;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(orig_X,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-pre-gaussian-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
figure
volshow(orig_Y,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-pre-gaussian-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

for sigma=sigmas
X = orig_X;
Y = orig_Y;


% Gaussian filter pre-processing
X = imgaussfilt3(X, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
X = X ./ sum(X(:));
Y = imgaussfilt3(Y, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
Y = Y ./ sum(Y(:));

% Post-Gaussian-filter display
figure
volshow(X,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
figure
volshow(Y,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

% Print meta info
disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% A-B Comparison
[wd,v,w] = Sinkhorn(X,Y);
ed = sqrt(sum(  (orig_X(:) - orig_Y(:)).^2  ));
margR = SinkhornEvalR(v,w,ones(size(v)));
margL = SinkhornEvalL(v,w,ones(size(v)));
disp("A-B: WD = " + wd + ", ED = " + ed + ", sum of right marginal = " + sum(margR(:)) + ", sum of left marginal = " + sum(margL(:)) )

figure
volshow(margR,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

figure
volshow(margL,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

end

disp(toc)

%% Eternally repeating 2d Test of Interpolation using Wasserstein Barycenters

close all
clc

filename = 'animations/2d-barycenter-interpolation-example-01.gif';
frameRate = 60;
alphas = [linspace(0,0,2) linspace(0,1,6) linspace(1,1,2) linspace(1,0,6)];

dist1 = make_dist(im2double(rgb2gray(imread('../images/smile_s.png'))));
dist2 = make_dist(im2double(rgb2gray(imread('../images/dots_s.png'))));

fig = figure;
while(true)
for alpha = alphas
    barycenter = WassersteinBarycenter(dist1, dist2, alpha);
    img = barycenter;% ./ max(barycenter(:));
    imagesc(img);
    title(['sum = ', num2str(sum(barycenter(:)), 4)])
    drawnow
end
end


%% 2d Test of Interpolation using Wasserstein Barycenters

close all
clc

enable_entropic_sharpening = 1;
entropic_factor = 0.95;
write_gif = 0;
filename = 'animations/2d-barycenter-interpolation-example-01.gif';
frameRate = 60;
alphas = [linspace(0,0,20) linspace(0,1,60) linspace(1,1,20) linspace(1,0,60)];
% alphas = linspace(0,1,5);

dist1 = make_dist(im2double(rgb2gray(imread('../images/smile_s.png'))));
dist2 = make_dist(im2double(rgb2gray(imread('../images/dots_s.png'))));

fig = figure;
n = 1;
disp('Creating gif...')
for alpha = alphas
    if mod(n,10) == 0
        disp([num2str(round(100 * n / length(alphas))), '%...'])
    end
    
    barycenter = WassersteinBarycenterGeneralized({dist1, dist2}, [alpha, 1-alpha], entropic_factor);
    img = barycenter;% ./ max(barycenter(:));
    imagesc(img);
    title(['sum = ', num2str(sum(barycenter(:)), 4)])
    drawnow
    
    % Write to gif
    if (write_gif)
        frame = getframe(fig);
        [im,map] = frame2im(frame);

        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'DelayTime', 1/frameRate); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1/frameRate);
        end
    end
    n = n+1;
end
clc
close all
disp('Finished creating gif.')


%% 3d Test of Interpolation using Wasserstein Barycenters

format long g

close all
clc

renewVars = 0;
color = jet(256);
makeGif = 0;
frameRate = 30;
renderingType = 'Isosurface';
filename = ['animations/3d-barycenter-interpolation-example-',lower(renderingType),'-01.gif'];
alphas = [linspace(0,0,10) linspace(0,1,60) linspace(1,1,10)];

if(renewVars | exist('orig_A','var') ~= 1)
orig_A = rand([64,64,64]);
orig_A(orig_A > 0.99) = 1;
orig_A(orig_A < 0.99) = 0;

orig_B = zeros([64,64,64]);
orig_B(16,16,48) = 1;
end

A = orig_A;
B = orig_B;

% Gaussian filter pre-processing
A = imgaussfilt3(A, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
A = A ./ sum(A(:));
B = imgaussfilt3(B, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
B = B ./ sum(B(:));

if (~makeGif)
    figure
    volshow(A,...
        'Renderer', renderingType,...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
    figure
    volshow(B,...
        'Renderer', renderingType,...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);

end

fig = figure;
n = 1;
if (makeGif)
    disp('Creating gif...')
end
for alpha = alphas
    if mod(n,10) == 0
        disp([num2str(round(100 * n / length(alphas))), '%...'])
    end
    
    barycenter = WassersteinBarycenter(A, B, alpha);
    img = barycenter;% ./ max(barycenter(:));
    
    volshow(img,...
        'Renderer', renderingType,...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);

    title(['sum = ', num2str(sum(barycenter(:)), 4)])
    drawnow
    
    if (makeGif)
        % Write to gif
        frame = getframe(fig);
        [im,map] = frame2im(frame);

        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'DelayTime', 1/frameRate); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1/frameRate);
        end
    end
    n = n+1;
end
clc
close all
if(makeGif)
    disp('Finished creating gif.')
end


% Print meta info
% disp("Sigma: " + sigma + ", kernel size: "+ filter_size)

% A-B Comparison
% [wd,v,w] = Sinkhorn(A,B);
% ed = sqrt(sum(  (orig_A(:) - orig_B(:)).^2  ));
% marg = SinkhornEvalR(v,w,ones(size(v)));
% disp("A-B: WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))


%% 3d Test of Interpolation using Wasserstein Barycenters 2

format long g

close all
clc

renewVars = 0;
color = jet(256);
createGif = 1;

if(renewVars | exist('orig_A','var') ~= 1)
orig_A = rand([64,64,64]);
orig_A(orig_A > 0.99) = 1;
orig_A(orig_A < 0.99) = 0;
orig_A(1:64,1:64,1:32) = 0;
orig_A(32:64,1:64,1:64) = 0;
orig_A(1:64,32:64,1:64) = 0;

orig_B = rand([64,64,64]);
orig_B(orig_B > 0.99) = 1;
orig_B(orig_B < 0.99) = 0;
orig_B(1:64,1:64,32:64) = 0;
orig_B(1:32,1:64,1:64) = 0;
orig_B(1:64,1:32,1:64) = 0;
end

A = orig_A;
B = orig_B;

% Gaussian filter pre-processing
A = imgaussfilt3(A, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
A = A ./ sum(A(:));
B = imgaussfilt3(B, sigma, 'FilterSize', filter_size, 'Padding', filter_padding_value);
B = B ./ sum(B(:));


filename = 'animations/3d-barycenter-interpolation-example-02.gif';
frameRate = 60;
alphas = [linspace(0,0,20) linspace(0,1,60) linspace(1,1,20) linspace(1,0,60)];
% alphas = linspace(0,1,5);


if (~createGif)
figure
volshow(A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);

figure
volshow(B,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);
end


if (createGif)
fig = figure;
n = 1;
disp('Creating gif...')
for alpha = alphas
    if mod(n,10) == 0
        disp([num2str(round(100 * n / length(alphas))), '%...'])
    end
    
    barycenter = WassersteinBarycenter(A, B, alpha);
    img = barycenter;% ./ max(barycenter(:));
    
    volshow(img,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 1 0],...
            'CameraPosition',[2 1.5 2]);

    title(['sum = ', num2str(sum(barycenter(:)), 4)])
    drawnow
    
    % Write to gif
    frame = getframe(fig);
    [im,map] = frame2im(frame);
    
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'DelayTime', 1/frameRate); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1/frameRate);
    end
    n = n+1;
end
clc
close all
disp('Finished creating gif.')
end


%% Barycenter animation
close all
clc

sigma = 2.5;
enable_entropic_sharpening = 1;
entropic_factor = 0.99;
figure_1_name = 'bunny';
figure_2_name = 'duck';
makeGif = 1;
frameRate = 60;
renderingType = 'VolumeRendering';
alphas = [linspace(0,0,30) linspace(0,1,50) linspace(1,1,30) linspace(1,0,50)];
entropic_factor = 0.99;
grid_size = 64;
grid_shift = 4;
filename = ['animations/',figure_1_name,'-',figure_2_name,'-transformation-',lower(renderingType),'-grid-',num2str(grid_size),'.gif'];
color = [linspace(1,1,256)' linspace(0.8,0.8,256)' linspace(0.1,0.1,256)'];
bgcolor = [1 1 1];

A_data = importdata(['../data/objects/',figure_1_name,'_56.txt']);
B_data = importdata(['../data/objects/',figure_2_name,'_56.txt']);

orig_A = zeros([grid_size grid_size grid_size]);
for idx=1:size(A_data,1)
    x = A_data(idx,:)+grid_shift;
    orig_A(x(3), x(1), x(2)) = 1;
end
orig_A = flip(orig_A,3);
A = orig_A;

orig_B = zeros([grid_size grid_size grid_size]);
for idx=1:size(B_data,1)
    x = B_data(idx,:)+grid_shift;
    orig_B(x(3), x(2), x(1)) = 1;
end
orig_B = flip(orig_B,3);
B = orig_B;


A = filt3(A);
A = A ./ sum(A(:));
B = filt3(B);
B = B ./ sum(B(:));

fig = figure;
color1 = [linspace(1,1,256)' linspace(0.65,0.65,256)' linspace(1,1,256)'];
color2 = [linspace(1,1,256)' linspace(0.8,0.8,256)' linspace(0.1,0.1,256)'];
n = 1;
if (makeGif)
    disp('Creating gif...')
end
for alpha = alphas
    
    img = WassersteinBarycenterGeneralized({A, B}, [alpha 1-alpha], entropic_factor);
    color = alpha*color1 + (1-alpha)*color2;
    
    volshow(img,...
        'Renderer', renderingType,...
        'Colormap', color,...
        'BackgroundColor', bgcolor,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 110 0],...
        'CameraPosition',[2 1.5 2]);

    disp(['Alpha = ', num2str(alpha), ', sum = ', num2str(sum(img(:)))])
    drawnow
    
    
    if (makeGif)
        % Write to gif
        frame = getframe(fig);
        [im,map] = frame2im(frame);

        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'DelayTime', 1/frameRate); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 1/frameRate);
        end
    end
    
    if mod(n,10) == 0
        disp([num2str(round(100 * n / length(alphas))), '%...'])
    end
    
    n = n+1;
end

% close all
if(makeGif)
    disp('Finished creating gif.')
end

%% Barycenter 2d figure progression

close all
clc

sigma = 1.05;
color = jet(256);
renderingType = 'VolumeRendering';
filename = ['animations/teapot-duck-transformation-',lower(renderingType),'-02.gif'];
alphas = linspace(0,1,5);


A_data = importdata('../data/objects/bunny.txt');
B_data = importdata('../data/objects/duck.txt');

orig_A = zeros([32 32 32]);
for idx=1:size(A_data,1)
    x = A_data(idx,:)+8;
    orig_A(x(3), x(1), x(2)) = 1;
end
orig_A = flip(orig_A,3);
A = orig_A;

orig_B = zeros([32 32 32]);
for idx=1:size(B_data,1)
    x = B_data(idx,:)+8;
    orig_B(x(3), x(2), x(1)) = 1;
end
orig_B = flip(orig_B,3);
B = orig_B;


A = filt3(A);
A = A ./ sum(A(:));
B = filt3(B);
B = B ./ sum(B(:));


height = 250;
width = height * length(alphas)+1;

fig = figure('Position', [0 0 width height]);
hold on
n = 0;
for alpha=alphas
    barycenter = WassersteinBarycenter(A, B, alpha);
    img = barycenter;% ./ max(barycenter(:));
    
    width = 1/length(alphas);
    p = uipanel(fig, 'Position', [width*n, 0, width, 1], 'BorderType', 'none');
    
    volshow(img,...
        'Parent', p,...
        'Renderer', renderingType,...
        'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 110 0],...
            'CameraPosition',[1.7 1.5 1.7]);
    n = n + 1;
end

% if enable_prints
%     set(gcf,'PaperPositionMode','auto')
%     for type = print_types
%         print("prints/duck-teapot-barycenters"+type{1}{2}, type{1}{1}, "-r0")
%     end
% end




%% 3-way Barycenter test


close all
clc

sigma = 2.5;
enable_prints = false;
enable_entropic_sharpening = true;
entropic_factor = 0.99;
renderingType = 'VolumeRendering';
isosurface_value = 0.1;
element_height = 120;
color = jet(256);
% bgcolor = [0.3 0.75 0.93]; % blue
bgcolor = [1 1 1]; % white
alphamap = linspace(0,1,256)';
alphas = 5;
num_rows = 5;
grid_size = 64;
grid_displace = 4;
filename = "3-way-barycenter-interpolation-duck-teapot-bunny2-"+renderingType+"-sharpening-"+enable_entropic_sharpening+"-";

num_elements = 0;
for i = 1:num_rows
    num_elements = num_elements + i;
end

A_data = importdata('teapot_56.txt');
B_data = importdata('duck_56.txt');
C_data = importdata('bunny_56.txt');

orig_A = zeros([grid_size grid_size grid_size]);
for idx=1:size(A_data,1)
    x = A_data(idx,:)+grid_displace;
    orig_A(x(2), x(3), x(1)) = 1;
end
A = orig_A;

orig_B = zeros([grid_size grid_size grid_size]);
for idx=1:size(B_data,1)
    x = B_data(idx,:)+grid_displace;
    orig_B(x(3), x(2), x(1)) = 1;
end
orig_B = flip(orig_B,3);
B = orig_B;

orig_C = zeros([grid_size grid_size grid_size]);
for idx=1:size(C_data,1)
    x = C_data(idx,:)+grid_displace;
    orig_C(x(3), x(1), x(2)) = 1;
end
orig_C = flip(orig_C,3);
C = orig_C;

A = filt3(A);
A = A ./ sum(A(:));
B = filt3(B);
B = B ./ sum(B(:));
C = filt3(C);
C = C ./ sum(C(:));

width = element_height * num_elements;
height = width;

fig = figure('Position', [0 0 width height], 'Color', bgcolor);
set(gca,'visible','off')
hold on
rows = linspace(5,1,num_rows);

fig.WindowState = 'maximized';

for i=rows
    
    weight_towards_3rd = (length(rows)-i)/(length(rows)-1);
    alphas = linspace(0,1,i);
    
    for j = 1:i
        alpha = [alphas(j) 1.0-alphas(j)] .* (1.0-weight_towards_3rd);
        alpha = [alpha weight_towards_3rd];
        
        if (renderingType == "VolumeRendering")
            color = [linspace(alpha(1),alpha(1),256)' linspace(alpha(2),alpha(2),256)' linspace(alpha(3),alpha(3),256)'];
        end
        
        barycenter = WassersteinBarycenterGeneralized({A, B, C}, alpha, entropic_factor);
        img = barycenter;% ./ sum(barycenter(:));
        
        disp(['Sum = ', num2str(sum(img(:)))])
        
        w = 1/num_rows;
        h = 1/num_rows;
        p = uipanel(fig, 'Position', [w*(j-1)+(num_rows-i)/(num_rows*2), (num_rows-i)/num_rows, w, h], 'BorderType', 'none');

        volshow(img,...
            'Parent', p,...
            'Renderer', renderingType,...
            'Isovalue', isosurface_value,...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 110 0],...
            'CameraPosition',[1.3 1.1 1.3],...
            'BackgroundColor', bgcolor,...
            'Alphamap', alphamap);
    end
end

fig.WindowState = 'normal';
fig.Position = [0 0 width height];


if enable_prints
    set(gcf,'InvertHardCopy','off') % preserve background color
    for type = print_types
        print("prints/"+filename+renderingType+type{1}{2}, type{1}{1}, "-r0")
    end
end


%% Show cochlea(s)
format long g
close all
clc

sigma = 1.05;
enable_prints = true;
smooth = true;
renderingType = 'VolumeRendering';
file_resolution = 64;
shift = 8;
resolution = file_resolution+2*shift;
color = [linspace(0.8,0.8,256)', linspace(0.2,0.2,256)', linspace(1,1,256)'];
bgcolor = [1 1 1];
CameraViewAngle = 55;
CameraTarget = [0 0 0];
CameraUpVector = [0 -1 0];
CameraPosition = [-0.2 0.2 -1];
CameraPosition = CameraPosition./norm(CameraPosition)*1;
width = 400;
height = 400;

shapes = {
    "05",...
    "06",...
    "08",...
    "09",...
    "10",...
    "11",...
    "12",...
    "15",...
    "16",...
    "18",...
    "19",...
    "20",...
    "21",...
    "22",...
    "23",...
    "24",...
    "5876",...
    "6317"...
};


for shape = shapes
    
    print_position = [0 0 0 0];
    is_printing_on = true;
    while (print_position(3) ~= width && print_position(4) ~= height && is_printing_on)
        
        is_printing_on = enable_prints;
        
        shape_number = shape{1};

        A = [];
        A = loadVoxelGridFromDistanceField("shape"+shape_number+"_pca_"+num2str(file_resolution)+".txt", resolution, shift);

        if (smooth)
            A = filt3(A);
            A = A ./ sum(A(:));
        end

        fig = figure('Position', [0 0 width height], 'Color', [1 1 1]);

        volshow(A,...
            'Renderer', renderingType,...
            'Colormap', color,...
            'BackgroundColor', bgcolor,...
            'CameraViewAngle',CameraViewAngle,...
            'CameraTarget',CameraTarget,...
            'CameraUpVector',CameraUpVector,...
            'CameraPosition',CameraPosition);

        if enable_prints
%             set(gcf,'InvertHardCopy','off') % preserve background color
            set(fig,'PaperPositionMode','auto')
            set(fig, 'Position', [0 0 width height])
            set(fig, 'visible', 'off')
            for type = print_types
                if (type{1}{2} ~= '.eps')
                    set(fig, 'Position', [0 0 width height])
                    filename = "cochlea-"+num2str(file_resolution);
                    if smooth
                        filename = filename + "-smooth";
                    else
                        filename = filename + "-voxelized";
                    end
                    filename = filename +"-shape_"+num2str(shape_number);
                    print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
                    print_position = get(fig, 'Position');
                end
            end
            close all;
        end
    end
    if enable_prints
        disp('Finished printing "'+filename+'"');
    end
end

%% Compute Wasserstein distance between 2 cochlears
format long g
close all
clc

renderingType = 'VolumeRendering';
resolution = 80;
shift = 8;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

A_orig = loadVoxelGridFromDistanceField('shape10_pca_64.txt', resolution, shift);
B_orig = loadVoxelGridFromDistanceField('shape11_pca_64.txt', resolution, shift);

A = filt3(A_orig);
A = A ./ sum(A(:));
B = filt3(B_orig);
B = B ./ sum(B(:));

figure
volshow(A,...
    'Renderer', renderingType,...
    'Colormap', color,...
    'BackgroundColor', bgcolor,...
    'CameraTarget',[0 0 0],...
    'CameraViewAngle',30,...
    'CameraUpVector',[0 1 0],...
    'CameraPosition',[2 1.5 2]);

figure
volshow(B,...
    'Renderer', renderingType,...
    'Colormap', color,...
    'BackgroundColor', bgcolor,...
    'CameraTarget',[0 0 0],...
    'CameraViewAngle',30,...
    'CameraUpVector',[0 1 0],...
    'CameraPosition',[2 1.5 2]);


% Wasserstein distance
[wd,v,w] = Sinkhorn(A,B);
ed = sqrt(sum(  (A_orig(:) - B_orig(:)).^2  ));
marg = SinkhornEvalR(v,w,ones(size(v)));
disp("WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))


%% Compute Wasserstein distance between all paris of cochlears
format long g
close all
clc

sigma = 1.05;
renderingType = 'VolumeRendering';
file_resolution = 64;
shift = 8;
resolution = file_resolution + 2*shift;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

cochlea_folder = "../data/cochlears/";
cochlea_files = [
    "shape05_pca",
    "shape06_pca",
    "shape08_pca",
    "shape09_pca",
    "shape10_pca",
    "shape11_pca",
    "shape12_pca",
    "shape15_pca",
    "shape16_pca",
    "shape18_pca",
    "shape19_pca",
    "shape20_pca",
    "shape21_pca",
    "shape22_pca",
    "shape23_pca",
    "shape24_pca",
    "shape5876_pca",
    "shape6317_pca"
];

while (length(cochlea_files) > 1)
    A_file_name = cochlea_files{1};
    cochlea_files(1) = [];
    A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
    
    for k = 1:length(cochlea_files)
        B_file_name = cochlea_files{k};
        B_orig = loadVoxelGridFromDistanceField(B_file_name+"_"+file_resolution+".txt", resolution, shift);
        
        % Preprocess
        A = filt3(A_orig);
        A = A ./ sum(A(:));
        B = filt3(B_orig);
        B = B ./ sum(B(:));

        % Wasserstein distance
        [wd,v,w] = Sinkhorn(A,B);
        ed = sqrt(sum(  (A_orig(:) - B_orig(:)).^2  ));
        marg = SinkhornEvalR(v,w,ones(size(v)));
        
        disp(['("', A_file_name, '", "', B_file_name, '")'])
        disp("WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))
        disp('------------------------------')
    end
    disp('------------------------------')
end


%% Compute dissimilarity matrix of cochlea data
format long g
close all
clc

renderingType = 'VolumeRendering';
file_resolution = 64;
shift = 8;
resolution = file_resolution + 2*shift;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

cochlea_folder = "../data/cochlears/";
cochlea_files = [
    "shape05_pca",
    "shape06_pca",
    "shape08_pca",
    "shape09_pca",
    "shape10_pca",
    "shape11_pca",
    "shape12_pca",
    "shape15_pca",
    "shape16_pca",
    "shape18_pca",
    "shape19_pca",
    "shape20_pca",
    "shape21_pca",
    "shape22_pca",
    "shape23_pca",
    "shape24_pca",
    "shape5876_pca",
    "shape6317_pca"
];

dissimilarity_matrix = [];
dissimilarity_matrix_ed = [];

for i = 1:length(cochlea_files)
    A_file_name = cochlea_files{i};
    A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
    row = [];
    row_ed = [];
    for j = 1:length(cochlea_files)
        B_file_name = cochlea_files{j};
        B_orig = loadVoxelGridFromDistanceField(B_file_name+"_"+file_resolution+".txt", resolution, shift);
        
        % Preprocess
        A = filt3(A_orig);
        A = A ./ sum(A(:));
        B = filt3(B_orig);
        B = B ./ sum(B(:));

        % Wasserstein distance
        [wd,v,w] = Sinkhorn(A,B);
        ed = sqrt(sum(  (A_orig(:) - B_orig(:)).^2  ));
        marg = SinkhornEvalR(v,w,ones(size(v)));
        
        row = [row wd];
        row_ed = [row_ed ed];
        
        disp(['("', A_file_name, '", "', B_file_name, '")'])
        disp("WD = " + wd + ", ED = " + ed + ", sum of marginals = " + sum(marg(:)))
        disp('------------------------------')
    end
    dissimilarity_matrix = [dissimilarity_matrix; row];
    dissimilarity_matrix_ed = [dissimilarity_matrix_ed; row_ed];
    disp('------------------------------')
end

save("dissimilarity_matrix_"+file_resolution+".mat", 'dissimilarity_matrix', 'dissimilarity_matrix_ed');


%% Perform MDS on the dissimilarity matrix **unmodified**
format long g
close all
clc

enable_prints = false;
file_resolution = 64;

load("dissimilarity_matrix_"+file_resolution+".mat");

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

cochlea_labels

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

% Perform MDS
Y = mdscale(D, 2);

x = Y(:,1);
y = Y(:,2);

% Plot xy scatter plot
fig1 = figure
scatter(x,y)
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

dx = 0.1;
dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, cochlea_labels)



% Perform MDS on the ED matrix
D_ed = dissimilarity_matrix_ed;

% Perform MDS
Y = mdscale(D_ed, 2);

x = Y(:,1);
y = -Y(:,2);

% Plot xy scatter plot
fig2 = figure
scatter(x,y)
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

dx = 1.2;
dy = 1.2; % displacement so the text does not overlay the data points
text(x+dx, y+dy, cochlea_labels);

if enable_prints
%     set(gcf,'InvertHardCopy','off') % preserve background color
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print(fig1, "prints/"+"cochlea-mds-WD-64-raw"+type{1}{2}, type{1}{1}, "-r400");
        print(fig2, "prints/"+"cochlea-mds-ED-64-raw"+type{1}{2}, type{1}{1}, "-r400");
    end
end

%% Perform MDS on the dissimilarity matrix **heavily manipulated**
format long g
close all
clc

enable_prints = true;

load('dissimilarity_matrix.mat');
axis_limits = [-0.4 0.8 -0.5 0.4];

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

D = D ./ max(D(:));

% Perform MDS
Y = mdscale(D, 2);

x = Y(:,1);
y = Y(:,2);

% Plot xy scatter plot
fig1 = figure;
scatter(x,y);
axis(axis_limits);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

dx = 0.015;
dy = 0.015; % displacement so the text does not overlay the data points
text(x+dx, y+dy, cochlea_labels)



% Perform MDS on the ED matrix
D_ed = dissimilarity_matrix_ed;
D_ed = D_ed ./ max(D_ed(:));

% Perform MDS
Y = mdscale(D_ed, 2);

x = Y(:,1);
y = -Y(:,2)*2;

% Plot xy scatter plot
fig2 = figure;
scatter(x,y);
axis(axis_limits);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

dx = 0.015;
dy = 0.015; % displacement so the text does not overlay the data points
text(x+dx, y+dy, cochlea_labels);



% Difference plot
X = mdscale(D, 2);
Y = mdscale(D_ed, 2);
x = X(:,1);
y = X(:,2);
u = Y(:,1)-X(:,1);
v = -Y(:,2)*2-X(:,2);

fig3 = figure;
hold on
for i=1:size(D)
%     ax = plot([X(i,1) Y(i,1)], [X(i,2) -Y(i,2)*2]);
    ax = quiver(x, y, u, v);
end
axis(axis_limits);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);


% dx = 0.015;
% dy = 0.015; % displacement so the text does not overlay the data points
% text(x+dx, y+dy, cochlea_labels);


if enable_prints
%     set(gcf,'InvertHardCopy','off') % preserve background color
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print(fig1, "prints/"+"cochlea-mds-WD"+type{1}{2}, type{1}{1}, "-r400");
        print(fig2, "prints/"+"cochlea-mds-ED"+type{1}{2}, type{1}{1}, "-r400");
        print(fig3, "prints/"+"cochlea-mds-diff"+type{1}{2}, type{1}{1}, "-r400");
    end
end


%% Wasserstein Barycenters of two cochleas
close all
clc

sigma = 1.25;
color = jet(256);
renderingType = 'VolumeRendering';
filename = ['animations/teapot-duck-transformation-',lower(renderingType),'-02.gif'];
alphas = linspace(0.0001,0.9999,5);
resolution = 80;
shift = 8;

A_data = importdata(cochlea_folder+cochlea_files(6)+"_64.txt");
B_data = importdata(cochlea_folder+cochlea_files(12)+"_64.txt");

orig_A = zeros([resolution resolution resolution]);
for idx=1:size(A_data,1)
    x = A_data(idx,:)+shift;
    orig_A(x(1), x(2), x(3)) = 1;
end
A = orig_A;

orig_B = zeros([resolution resolution resolution]);
for idx=1:size(B_data,1)
    x = B_data(idx,:)+shift;
    orig_B(x(1), x(2), x(3)) = 1;
end
B = orig_B;


A = filt3(A);
A = A ./ sum(A(:));
B = filt3(B);
B = B ./ sum(B(:));

height = 250;
width = height * length(alphas)+1;

fig = figure('Position', [0 0 width height]);
hold on
n = 0;
for alpha=alphas
    barycenter = WassersteinBarycenter(A, B, alpha);
    img = barycenter;% ./ max(barycenter(:));
    
    disp("alpha = "+num2str(alpha)+ ", sum = " + num2str(sum(barycenter(:))))
    
    width = 1/length(alphas);
    p = uipanel(fig, 'Position', [width*n, 0, width, 1], 'BorderType', 'none');
    
    volshow(img,...
        'Parent', p,...
        'Renderer', renderingType,...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 110 0],...
        'CameraPosition',[-1.7 1.5 1.7]);
    n = n + 1;
end



%% Save barycenters of toy-data to perform MDS on


close all
clc

sigma = 2.5;
enable_prints = false;
enable_entropic_sharpening = true;
entropic_factor = 0.99;
renderingType = 'VolumeRendering';
isosurface_value = 0.1;
element_height = 120;
color = jet(256);
% bgcolor = [0.3 0.75 0.93]; % blue
bgcolor = [1 1 1]; % white
alphamap = linspace(0,1,256)';
alphas = 5;
num_rows = 5;
grid_size = 64;
grid_displace = 4;
filename = "3-way-barycenter-interpolation-duck-teapot-bunny2-"+renderingType+"-sharpening-"+enable_entropic_sharpening+"-";

num_elements = 0;
for i = 1:num_rows
    num_elements = num_elements + i;
end

A_data = importdata('teapot_56.txt');
B_data = importdata('duck_56.txt');
C_data = importdata('bunny_56.txt');

orig_A = zeros([grid_size grid_size grid_size]);
for idx=1:size(A_data,1)
    x = A_data(idx,:)+grid_displace;
    orig_A(x(2), x(3), x(1)) = 1;
end
A = orig_A;

orig_B = zeros([grid_size grid_size grid_size]);
for idx=1:size(B_data,1)
    x = B_data(idx,:)+grid_displace;
    orig_B(x(3), x(2), x(1)) = 1;
end
orig_B = flip(orig_B,3);
B = orig_B;

orig_C = zeros([grid_size grid_size grid_size]);
for idx=1:size(C_data,1)
    x = C_data(idx,:)+grid_displace;
    orig_C(x(3), x(1), x(2)) = 1;
end
orig_C = flip(orig_C,3);
C = orig_C;

A = filt3(A);
A = A ./ sum(A(:));
B = filt3(B);
B = B ./ sum(B(:));
C = filt3(C);
C = C ./ sum(C(:));

width = element_height * num_elements;
height = width;

fig = figure('Position', [0 0 width height], 'Color', bgcolor);
set(gca,'visible','off')
hold on
rows = linspace(5,1,num_rows);

fig.WindowState = 'maximized';

barycenters = {};

for i=rows
    
    weight_towards_3rd = (length(rows)-i)/(length(rows)-1);
    alphas = linspace(0,1,i);
    
    for j = 1:i
        alpha = [alphas(j) 1.0-alphas(j)] .* (1.0-weight_towards_3rd);
        alpha = [alpha weight_towards_3rd];
        
        if (renderingType == "VolumeRendering")
            color = [linspace(alpha(1),alpha(1),256)' linspace(alpha(2),alpha(2),256)' linspace(alpha(3),alpha(3),256)'];
        end
        
        barycenter = WassersteinBarycenterGeneralized({A, B, C}, alpha, entropic_factor);
        img = barycenter;% ./ sum(barycenter(:));
        barycenters{end+1} = barycenter;
        
        disp(['Sum = ', num2str(sum(img(:)))])
        
        w = 1/num_rows;
        h = 1/num_rows;
        p = uipanel(fig, 'Position', [w*(j-1)+(num_rows-i)/(num_rows*2), (num_rows-i)/num_rows, w, h], 'BorderType', 'none');

        volshow(img,...
            'Parent', p,...
            'Renderer', renderingType,...
            'Isovalue', isosurface_value,...
            'Colormap', color,...
            'CameraTarget',[0 0 0],...
            'CameraViewAngle',30,...
            'CameraUpVector',[0 110 0],...
            'CameraPosition',[1.3 1.1 1.3],...
            'BackgroundColor', bgcolor,...
            'Alphamap', alphamap);
    end
end

fig.WindowState = 'normal';
fig.Position = [0 0 width height];


%% Compute dissimilarity matrix from toy-data barycenters from above
format long g
close all
clc

sigma = 1.05;

% Compute dissimilarity matrices
dissimilarity_matrix = [];

for i = 1:length(barycenters)
    row = [];
    A = barycenters{i};
%     A = filt3(A);
%     A = A ./ sum(A(:));
    for j = 1:length(barycenters)
        % Wasserstein distance
        B = barycenters{j};
%         B = filt3(B);
%         B = B ./ sum(B(:));
        [wd,v,w] = Sinkhorn(A,B);
        marg1 = SinkhornEvalR(v,w,ones(size(v)));
        marg2 = SinkhornEvalL(v,w,ones(size(v)));
        
        row = [row wd];
        
        disp("WD = " + wd + ", sum of marginals = [" + num2str(sum(marg1(:))) + ", " + num2str(sum(marg2(:)))+"]")
        disp('------------------------------')
    end
    dissimilarity_matrix = [dissimilarity_matrix; row];
%     disp('------------------------------')
    disp(num2str(round(i/length(barycenters)*100)) + "% done...");
end

toy_matrix = dissimilarity_matrix;

%% Do MDS on dissimilarity matrix from toy-data above
format long g
close all
clc

enable_prints = false;
use_nonmetric_mds = false;
renderingType = "VolumeRendering";
show_scatter_plot = false;
file_resolution = 64;
resolution = 64;
shift = 0;
% color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];
width = 1920;
height = 1080;
x_shift_factor = 0.25;
y_shift_factor = 0.38;
w = 0.1;
% h = 0.1;
h = w*width/height;
window_offset = [0.05,0.1,0.91,0.82];
CameraViewAngle = 30;
CameraTarget = [0 0 0];
CameraUpVector = [0 1 0];
CameraPosition = [1.3 1.1 1.3];
CameraPosition = CameraPosition + 0.08 * cross(CameraPosition, CameraUpVector);
CameraPosition = CameraPosition./norm(CameraPosition)*2.2;
CameraTarget = CameraTarget + 0.08 * cross(CameraPosition, CameraUpVector);

% Preprocess WD dissimilarity matrix
D = toy_matrix;
D = (D+D')./2; % Symmetricize
D = D - mean(diag(D)); % Subtract the diagonal mean (shift towards 0)

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v; % Set diagonal to 0

% Perform MDS
if (use_nonmetric_mds)
    MDS = mdscale(D, 2);
else
    MDS = cmdscale(D, 2);
end
x = MDS(:,1);
y = MDS(:,2);

% % Plot xy scatter plot
% 
% fig1 = figure;
% scatter(x,y);
% axis();
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

% dx = 0.015;
% dy = 0.015; % displacement so the text does not overlay the data points
% text(x+dx, y+dy, cochlea_labels)


% Plot volume rendering

x_shift = x_shift_factor*min(x);
y_shift = y_shift_factor*min(y);

fig = figure('Position', [0 0 width height], 'Color', bgcolor);
% grid on
set(gca,'visible','on')
set(gca,'xlim', [min(x)+x_shift max(x)])
set(gca,'ylim', [min(y)+y_shift max(y)])
hold on
fig.WindowState = 'maximized';

set(gca,'Position', window_offset);

for i = 1:length(x)
    i = length(x)+1-i;
    m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
    n = (y(i)+abs(y_shift)-min(y)) / (max(y)-min(y)-y_shift) * window_offset(4)+window_offset(2)-h/2;
    color = [linspace(m,m,256)' linspace((2-m-n)/2,(2-m-n)/2,256)' linspace(n,n,256)'];
    p = uipanel(fig, 'Position', [m, n, w, h], 'BorderType', 'none');
    volshow(barycenters{i},...
        'Parent', p,...
        'Renderer', renderingType,...
        'Isovalue', isosurface_value,...
        'Colormap', color,...
        'CameraTarget',CameraTarget,...
        'CameraViewAngle',CameraViewAngle,...
        'CameraUpVector',CameraUpVector,...
        'CameraPosition',CameraPosition,...
        'BackgroundColor', bgcolor);
    disp(num2str(round(i/length(x)*100)) + "% done...")
end
disp('------------------------')

fig.WindowState = 'normal';
fig.Position = [0 0 width height];

xlabel("First MDS coordinate")
ylabel("Second MDS Coordinate")
set(gca,'FontSize',20)

if enable_prints
    set(gcf,'InvertHardCopy','off') % preserve background color
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        if (type{1}{2} ~= '.eps')
            filename = "prints/"+"mds-toy-example-VolumeRendering";
            if (~use_nonmetric_mds)
                filename = filename+"-cmds";
            end
            set(fig, 'visible', 'off');
            print(fig, filename+type{1}{2}, type{1}{1}, "-r0");
            set(gca, 'visible', 'off')
            print(fig, filename+"-wallpaper"+type{1}{2}, type{1}{1}, "-r0");
            set(gca, 'visible', 'on')
            set(fig, 'visible', 'on');
        end
    end
end

% close all

%% Perform 2-d MDS on the cochlea dissimilarity matrix and plot with volume points
format long g
close all
clc

sigma = 1.05;
enable_prints = false;
renderingType = "VolumeRendering";
prevent_layering = true;
show_scatter_plot = false;
file_resolution = 128;
resolution = 144;
shift = 8;
bgcolor = [1 1 1];
width = 1920;
height = 1080;
x_shift_factor = 0.18;
y_shift_factor = 0.22;
w = 0.05;
% h = w*width/height;
h = 0.12
window_offset = [0.07, 0.09, 0.89, 0.83];
camera_settings = {
    {
        [0 0 0],...
        [0 -1 0],...
        [-1 -1 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [-1 0 0],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [1 -0.3 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [1 -0.5 1],...
        55,...
        0.8
    }
};

load("dissimilarity_matrix_"+file_resolution+".mat");

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

% D = D ./ max(D(:));

% Perform MDS
Y = mdscale(D, 2);

x = Y(:,1);
y = Y(:,2);

if (show_scatter_plot)
    % Plot xy scatter plot
    fig1 = figure;
    scatter(x,y);
    axis fill;
    % set(gca, 'XTick', []);
    % set(gca, 'YTick', []);

    dx = 0.015;
    dy = 0.015; % displacement so the text does not overlay the data points
    text(x+dx, y+dy, cochlea_labels)
end

% Plot volume plot

x_shift = x_shift_factor*min(x);
y_shift = y_shift_factor*min(y);

for camera_idx = 1:length(camera_settings)
    camera = camera_settings{camera_idx};
    camera{3} = camera{5}*camera{3}/norm(camera{3});
    
    fig2 = figure('Position', [0 0 width height], 'Color', bgcolor);
    set(gca,'FontSize',20)
    xlabel("First MDS coordinate")
    ylabel("Second MDS coordinate")
    set(gca,'visible','on')
    set(gca,'xlim', [min(x)+x_shift max(x)])
    set(gca,'ylim', [min(y)+y_shift max(y)])
    hold on
    fig2.WindowState = 'maximized';

    set(gca,'Position', window_offset);
    taken_positions = [];
    for i = 1:length(x)
%         i = length(x)+1-i;
        % Load and preprocess voxel data
        A_file_name = cochlea_files{i};
        A_orig = [];
        A = [];
        A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
        A = A_orig;
        A = filt3(A_orig);
        A = A ./ sum(A(:));

        m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
        n = (y(i)+abs(y_shift)-min(y)) / (max(y)-min(y)-y_shift) * window_offset(4)+window_offset(2)-h/2;

        color = [linspace(m,m,256)' linspace((2-m-n)/2,(2-m-n)/2,256)' linspace(n,n,256)'];
        
        if length(taken_positions) < 1
            position_taken = false;
        else
            position_taken = false;
            position_taken = sum(taken_positions(...
                and(...
                    abs(taken_positions(:,1)-repmat(m,length(taken_positions(:,1)),1)) < w,...
                    abs(taken_positions(:,2)-repmat(n,length(taken_positions(:,2)),1)) < h...
                )...
            )) > 0;
        end
        
        if ~position_taken && prevent_layering
            p = uipanel(fig2, 'Position', [m, n, w, h], 'BorderType', 'none', 'ShadowColor', [0 1 0]);
            volshow(A,...
                'Parent', p,...
                'Renderer', renderingType,...
                'Colormap', color,...
                'CameraViewAngle',camera{4},...
                'CameraTarget',camera{1},...
                'CameraUpVector',camera{2},...
                'CameraPosition',camera{3},...
                'BackgroundColor', bgcolor,...
                'Alphamap', linspace(0,1,256)'...
            );
        end
        taken_positions = [taken_positions; [m, n, w, h]];
        disp(num2str(round(i/length(x)*100)) + "% done...")
    end
    disp('------------------------')

    fig2.WindowState = 'normal';
    fig2.Position = [0 0 width height];
    
    if enable_prints
        set(gcf,'InvertHardCopy','off') % preserve background color
%         set(gcf,'PaperPositionMode','auto')
        for type = print_types
            if (type{1}{2} ~= '.eps')
                set(fig2, 'visible', 'off')
                print(fig2, "prints/"+"cochlea-mds-VolumeRendering-"+file_resolution+"-camera-"+num2str(camera_idx)+type{1}{2}, type{1}{1}, '-r0');
                set(fig2, 'visible', 'on')
            end
        end
    end

end


%% Perform 1-d MDS on the cochlea dissimilarity matrix and plot with volume points
format long g
close all
clc

sigma = 1.05;
enable_prints = true;
renderingType = "VolumeRendering";
prevent_layering = false;
show_scatter_plot = false;
mds_dimension = 3;
cameras = [mds_dimension];
file_resolution = 64;
shift = 8;
resolution = file_resolution+2*shift;
bgcolor = [1 1 1];
width = 1920;
height = 280;
x_shift_factor = 0.18;
y_shift_factor = 0;
w = 0.05;
h = w*width/height*1.4;
% h = 0.12;
window_offset = [0.07, 0.32, 0.89, 0.83];
camera_settings = {
    {
        [0 0 0],...
        [0 -1 0],...
        [0 0 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [1 -0.3 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [-1 1 1],...
        55,...
        0.8
    }
};

coordinate_label = [...
    "First MDS Coordinate",...
    "Second MDS Coordinate",...
    "Third MDS Coordinate",...
    "Fourth MDS Coordinate"
];

load("dissimilarity_matrix_"+file_resolution+".mat");

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

% D = D ./ max(D(:));

% Perform MDS
Y = mdscale(D, 3);

x = Y(:,mds_dimension);
y = linspace(0,0,length(x));

if (show_scatter_plot)
    % Plot xy scatter plot
    fig1 = figure;
    scatter(x,y);
    axis fill;
    % set(gca, 'XTick', []);
    % set(gca, 'YTick', []);

    dx = 0.015;
    dy = 0.015; % displacement so the text does not overlay the data points
    text(x+dx, y+dy, cochlea_labels)
else

    % Plot volume plot

    x_shift = x_shift_factor*min(x);
    y_shift = y_shift_factor*min(y);

    for camera_idx = cameras
        camera = camera_settings{camera_idx};
        camera{3} = camera{5}*camera{3}/norm(camera{3});

        fig2 = figure('Position', [0 0 width height], 'Color', bgcolor);
        axes('Color','none','YColor','none');
        set(gca,'FontSize',20)
        xlabel(coordinate_label(mds_dimension))
        set(gca,'visible','on')
        set(gca,'xlim', [min(x)+x_shift max(x)])
        hold on
%         fig2.WindowState = 'maximized';

        set(fig2, 'visible', 'off')
        set(gca,'Position', window_offset);
        
        taken_positions = [];
        
        for i = 1:length(x)
    %         i = length(x)+1-i;
            % Load and preprocess voxel data
            A_file_name = cochlea_files{i};
            A_orig = [];
            A = [];
            A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
            A = A_orig;
            A = filt3(A_orig);
            A = A ./ sum(A(:));

            m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
            n = 0.4;

%             color = [linspace(m,m,256)' linspace((2-m-n)/2,(2-m-n)/2,256)' linspace(n,n,256)'];
            color = [linspace(0.8,0.8,256)' linspace(0.2,0.2,256)' linspace(1,1,256)'];

            if length(taken_positions) < 1 || ~prevent_layering
                position_taken = false;
            else
                position_taken = sum(taken_positions(abs(taken_positions(:,1)-repmat(m,length(taken_positions(:,1)),1)) < w)) > 0;
            end
            
            if ~position_taken
                p = uipanel(fig2, 'Position', [m, n, w, h], 'BorderType', 'none', 'ShadowColor', [0 1 0]);
                volshow(A,...
                    'Parent', p,...
                    'Renderer', renderingType,...
                    'Colormap', color,...
                    'CameraViewAngle',camera{4},...
                    'CameraTarget',camera{1},...
                    'CameraUpVector',camera{2},...
                    'CameraPosition',camera{3},...
                    'BackgroundColor', bgcolor,...
                    'Alphamap', linspace(0,1,256)'...
                );
                
                t = cochlea_labels(i);
                t = t{1}(6:end);
                lw = 0.007*length(t);
                lh = 0.06;
                dx = w/2-lw;
                dy = h;
                label = uicontrol(fig2, 'Style', 'text', 'String',t);
                label.Units = 'normalized';
                label.BackgroundColor = bgcolor;
                label.FontSize = 12;
                label.Position = [m+dx n+dy lw lh];
                
                
                
                taken_positions = [taken_positions; [m, n, w, h]];

                disp(num2str(round(i/length(x)*100)) + "% done...")
            end
        end
        disp('------------------------')

        fig2.WindowState = 'normal';
        fig2.Position = [0 0 width height];

        if enable_prints
            set(gcf,'InvertHardCopy','off') % preserve background color
    %         set(gcf,'PaperPositionMode','auto')
            for type = print_types
                if (type{1}{2} ~= '.eps')
                    set(fig2, 'visible', 'off')
                    filename = "prints/"+"cochlea-mds-VolumeRendering-1d-mds"+mds_dimension+"-"+file_resolution+"-camera-"+num2str(camera_idx);
                    if ~prevent_layering
                        filename = filename + "-layering";
                    end
                    filename = filename + type{1}{2};
                    print(fig2, filename, type{1}{1}, '-r0');
                    set(fig2, 'visible', 'on')
                end
            end
        end

    end


end



%% Perform 3d MDS on the cochlea dissimilarity matrix and plot cochleae in figure, including color coded surface/volume chart

%%%%% Todo: rewrite certain sections of the thesis using new vs figures
%%%%%       add exponential preprocessing to this section as an option
%%%%%       add new isosurface rendering as default when using exp prep
%%%%%       produce figures using exp prep

format long g
close all
clc

% General settings
enable_prints = false;
preprocessing_method = "gauss"; % posibilities: gauss, exp, exp2
prevent_layering = [true];
color_coordinate = {'vs-ratio'}; % posibilies: volume, surface, vs, vs-ratio
coordinate_permutations = {[3 1 2]};
show_scatter_plot = false;
bgcolor = [1 1 1];
cmap = parula(18);
% cmap = [linspace(0.5,1,18)' linspace(0.2,0.8,18)' linspace(1,0.2,18)'];
file_resolution = 64;
width = 1920;
height = 1080;
x_shift_factor = 0.18;
y_shift_factor = 0.22;
w = 0.04;
h = w*width/height*1.4;
window_offset = [0.07, 0.09, 0.82, 0.83];
coordinate_labels = ["First MDS Coordinate", "Second MDS Coordinate", "Third MDS Coordinate"];
cameras = [1];

% Gaussian preprocessing settings
sigma = 1.05;
renderingType = "VolumeRendering";
shift = 8;
resolution = file_resolution + 2*shift;

% Exponential Preprocessing settings
alpha = 30;
beta = 0;
isoValue = 1;
shift = 0;
resolution = file_resolution + 2*shift;
iso_cam = {
        [0 0 0],...
        [0 -1 0],...
        [0.2 -0.2 1],...
        55,...
        1};

% Camera settings
camera_settings = {
    {
        [0 0 0],...
        [0 -1 0],...
        [0 0 -1],...
        55,...
        1
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [1 0 -1],...
        55,...
        1
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [-1 1 1],...
        55,...
        1
    }
};

if strcmp(preprocessing_method, "gauss")
    load("dissimilarity_matrix_"+file_resolution+".mat");
elseif strcmp(preprocessing_method, "exp")
    load("dissimilarity_matrix_"+file_resolution+"_exp_alpha"+alpha+".mat");
elseif strcmp(preprocessing_method, "exp2")
    load("dissimilarity_matrix_"+file_resolution+"_exp2_alpha"+alpha+".mat");
end

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

% D = D ./ max(D(:));

for coordinate_permutation = coordinate_permutations
    
    coordinate_permutation = coordinate_permutation{1};
    
    % Perform MDS
    Y = mdscale(D, 3);

    x = Y(:,coordinate_permutation(1));
    y = Y(:,coordinate_permutation(2));
    z = Y(:,coordinate_permutation(3));

    if (show_scatter_plot)
        % Plot xy scatter plot
        fig1 = figure;
        scatter(x,y);
        axis fill;
        % set(gca, 'XTick', []);
        % set(gca, 'YTick', []);

        dx = 0.015;
        dy = 0.015; % displacement so the text does not overlay the data points
        text(x+dx, y+dy, cochlea_labels)
    else

        % Plot volume plot

        x_shift = x_shift_factor*min(x);
        y_shift = y_shift_factor*min(y);


        for cCoordinate = color_coordinate

            % Compute volumes
            V = [];
            S = importdata("cochleae_surfaces.txt");
            for i = 1:length(x)
                file_name = cochlea_files{i};
                vox = loadVoxelGridFromDistanceField(file_name+"_"+file_resolution+".txt", resolution, shift);
                V = [V; sum(vox(:)) * 24^3/file_resolution^3];
            end
            if strcmp(cCoordinate, 'volume')
                VS = V;
                color_label = "Volume";
            elseif strcmp(cCoordinate, 'surface')
                VS = S;
                color_label = "Surface";
            elseif strcmp(cCoordinate, 'vs')
                VS = V.*S/10000;
                color_label = "Volume \times Surface (1e4)";
            elseif strcmp(cCoordinate, 'vs-ratio')
                VS = S./V*100;
                color_label = "Volume \times Surface (1e4)";
            end


            for camera_idx = cameras
            for prevLayering = prevent_layering

                camera = camera_settings{camera_idx};
                camera{3} = camera{5}*camera{3}/norm(camera{3});

                fig2 = figure('Position', [0 0 width height], 'Color', bgcolor);
                set(gca,'FontSize',20)
                xlabel(coordinate_labels(coordinate_permutation(1)))
                ylabel(coordinate_labels(coordinate_permutation(2)))
                set(gca,'visible','on')
                set(gca,'xlim', [min(x)+x_shift max(x)])
                set(gca,'ylim', [min(y)+y_shift max(y)])

                hold on
                fig2.WindowState = 'maximized';

                colormap(cmap);
                c = colorbar;
                set(c,'FontSize',20)
                ylabel(c, color_label)
                c.Position = [0.93 0.1 0.01 0.8];
                c.Limits = [min(VS) max(VS)];
                c.Ticks = linspace(min(VS), max(VS),7);
                caxis([min(VS) max(VS)]);
                L=cellfun(@(x)sprintf('%1.0f',x),num2cell(get(c,'xtick')),'Un',0);
                set(c,'xticklabel',L)

                set(gca,'Position', window_offset);
                taken_positions = [];


                for i = 1:length(x)
            %         i = length(x)+1-i;
                    % Load and preprocess voxel data
                    A_file_name = cochlea_files{i};
                    A_orig = [];
                    A = [];
                    if strcmp(preprocessing_method, "gauss")
                        A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
                        A = filt3(A_orig);
                        A = A ./ sum(A(:));
                        isoValueScaled = 1;
                    elseif strcmp(preprocessing_method, "exp")
                        A_orig = loadVoxelGridFromDistanceFieldExp(A_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha, beta);
                        A = A_orig ./ sum(A_orig(:));
                        isoValueScaled = isoValue/sum(A_orig(:)) / (max(A(:)) - min(A(:)));
                        renderingType="Isosurface";
                    elseif strcmp(preprocessing_method, "exp2")
                        A_orig = loadVoxelGridFromDistanceFieldExp2(A_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha);
                        A = A_orig ./ sum(A_orig(:));
                        isoValueScaled = isoValue/sum(A_orig(:)) / (max(A(:)) - min(A(:)));
                        renderingType="Isosurface";
                    end
                    
                    m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
                    n = (y(i)+abs(y_shift)-min(y)) / (max(y)-min(y)-y_shift) * window_offset(4)+window_offset(2)-h/2;

        %             color = repmat(cmap(sort(z) == z(i),:),256,1); % 3rd dimension coloring
                    color = repmat(cmap(sort(VS) == VS(i),:),256,1);


                    if length(taken_positions) < 1 || ~prevLayering
                        position_taken = false;
                    else
                        position_taken = false;
                        position_taken = sum(taken_positions(...
                            and(...
                                abs(taken_positions(:,1)-repmat(m,length(taken_positions(:,1)),1)) < w,...
                                abs(taken_positions(:,2)-repmat(n,length(taken_positions(:,2)),1)) < h...
                            )...
                        )) > 0;
                    end

                    if ~position_taken
                        p = uipanel(fig2, 'Position', [m, n, w, h], 'BorderType', 'none', 'ShadowColor', [0 1 0]);
                        if strcmp(renderingType, "Isosurface")
                            camera = iso_cam;
                        end
                        volshow(A,...
                            'Parent', p,...
                            'Renderer', renderingType,...
                            'IsoValue',isoValueScaled,...
                            'IsosurfaceColor',color(1,:),...
                            'Colormap', color,...
                            'CameraViewAngle',camera{4},...
                            'CameraTarget',camera{1},...
                            'CameraUpVector',camera{2},...
                            'CameraPosition',camera{3},...
                            'BackgroundColor', bgcolor,...
                            'Alphamap', linspace(0,1,256)'...
                        );
                        
                        t = cochlea_labels(i);
                        t = t{1}(6:end);
                        t = [t, ', V=', num2str(round(V(i))), ', S=', num2str(round(S(i)))];
                        lw = 0.0042*length(t);
                        lh = 0.015;
                        dx = w/2-lw/2;
                        dy = -lh*0.3;
                        label = uicontrol(fig2, 'Style', 'text', 'String',t);
                        label.Units = 'normalized';
                        label.BackgroundColor = bgcolor;
                        label.FontSize = 10;
                        label.Position = [m+dx n+dy lw lh];
                        taken_positions = [taken_positions; [m, n, w, h]];
                    end
                    disp(num2str(round(i/length(x)*100)) + "% done...")
                end
                disp('------------------------')

                fig2.WindowState = 'normal';
                fig2.Position = [0 0 width height];

                if enable_prints
                    set(gcf,'InvertHardCopy','off') % preserve background color
            %         set(gcf,'PaperPositionMode','auto')
                    for type = print_types
                        if (type{1}{2} ~= '.eps')
                            set(fig2, 'visible', 'off')
                            filename = "prints/"+"cochlea-mds-VolumeRendering-3d-"+file_resolution;
                            if ~prevLayering
                                filename = filename + "-layering";
                            end
                            if strcmp(preprocessing_method, 'exp')
                                filename = filename + "-exp-alpha" + num2str(alpha);
                            elseif strcmp(preprocessing_method, 'exp2')
                                filename = filename + "-exp2-alpha" + num2str(alpha);
                            end
                            filename = filename + "-mdscoords-"+num2str(coordinate_permutation(1))+num2str(coordinate_permutation(2))+"-"+"-camera-"+num2str(camera_idx)+"-"+cCoordinate + type{1}{2};
                            print(fig2, filename, type{1}{1}, '-r0');
                            set(fig2, 'visible', 'on')
                        end
                    end
                    close all;
                end

            end
            end
        end

    end
end

%% Scatter plots of the 3d MDS cochlea analysis

format long g
close all
clc

enable_prints = true;
renderingType = "VolumeRendering";
coordinate_permutation = [1 2 3];
width = 1920;
height = 1080;
x_shift_factor = 0.18;
y_shift_factor = 0.22;
w = 0.04;
h = w*width/height*1.4;
window_offset = [0.07, 0.09, 0.92, 0.88];
coordinate_labels = ["First MDS Coordinate", "Second MDS Coordinate", "Third MDS Coordinate"];
cameras = [1 2 3]
camera_settings = {
    {
        [0 0 0],...
        [0 -1 0],...
        [0 0 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [1 0 -1],...
        55,...
        0.8
    },...
    {
        [0 0 0],...
        [0 -1 0],...
        [-1 1 1],...
        55,...
        0.8
    }
};


load("dissimilarity_matrix_"+file_resolution+".mat");

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v;

for i = 2:size(D, 1)
    for j = 1:i-1
        m = (D(i,j)+D(j,i)) / 2;
        D(i,j) = m;
        D(j,i) = m;
    end
end

% D = D ./ max(D(:));

% Perform MDS
Y = mdscale(D, 3);

x = Y(:,coordinate_permutation(1));
y = Y(:,coordinate_permutation(2));
z = Y(:,coordinate_permutation(3));



modified_cochlea_labels = cochlea_labels;
for i = 1:length(cochlea_labels)
    modified_cochlea_labels{i} = cochlea_labels{i}(6:end);
end

% Plot xy scatter plot



fig1 = figure('Position', [0 0 width height], 'Color', bgcolor);
grid on
set(gca,'FontSize',20)
xlabel(coordinate_labels(coordinate_permutation(1)))
ylabel(coordinate_labels(coordinate_permutation(2)))
set(gca,'visible','on')
set(gca,'xlim', [min(x)+x_shift max(x)])
set(gca,'ylim', [min(y)+y_shift max(y)])

hold on

set(gca,'Position', window_offset);




scatter(x,y);
range_x = max(x)-min(x);
range_y = max(y)-min(y);
axis([...
    min(x)-0.05*range_x...
    max(x)+0.1*range_x...
    min(y)-0.1*range_y...
    max(y)+0.1*range_y...
]);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

dx = range_x * 0.005;
dy = range_y * 0.0;
text(x+dx, y+dy, modified_cochlea_labels,'FontSize',14 )






if enable_prints
    set(gcf,'InvertHardCopy','off') % preserve background color
%         set(gcf,'PaperPositionMode','auto')
    for type = print_types
%             if (type{1}{2} ~= '.eps')
            set(fig1, 'visible', 'off')
            filename = "prints/"+"cochlea-mds-ScatterPlot-3d-"+file_resolution+"-permute-"+num2str(coordinate_permutation(1))+num2str(coordinate_permutation(2));
            filename = filename + type{1}{2};
            print(fig1, filename, type{1}{1}, '-r0');
            set(fig1, 'visible', 'on')
%             end
    end
end


%% Simple Gaussian figure for illustritative purposes

format long g
close all
clc



si = 1;

[X,Y] = meshgrid(linspace(-2,2,50),linspace(-2,2,50));
Z = 1/(2*pi*si) * exp(-(X.^2+Y.^2)/(2*si^2));
fig1 = figure
s = surf(X,Y,Z)
ylabel("y")
xlabel("x")
zlabel("G(x,y), \sigma = 1")
% s.EdgeColor = 'none';
% s.FaceAlpha = 1;
% axis([-2 2 -2 2 -2 2])
% 
filename = "prints/"+"gaussian-demonstration-surface";
print(fig1, filename+".eps", '-depsc', '-r0');
print(fig1, filename+".png", '-dpng', '-r0');






S = linspace(0.5,2,4)

x = linspace(-5,5,1000);

fig2 = figure
hold on
l = {};
for i = 1:length(S)
    y = 1/(2*pi*si) * exp(-(x.^2)/(2*S(i)^2))
    plot(x,y)
    l{i} = "\sigma = "+S(i);
    % s.EdgeColor = 'none';
    % s.FaceAlpha = 1;
    % axis([-2 2 -2 2 -2 2])

end
xlabel("x")
ylabel("G(x)")
legend(l,'Location','nw')

filename = "prints/"+"gaussian-demonstration-line";
print(fig2, filename+".eps", '-depsc', '-r0');
print(fig2, filename+".png", '-dpng', '-r0');



%% Show cochleas using signed-distance exponential pre-processing


format long g
close all
clc

alpha = 30;
beta = 0;

enable_prints = false;
renderingType = 'Isosurface';
isoColor = [0.8 0.2 1];
isoValue = exp(beta);
file_resolution = 64;
shift = 0;
resolution = file_resolution+2*shift;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];
CameraViewAngle = 55;
CameraTarget = [0 0 0];
CameraUpVector = [0 -1 0];
CameraPosition = [-0.2 -0.8 1];
CameraPosition = CameraPosition./norm(CameraPosition)*1;
width = 400;
height = 400;

shapes = {
%     "05",...
%     "06",...
%     "08",...
%     "09",...
%     "10",...
%     "11",...
    "12",...
%     "15",...
%     "16",...
%     "18",...
%     "19",...
%     "20",...
%     "21",...
%     "22",...
%     "23",...
%     "24",...
%     "5876",...
%     "6317"...
};


for shape = shapes
    
    print_position = [0 0 0 0];
    is_printing_on = true;
    while (print_position(3) ~= width && print_position(4) ~= height && is_printing_on)
        
        is_printing_on = enable_prints;
        
        shape_number = shape{1};

        A = [];
        A = loadVoxelGridFromDistanceFieldExp2("shape"+shape_number+"_pca_"+num2str(file_resolution)+"_signed_distance.txt", resolution, shift, alpha);
        sum_A = sum(A(:));
        A = A ./ sum_A;
        isoValueRatio = (isoValue / sum_A) / (max(A(:)) - min(A(:)));
        
        fig = figure('Position', [0 0 width height], 'Color', [1 1 1]);
        volshow(A,...
            'Renderer', renderingType,...
            'IsosurfaceColor', isoColor,...
            'Isovalue', isoValueRatio,...
            'BackgroundColor', bgcolor,...
            'Colormap', color,...
            'CameraViewAngle',CameraViewAngle,...
            'CameraTarget',CameraTarget,...
            'CameraUpVector',CameraUpVector,...
            'CameraPosition',CameraPosition);
        
        
        if enable_prints
%             set(gcf,'InvertHardCopy','off') % preserve background color
            set(fig,'PaperPositionMode','auto')
            set(fig, 'Position', [0 0 width height])
            set(fig, 'visible', 'off')
            for type = print_types
                if (type{1}{2} ~= '.eps')
                    set(fig, 'Position', [0 0 width height])
                    filename = "cochlea-"+num2str(file_resolution);
                    filename = filename + "-exponential";
                    filename = filename +"-shape_"+num2str(shape_number);
                    filename = filename + "-alpha-"+num2str(alpha);
                    print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
                    print_position = get(fig, 'Position');
                end
            end
            close all;
        end
    end
    if enable_prints
        disp('Finished printing "'+filename+'"');
    end
end


%% Display cochleae using the exponential preprocessing method and isosurface+patch functions (for greater shading control)

format long g
close all
clc

enable_prints = true;
alpha = 0.1;
beta = exp(-alpha);
isoColor = [0.8 0.2 1];
isoValue = exp(0);
file_resolution = 64;
shift = 0;
resolution = file_resolution+2*shift;
bgcolor = [1 1 1];
zoomFactor = 130;
CameraViewAngle = 55;
CameraTarget = [0 0 0];
CameraUpVector = [0 -1 0];
CameraPosition = [-0.2 0.2 -1];
CameraPosition = CameraPosition./norm(CameraPosition)*zoomFactor;
width = 400;
height = 400;


shapes = {
    "05",...
    "06",...
    "08",...
    "09",...
    "10",...
    "11",...
    "12",...
    "15",...
    "16",...
    "18",...
    "19",...
    "20",...
    "21",...
    "22",...
    "23",...
    "24",...
    "5876",...
    "6317"...
};


for shape = shapes
    
    print_position = [0 0 0 0];
    is_printing_on = true;
    while (print_position(3) ~= width && print_position(4) ~= height && is_printing_on)
        
        is_printing_on = enable_prints;
        
        shape_number = shape{1};
        
        x = linspace(-12,12,64);
        [X,Y,Z] = meshgrid(x,x,x);
        V = loadVoxelGridFromDistanceFieldExp2("shape"+num2str(shape_number)+"_pca_"+num2str(file_resolution)+"_signed_distance.txt", resolution, shift, alpha);
        sum_V = sum(V(:));
        V = V ./ sum_V;
        isoValueRatio = (isoValue / sum_V);
        
        fig = figure('Position', [500 500 width height], 'Color', [1 1 1]);
        xlabel("x")
        ylabel("y")
        zlabel("z")
        set(fig,'Color',bgcolor);
        set(gca,'XColor',bgcolor,'YColor',bgcolor,'ZColor',bgcolor,'TickDir','out')
        p = patch(isosurface(X,Y,Z,V,isoValueRatio));
        p.FaceColor = [0.8 0.2 1];
        p.EdgeColor = 'none';
        axis image
        axis vis3d
        enableDefaultInteractivity(gca);
        ax = gca;
        ax.Interactions = [rotateInteraction zoomInteraction];
        campos(CameraPosition)
        camtarget(CameraTarget)
        camup(CameraUpVector)
        camlight('headlight')
        lighting gouraud
        material dull


        
        if enable_prints
%             set(gcf,'InvertHardCopy','off') % preserve background color
            set(fig,'PaperPositionMode','auto')
            set(fig, 'Position', [0 0 width height])
            set(fig, 'visible', 'off')
            for type = print_types
                if (type{1}{2} ~= '.eps')
                    set(fig, 'Position', [0 0 width height])
                    filename = "cochlea-patch-"+num2str(file_resolution);
                    filename = filename + "-exponential";
                    filename = filename +"-shape_"+num2str(shape_number);
                    filename = filename + "-alpha-"+num2str(alpha);
                    print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
                    print_position = get(fig, 'Position');
                end
            end
            close all;
        end
    end
    if enable_prints
        disp('Finished printing "'+filename+'"');
    end
end



%% Compute Wasserstein distance between all cochleas using exponential preprocessing

format long g
close all
clc

alphas = [1.1 1.5 2 3 4 5 7 10 12 15 17 20 25 30 50];
method = 1;

for alpha = alphas
    
    disp('------------------------------')
    disp("alpha = "+num2str(alpha))
    disp('------------------------------')
    
    beta = 0;
    renderingType = 'VolumeRendering';
    file_resolution = 64;
    shift = 0;
    resolution = file_resolution + 2*shift;
    color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
    bgcolor = [1 1 1];

    cochlea_folder = "../data/cochlears/";
    cochlea_files = [
        "shape05_pca",
        "shape06_pca",
        "shape08_pca",
        "shape09_pca",
        "shape10_pca",
        "shape11_pca",
        "shape12_pca",
        "shape15_pca",
        "shape16_pca",
        "shape18_pca",
        "shape19_pca",
        "shape20_pca",
        "shape21_pca",
        "shape22_pca",
        "shape23_pca",
        "shape24_pca",
        "shape5876_pca",
        "shape6317_pca"
    ];

    dissimilarity_matrix = [];
    marginals = [];

    for i = 1:length(cochlea_files)
        A_file_name = cochlea_files{i};
        if method == 1
            A_orig = loadVoxelGridFromDistanceFieldExp(A_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha, beta);
        elseif method == 2
            A_orig = loadVoxelGridFromDistanceFieldExp2(A_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha);
        elseif method == 3
            A_orig = loadVoxelGridFromDistanceFieldExp3(A_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift);
        end
        A = A_orig ./ sum(A_orig(:));
        row = [];
        row_marginals = [];
        for j = 1:length(cochlea_files)
            B_file_name = cochlea_files{j};
            if method == 1
                B_orig = loadVoxelGridFromDistanceFieldExp(B_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha, beta);
            elseif method == 2
                B_orig = loadVoxelGridFromDistanceFieldExp2(B_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift, alpha);
            elseif method == 3
                B_orig = loadVoxelGridFromDistanceFieldExp3(B_file_name+"_"+file_resolution+"_signed_distance.txt", resolution, shift);
            end
            B = B_orig ./ sum(B_orig(:));

            % Wasserstein distance
            [wd,v,w] = Sinkhorn(A,B);
            marg = SinkhornEvalR(v,w,ones(size(v)));

            row = [row wd];
            row_marginals = [row_marginals sum(marg(:))];
            
            disp(['("', A_file_name(6:end-4), '", "', B_file_name(6:end-4), '")'] + ", WD = " + wd + ", sum of marginals = " + sum(marg(:)))
        end
        dissimilarity_matrix = [dissimilarity_matrix; row];
        marginals = [marginals; row_marginals];
    end


    filename = "dissimilarity_matrix_"+file_resolution;
    if method == 1
        filename = filename + "_exp";
    elseif method == 2
        filename = filename + "_exp2";
    elseif method == 3
        filename = filename + "_exp3";
    end
    filename = filename + "_alpha"+num2str(alpha)+".mat";
    save(filename, 'dissimilarity_matrix', 'marginals');

    dissimilarity_matrix
    marginals
    
    disp('------------------------------')
end


%% Plot mass loss as function of alpha
close all
clc

alphas = [1.1 1.5 2 3 4 5 7 10 12 15 17 20 25 30 50];

enable_prints = true;
file_resolution = 64;
methods = [1 2];
width = 600;
height = 600;
yscale = 'linear';

fig = figure('Position', [500 500 width height], 'Color', [1 1 1]);
hold on

l = {};
i = 0;
for method = methods
    i = i+1;
    
    if method == 1
        l{i} = 'symmetric';
    elseif method == 2
        l{i} = 'asymmetric';
    end
    
    mass_loss = [];
    filename = "dissimilarity_matrix_"+file_resolution;
    if method == 1
        filename = filename + "_exp_alpha";
    elseif method == 2
        filename = filename + "_exp2_alpha";
    elseif method == 3
        filename = filename + "_exp3_alpha";
    end
    for alpha = alphas

       load(filename+num2str(alpha)+".mat");
       mass_loss = [mass_loss length(marginals(:))-sum(marginals(:))];

    end


    plot(alphas,mass_loss,'.-');
    
end

legend(l)
axis tight;
xlabel("\alpha");
ylabel("Mass loss");
set(gca, 'YScale', yscale)
set(gca,'FontSize',16)

if enable_prints
    set(fig,'PaperPositionMode','auto')
    set(fig, 'Position', [0 0 width height])
    set(fig, 'visible', 'off')
    for type = print_types
        set(fig, 'Position', [0 0 width height])
        filename = "cochlea-mass-loss-vs-alpha";
        print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
    end
    set(fig, 'visible', 'on')
end
% if strcmp(yscale, "linear")
%     set(gca,'YLim',[0, exp(1)])
% end
% 
% if strcmp(yscale, "log")
%     set(gca,'YTick',logspace(-20,20,11))
% end




%% Visualise clamped preprocessing function
close all
clc

enable_prints = false;
width = 600;
height = 600;
alphas = [1 2 5 10 15 30 50];
asymmetric_clamping = false;
yscale = 'linear';

fig = figure('Position', [500 500 width height], 'Color', [1 1 1]);
hold on

l = {};
i = 1;
for alpha = alphas

    x = linspace(-max(alphas)-3, max(alphas)+3,1000);
    y = exp(-x);
    if asymmetric_clamping
        y(x<-1) = exp(1);
        y(x>alpha) = exp(-alpha);
    else
        y(x<-alpha) = exp(alpha);
        y(x>alpha) = exp(-alpha);
    end
    
    l{i} = ['\alpha = ', num2str(alpha)];

    plot(x,y)

    i = i +1;
end

legend(l)
axis tight
xlabel('d(x)')
ylabel('P(x) = exp(-d_\alpha(x))')
set(gca, 'YScale', yscale)
set(gca,'FontSize',16)


if strcmp(yscale, "linear")
    set(gca,'YLim',[0, exp(1)])
end

if strcmp(yscale, "log")
    set(gca,'YTick',logspace(-20,20,11))
end

if enable_prints
    set(fig,'PaperPositionMode','auto')
    set(fig, 'Position', [0 0 width height])
    set(fig, 'visible', 'off')
    for type = print_types
        set(fig, 'Position', [0 0 width height])
        filename = "exponentials-clamped";
        if asymmetric_clamping
            filename = filename + "-asymmetric";
        end
        filename = filename + "-" + yscale;
        print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
    end
    set(fig, 'visible', 'on')
end




%% Number of guaranteed iterations per alpha for clamped preprocessing
close all
clc

alphas = linspace(1,60,48)

enable_prints = true;
width = 600;
height = 600;
yscale = 'linear';
bases =  [0 1];

fig = figure('Position', [500 500 width height], 'Color', [1 1 1]);
hold on

y = [];
l = {};
row_y = []
for alpha = alphas
    i = 0;
    v = 1;
    while(isfinite(v))
        i = i+1;
        v = exp(i*alpha) / exp(-i*alpha);
    end
    row_y = [row_y i];
end
y = [y; row_y]
plot(alphas,row_y,'-')
l{1} = 'symmetric';

j = 1;
for base = bases
    j = j+1;
    l{j} = "asymmetric, inner clamp = "+num2str(exp(base));
    row_y = []
    for alpha = alphas
        i = 0;
        v = 1;
        while(isfinite(v))
            i = i+1;
            v = exp(i*base) / exp(-i*alpha);
        end
        row_y = [row_y i];
    end
    y = [y; row_y]
    plot(alphas,row_y,'-')
end
y

legend(l)
axis tight
xlabel('\alpha')
ylabel('n')
set(gca, 'YScale', yscale)
set(gca,'FontSize',16)

if enable_prints
    set(fig,'PaperPositionMode','auto')
    set(fig, 'Position', [0 0 width height])
    set(fig, 'visible', 'off')
    for type = print_types
        set(fig, 'Position', [0 0 width height])
        filename = "exponential-prep-least-runs";
        filename = filename + "-" + yscale;
        print(fig, "prints/"+filename+type{1}{2}, type{1}{1}, '-r0');
    end
    set(fig, 'visible', 'on')
end




