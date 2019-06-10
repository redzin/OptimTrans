%% Global Settings

close all
clear all
clc

addpath('scripts')
addpath('libraries/WOBJ_toolbox_Version2b')
addpath('../data/objects')
addpath('../data/images')
addpath('../data/cochleas')
format long
clc
clear all
close all


global sigma;
global filter_size;
global filter_padding_value;
global enable_entropic_sharpening;
global enable_prints;

sigma = 1.25;
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
%     "shape05_pca",
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


%% Show a cochlea
format long g
close all
clc

renderingType = 'VolumeRendering';
resolution = 80;
shift = 8;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

A = loadVoxelGridFromDistanceField("shape11_pca_64.txt", resolution, shift);

% A = filt3(A_orig);
% A = A ./ sum(A(:));


figure
volshow(A,...
    'Renderer', renderingType,...
    'Colormap', color,...
    'BackgroundColor', bgcolor,...
    'CameraTarget',[0 0 0],...
    'CameraViewAngle',30,...
    'CameraUpVector',[0 1 0],...
    'CameraPosition',[2 1.5 2]);

%% Compute Wasserstein distance between 2 cochlears
format long g
close all
clc

renderingType = 'VolumeRendering';
resolution = 80;
shift = 8;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

A_orig = loadVoxelGridFromDistanceField('shape05_pca_64.txt', resolution, shift);
B_orig = loadVoxelGridFromDistanceField('shape08_pca_64.txt', resolution, shift);

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

renderingType = 'VolumeRendering';
resolution = 80;
file_resolution = 64;
shift = 8;
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


%% Compute dissimilarity matrix
format long g
close all
clc

renderingType = 'VolumeRendering';
resolution = 80;
file_resolution = 64;
shift = 8;
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];

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

save('dissimilarity_matrix.mat', 'dissimilarity_matrix', 'dissimilarity_matrix_ed');


%% Perform MDS on the dissimilarity matrix **unmodified**
format long g
close all
clc

load('dissimilarity_matrix.mat');

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

barycenters = {}

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

% Compute dissimilarity matrices
dissimilarity_matrix = [];

for i = 1:length(barycenters)
    row = [];
    A = barycenters{i};
    for j = 1:length(barycenters)
        % Wasserstein distance
        B = barycenters{j};
        [wd,v,w] = Sinkhorn(A,B);
        marg = SinkhornEvalR(v,w,ones(size(v)));
        
        row = [row wd];
        
%         disp("WD = " + wd + ", sum of marginals = " + sum(marg(:)))
%         disp('------------------------------')
    end
    dissimilarity_matrix = [dissimilarity_matrix; row];
%     disp('------------------------------')
    disp(num2str(round(i/length(barycenters)*100)) + "% done...");
end

%% Do MDS on dissimilarity matrix from toy-data above
format long g
close all
clc

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
D = (D+D')./2; % Symmetricize
D = D - mean(diag(D)); % Subtract the diagonal mean (shift towards 0)

c = size(D, 1);
idx = 1:c+1:numel(D);
v = D(idx);
v(:) = 0;
D(idx) = v; % Set diagonal to 0

% Perform MDS
MDS = mdscale(D, 2);
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

renderingType = 'VolumeRendering';
color = [linspace(0.8,0.8,256)', linspace(0,0.2,256)', linspace(0,1,256)'];
bgcolor = [1 1 1];
element_height = 120;

% Plot volume rendering

width = element_height * length(barycenters);
height = width;

x_shift = 0.2*min(x);
y_shift = 0.2*min(y);

fig = figure('Position', [0 0 width height], 'Color', bgcolor);
% grid on
set(gca,'visible','on')
set(gca,'xlim', [min(x)+x_shift max(x)])
set(gca,'ylim', [min(y)+y_shift max(y)])
hold on
fig.WindowState = 'maximized';

window_offset = [0.1,0.1,0.8,0.8];
set(gca,'Position', window_offset);

for i = 1:length(x)
    w = 0.1;
    h = 0.1;
    m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
    n = (y(i)+abs(y_shift)-min(y)) / (max(y)-min(y)-y_shift) * window_offset(4)+window_offset(2)-h/2;
    color = [linspace(m,m,256)' linspace(0,1,256)' linspace(n,n,256)'];
    p = uipanel(fig, 'Position', [m, n, w, h], 'BorderType', 'none');
    volshow(barycenters{i},...
        'Parent', p,...
        'Renderer', renderingType,...
        'Isovalue', isosurface_value,...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',30,...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[1.3 1.1 1.3],...
        'BackgroundColor', bgcolor);
    disp(num2str(round(i/length(x)*100)) + "% done...")
end
disp('------------------------')

fig.WindowState = 'normal';
fig.Position = [0 0 width height];


if enable_prints
    set(gcf,'InvertHardCopy','off') % preserve background color
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        if (type{1}{2} ~= '.eps')
            print(fig2, "prints/"+"mds-toy-example-VolumeRendering-"+type{1}{2}, type{1}{1}, "-r0");
        end
    end
end


%% Perform MDS on the cochlea dissimilarity matrix nad plot with volume points
format long g
close all
clc

enable_prints = true;
show_scatter_plot = true;
file_resolution = 64;
resolution = 64;
shift = 0;


load('dissimilarity_matrix.mat');

% Preprocess labels for plots
cochlea_labels = cochlea_files;
for i = 1:size(cochlea_labels)
    cochlea_labels{i} = cochlea_labels{i}(1:end-4);
end

% Preprocess WD dissimilarity matrix
D = dissimilarity_matrix;
dis_m_mean = mean(diag(D));
D = D - dis_m_mean;
element_height = 120;

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
width = element_height * length(x);
height = width;

x_shift = 0.2*min(x);
y_shift = 0.2*min(y);

camera_settings = {
    {
        [0 0 0],...
        [0 1 0],...
        [0 0.75 0.75],...
        55
    },...
    {
        [0 0.2 0],...
        [0 1 0],...
        [-0.7 0.7 0.75]*1.2,...
        55
    }
};

for camera_idx = 1:length(camera_settings)
    camera = camera_settings{camera_idx}
    fig2 = figure('Position', [0 0 width height], 'Color', bgcolor);
    set(gca,'visible','on')
    set(gca,'xlim', [min(x)+x_shift max(x)])
    set(gca,'ylim', [min(y)+y_shift max(y)])
    hold on
    fig2.WindowState = 'maximized';

    window_offset = [0.1,0.1,0.8,0.8];
    set(gca,'Position', window_offset);

    for i = 2:length(x)
        % Load and preprocess voxel data
        A_file_name = cochlea_files{i};
        A_orig = loadVoxelGridFromDistanceField(A_file_name+"_"+file_resolution+".txt", resolution, shift);
        A = A_orig;
        A = filt3(A_orig);
        A = A ./ sum(A(:));

        w = 0.05;
        h = 0.1;
        m = (x(i)+abs(x_shift)-min(x)) / (max(x)-min(x)-x_shift) * window_offset(3)+window_offset(1)-w/2;
        n = (y(i)+abs(y_shift)-min(y)) / (max(y)-min(y)-y_shift) * window_offset(4)+window_offset(2)-h/2;

        color = [linspace(m,m,256)' linspace(0,1,256)' linspace(n,n,256)'];

        p = uipanel(fig2, 'Position', [m, n, w, h], 'BorderType', 'none');
        volshow(A,...
            'Parent', p,...
            'Renderer', renderingType,...
            'Isovalue', isosurface_value,...
            'Colormap', color,...
            'CameraViewAngle',camera{4},...
            'CameraTarget',camera{1},...
            'CameraUpVector',camera{2},...
            'CameraPosition',camera{3},...
            'BackgroundColor', bgcolor);
        disp(num2str(round(i/length(x)*100)) + "% done...")
    end
    disp('------------------------')

    fig2.WindowState = 'normal';
    fig2.Position = [0 0 width height];


    if enable_prints
        set(gcf,'InvertHardCopy','off') % preserve background color
        set(gcf,'PaperPositionMode','auto')
        for type = print_types
            if (type{1}{2} ~= '.eps')
                print(fig2, "prints/"+"cochlea-mds-VolumeRendering-camera-"+num2str(camera_idx)+"-"+type{1}{2}, type{1}{1}, "-r0");
            end
        end
    end

end
