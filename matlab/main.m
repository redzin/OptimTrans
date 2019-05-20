%% Global Settings

close all
clear all
clc

addpath('scripts')
format long
clc
clear all
close all

enable_prints = true;

global sigma;
global filter_size;
global filter_padding_value;

sigma = 5.0;
filter_size = 129;
filter_padding_value =  0.0;

sigmas = [4.1]; %[1.0 2.5 4.1 5.0];
filter_sizes = [101]; %[5 41 45 101];
print_types = {{"-depsc", ".eps"}, {"-dpng", ".png"}};


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

A(1:64,1:64,1:48) = 0;
A(12:64,1:64,1:64) = 0;
A(1:64,12:64,1:64) = 0;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(A,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',[30],...
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
        'CameraViewAngle',[30],...
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
        'CameraViewAngle',[30],...
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
            'CameraViewAngle',[30],...
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

X = rand([64,64,64]);
X(X > 0.99) = 1;
X(X < 0.99) = 0;
orig_X = X;
Y = rand([64,64,64]);
Y(Y > 0.99) = 1;
Y(Y < 0.99) = 0;
orig_Y = Y;
% B = zeros([64,64,64]);
% B(32,32,32) = 1;
% orig_B = B;

% FIGURE 1; Pre-Gaussian filter
color = jet(256);
figure
volshow(orig_X,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',[30],...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-A-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end
figure
volshow(orig_Y,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',[30],...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test3-B-pre-gaussian-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

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
        'CameraViewAngle',[30],...
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
        'CameraViewAngle',[30],...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test-B-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
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
        'CameraViewAngle',[30],...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test-A-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end

figure
volshow(margL,...
        'Renderer', 'MaximumintensityProjection',...
        'Colormap', color,...
        'CameraTarget',[0 0 0],...
        'CameraViewAngle',[30],...
        'CameraUpVector',[0 1 0],...
        'CameraPosition',[2 1.5 2]);
if enable_prints
    set(gcf,'PaperPositionMode','auto')
    for type = print_types
        print("prints/3d-test-B-marginal-"+sigma+"-sigma-"+filter_size+"-filter-size"+type{1}{2}, type{1}{1}, "-r0")
    end
end







gaussian