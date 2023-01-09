%% 3D reconstruction in computerized tomography
% SABIR Ilyass.

clear; close all;
%% 1 Firt analysis of the data
% Question 1:

% show some images from different angles:
im0 = double(imread(sprintf("data2/proj_%d.png", 0)));
im45 = double(imread(sprintf("data2/proj_%d.png", 45)));
im90 = double(imread(sprintf("data2/proj_%d.png", 90)));

figure;
subplot(1,3,1); imshow(im0, []); title("original image, angle = 0");
subplot(1,3,2); imshow(im45, []); title("original image, angle = 45");
subplot(1,3,3); imshow(im90, []); title("original image, angle = 90");


% Question 2:
f = 0.5;
line = 128;

sino = sinogram(floor(line * f), f, 0);

figure;
imshow(sino,[]); title(sprintf('Sinogram of the object on line = %d', line)); ylabel('Displacement of projection'); xlabel('Angle of projection');

%% 2 Modeling the image formation
% Question 3:
I = zeros(100,100);

radius1 = 10;
center1 = [50, 25];

for i = center1(1) - radius1: center1(1) + radius1
    for j = center1(2) - radius1: center1(2) + radius1
        if ((i - center1(1)) ^2 + (j - center1(2)) ^2 < radius1 ^2)
            I(i, j) = 1;
        end
    end
end

radius2 = 15;
center2 = [20, 75];

for i = center2(1) - radius2: center2(1) + radius2
    for j = center2(2) - radius2: center2(2) + radius2
        if ((i - center2(1)) ^2 + (j - center2(2)) ^2 < radius2 ^2)
            I(i, j) = 1;
        end
    end
end

radius3 = 5;
center3 = [90, 50];

for i = center3(1) - radius3: center3(1) + radius3
    for j = center3(2) - radius3: center3(2) + radius3
        if ((i - center3(1)) ^2 + (j - center3(2)) ^2 < radius3 ^2)
            I(i, j) = 1;
        end
    end
end


figure;
subplot(1,2,1); imshow(I,[]); title("image with 3 disks");

[R,xp] = radon(I);
subplot(1,2,2); imshow(R,[]); title("sinogram of the image");
size(R)

% Question 4:
figure;
for tetha = 0:3
    x = 1 + 16 * tetha;
    y = 1 + 16 * tetha;
    I = zeros(50,50);
    I(x, y) = 1;
    subplot(4,2,2 * tetha + 1); imshow(I,[]); title(sprintf("image%d", tetha + 1));

    [R,xp] = radon(I);
    subplot(4,2,2 * tetha + 2); imshow(R,[]); title(sprintf("sinogram of the image %d",tetha + 1));
end

%% 3 First attempt at reconstruction
% Question 6:
sinogram_line = sinogram(175, 1, 0);
plan_175 = iradon(sinogram_line, 0:179);
figure;
imshow(plan_175, []);

%{
[w, h] = size(double(imread(sprintf("data2/proj_%d.png", 0))));

object3D = zeros(180, 180, w);

for line = 1:w
    sinogram_line = sinogram(line, 1);
    plan = iradon(sinogram_line, 0:179);
    object3D(:, :, line) = plan(:, :);
end

figure;

griddata3()
%}   


% Question 7:
sinogram_radon = radon(plan_175);

disp(size(sinogram_radon));
disp(size(sinogram_line)); 

%% 4 System self-calibration
% Question 8:

% image d'angle = 0 et image d'angle = 180
image0 = imread("data2/proj_0.png");
image180 = imread("data2/proj_180.png");

figure;
subplot(1,3,1); imshow(image0,[]), title("image of the object, angle = 0");
subplot(1,3,2); imshow(image180,[]); title("image of the object, angle = 180");

% determine l'image symétrique de image180
[n, m] = size(image0);
imSymetric = zeros(n, m);
for i = 1:n
    for j = 1:m
        imSymetric(i,j) = image180(i, m + 1 - j);
    end
end
subplot(1,3,3); imshow(imSymetric,[]); title("symetric of the image of the object, angle = 180");

% Pattern
[~, h, ~] = size(image0);
x = floor(h / 2);

motif = image0(:, x - 30:30 + x);
maximum = -inf;
index = 0;
for j = 31: h - 30
    im = imSymetric(:, j - 30:j + 30);
    corr = corr2(im, motif);
    if (corr > maximum)
        index = j;
        maximum = corr;
    end
end

decalage = x - index;
%{
% Sinogram for which the axis of rotation is at the center of the image.

Sinogram3D = zeros(m - abs(decalage), 180, n);
for line = 1:n
    Sinogram3D(:, :, line) = sinogram(line, 1, decalage);
end

figure;
imshow(Sinogram3D(:, :, 128), []); title('Sinogram of the object  for which the axis of rotation is at the center on line = 128');

%% 5 Reconstruction by the filtered back-projections algorithm
% Question 9:

object3D = zeros(178, 178, n);

for line = 1:n
    plan = iradon(Sinogram3D(:, :, line), 0:179);
    object3D(:, :, line) = plan(:, :);
end


figure;
isosurface(object3D);

%% 6 Understanding of the principle of the filtered back-projection algorithm
% Question 10:
figure;
for tetha = 0:3
    x = 1 + 16 * tetha;
    y = 1 + 16 * tetha;
    I = zeros(50,50);
    I(x, y) = 1;
    subplot(4,2,2 * tetha + 1); imshow(I,[]); title(sprintf("image%d", tetha + 1));

    [R,xp] = radon(I);
    iradon_I = iradon(R, 0:179, "linear", "Ram-Lak");

    subplot(4,2,2 * tetha + 2); imshow(iradon_I,[]); title(sprintf("iradon of the sinogram of the image %d",tetha + 1));
end
%}
% Question 11:

I = sinogram(175, 1, decalage);
[n, m] = size(I);
line = floor(n / 2);

iw = 2 * floor(n / (2 * sqrt(2)));
iw_n = iw / 2;

% Filtered Back Projection, Compute the Ramlak filter
g = [0:line, line - 1:-1:1];

H = 2 * g / n;

% Compute Fourier transformation of sinogram
gf = fft(I, [], 1);

% Multiply the filter in frequency domain
gff = bsxfun(@times, gf, H');

% Do inverse Fourier transformation
gffi = real(ifft(gff, [], 1));

% Initialize the reconstructed image
img = zeros(iw);

% Compute some arguments for back projection
% Positions map of reconstructed image
[posX, posY] = meshgrid((1:iw) - iw_n);

for t = 1:180

    % Calculate the position in sinogram
    pos = posX * cosd(t) + posY * sind(t) + line;
    % Accumulate projection of each degree sinogram
    img = img + interp1(1:n, gffi(:, t), pos);
    
end
N = 180;
% Multiply the factor
img = img * (pi / (2 * N));

figure;
imshow(img, []); title("reconstruct the plan 175");




