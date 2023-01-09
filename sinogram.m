function sino = sinogram(line, f, decalage)

im0 = imresize(double(imread(sprintf("data2/proj_%d.png", 0))), f);

[h, c] = size(im0);

sino = zeros(c - abs(decalage), 180);
decal = floor(decalage * f);

if decal > 0
    for tetha = 0:179
        im = imresize(double(imread(sprintf("data2/proj_%d.png", tetha))), f);
        im = im(:, decal + 1 : c);
        sino(:, tetha + 1) = im(line, :);
    end
else
    for tetha = 0:179
        im = imresize(double(imread(sprintf("data2/proj_%d.png", tetha))), f);
        im = im(:, 1 : decal + c);
        sino(:, tetha + 1) = im(line, :);
    end
end

