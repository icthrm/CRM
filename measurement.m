function [ psnrval, ssimval, nccval ] = measurement( path1, path2, r, length, flag_reshape)%, peakval)
%MEASUREMENT Summary of this function goes here
%   path2 ref
currentFolder = pwd;
addpath(genpath(currentFolder));
image1= tom_mrcreadimod(path1);
image2= tom_mrcreadimod(path2);
image1= image1.Value;
image2= image2.Value;
if flag_reshape
    image1 = reshape(image1, length, length);
    image2 = reshape(image2, length, length);
end

w=size(image1)
MASK = mask_ring( 0, r, w );
image1 = image1 .* MASK;
image2 = image2 .* MASK;

%mapminmax(image1，0，1);
%mapminmax(image2，0，1);

psnrval = psnr(image1, image2, 1);%,peakval);
ssimval = ssim(image1, image2);

nccval_all = normxcorr2(image2,image1);
w2 = floor((size(nccval_all)+1)/2);
nccval = nccval_all(w2(1),w2(2));

x=[num2str(psnrval),' ', num2str(ssimval), ' ', num2str(nccval)];
disp(x);
%disp(num2str(psnrval));disp(num2str(ssimval));disp(num2str(nccval));
end
