function [ MASK ] = mask_ring( r1, r2, w )
%MASK_RING Summary of this function goes here
%   Detailed explanation goes here
%r1: small radiam of the ring
%r2: large radiam of the ring
%w : size of image
%mask: 1:ring; 0:not ring

    MASK = zeros(w);
    for i = 1:w(1)
        for j = 1:w(2)
            distance_2 = ((1+w(1))/2-i)^2 + ((1+w(2))/2-j)^2;
            if distance_2 <= r2^2 & distance_2 >= r1^2
                MASK(i,j)=1;
            end
        end
    end
end

