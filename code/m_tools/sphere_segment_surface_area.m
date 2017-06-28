function [area] = sphere_segment_surface_area(radius,zprime)

% Compute the surface area of a sphere segment as a function of radius.
% zprime    cutoff radius. For 0<radius<zprime, surface area is that of
%           sphere, for larger radii, area grows linear with radius
%
% menandrin@gmail.com, 150102

area = zeros(size(radius)); 

lgc        = radius<zprime;
area( lgc) = 2*pi*radius( lgc).^2;
area(~lgc) = 2*pi*radius(~lgc)*zprime;