function [ fig ] = plot_image( image_in, plot_name, clim )
%Plots image_in with a colorbar
%   image_in should be a 2D array
%   plot_name (optional) should be a string for the plot title
%   clim (optional) is a 2-element [cmin, cmax] giving the upper and lower
%       limits of the color scale
fig=figure();
image(image_in,'CDataMapping','scaled');
if nargin>1
    title(plot_name,'Interpreter','none');
end
colorbar();
if nargin>2
    set(gca,'clim',clim);
end
end