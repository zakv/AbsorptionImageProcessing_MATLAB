function [ fig ] = plot_image( image_in, plot_name, clim )
%Plots image_in with a colorbar
%   === Inputs ===
%   image_in should be a 2D array
%
%   plot_name (optional) should be a string for the plot title
%
%   clim (optional) is a 2-element array [cmin, cmax] giving the upper and
%   lower limits of the color scale
%
%   === Example Usage ===
%   >> filename = fullfile('20170412','Savefile_45_back.ascii');
%   >> image_out = load_image(filename);
%   >> plot_image(image_out,'Test of plot_image()');

fig=figure();
image(image_in,'CDataMapping','scaled');
if nargin>1
    %Prevent interpreter from using "_" to indicate subscripts
    title(plot_name,'Interpreter','none');
end
colorbar();
if nargin>2
    set(gca,'clim',clim);
end
end