function [ fig ] = plot_cross_sections( image_in, DisplayName, limits, fig )
%Plots the integrated horizontal and vertical cross sections of image_in
%   image_in should be a 2D array
%
%   DisplayName should be a string giving the name for the legend
%
%   limits should be [row_min,row_max,col_min,col_max] where each refers
%   to the limis of the horizontal or vertical integration directions.
%   For example, the vertical cross section will be the sum of the rows
%   where only indices row_min:row_max are included in the integration
%
%   fig is an optional argument.  If desired the data can be plotted on
%   a figure previously output by this function.  It will be assumed
%   that the figure has the same subplot structure as this function
%   outputs.  If fig is not provided, a new figure will be created and
%   returned

%open a new figure if necessary
if nargin<4
    fig=figure();
end

%Unpack limits
row_min=limits(1);
row_max=limits(2);
col_min=limits(3);
col_max=limits(4);

%Start with horizontal cross section (i.e. Integrate columns)
subplot(2,1,1);
hold on
title('Horizontal Cross Section');
plot(sum(image_in(row_min:row_max,:),1),'DisplayName',DisplayName);
legend('show');
hold off

%Vertical cross section (i.e. Integrate rows)
subplot(2,1,2);
hold on
title('Vertical Cross Section');
plot(sum(image_in(:,col_min:col_max),2),'DisplayName',DisplayName);
hold off
end