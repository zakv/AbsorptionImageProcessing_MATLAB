function [ fig ] = plot_cross_sections( image_in, DisplayName, limits, fig )
%Plots the integrated horizontal and vertical cross sections of image_in
%   === Inputs ===
%   image_in should be a 2D array.  This function is useful for both pixel
%   counts and optical depths.
%
%   DisplayName should be a string giving the name of this cross section
%   for the plot legend.
%
%   limits should be [row_min,row_max,col_min,col_max] where each refers
%   to the limis of the horizontal or vertical integration directions.
%   For example, the vertical cross section will be the sum of the rows
%   where only indices row_min:row_max are included in the integration.
%
%   fig (optional) is the figure onto which this plot should be added.  If
%   desired the data can be plotted on a figure previously output by this
%   function.  This allows the user to plot multiple cross sections from
%   different images on the same figure.  It will be assumed that the
%   figure has the same subplot structure as this function outputs.  If fig
%   is not provided, a new figure will be created and returned.
%
%   === Outputs ===
%   fig is the figure on which the data is plotted.  It may be used as an
%   input in subsequent calls to this function, which allows the user to
%   overlay cross sections from different images onto the same figure.
%
%   === Example Usage ===
%   >> %Get the image we'd like to analyze
%   >> filename = fullfile('20170405','Cool100d100d80PGCZ4.4_1_raw.ascii');
%   >> image_in = load_image(filename);
%   >> 
%   >> %Select a background region
%   >> row_min=40; row_max=60; col_min=50; col_max=80;
%   >> back_region = make_back_region(image_in,row_min,row_max,col_min,col_max);
%   >> 
%   >> %Make a basis for eigenfaces and svd algorithms
%   >> max_vectors = 20; %20 is typically a good number for this
%   >> ls_pattern = fullfile('20170405','*_back.ascii');
%   >> file_list = get_file_list(ls_pattern);
%   >> [basis_eig, mean_back] = make_basis_eig(file_list,back_region,max_vectors);
%   >> basis_svd = make_basis_svd(file_list,back_region,max_vectors);
%   >> 
%   >> %Use the basis to get the atomic cloud's optical depth
%   >> OD_eig = get_OD_eig(image_in,basis_eig,mean_back,back_region);
%   >> OD_svd = get_OD_svd(image_in,basis_svd,back_region);
%   >> 
%   >> %Plot the results
%   >> limits = [row_min,row_max,col_min,col_max];
%   >> %Generate new figure and plot results from using eigenfaces
%   >> fig = plot_cross_sections(OD_eig,'Eigenfaces Result',limits);
%   >> %Now overlay the results from using svd for comparison
%   >> plot_cross_sections(OD_svd,'SVD Result',limits,fig);
%   >> %They're very similar, let's look at the difference in a new figure
%   >> plot_cross_sections(OD_eig-OD_svd,'Difference (Eig-SVD)',limits);

%open a new figure if necessary
if nargin<4
    fig=figure();
end

%Make the correct figure active
figure(fig);

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