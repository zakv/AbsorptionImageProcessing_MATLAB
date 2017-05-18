function plot_basis_vectors( basis, indices, image_size)
%Plots the basis images specified by indices
%   === Inputs ===
%   basis should be a 2D array of basis vectors from get_basis_eig() or
%   get_basis_svd().
%
%   indices should be a 1D array of indices of the basis vectors to plot
%   (i.e. 1:3 will plot the first three basis vectors).
%
%   image_size (optional) should be a vector giving the dimensions of the
%   image.  If it is not given, the images are assumed to be 121x121.
%
%   === Example Usage ===
%   >> %Get an image we'd like to analyze in the future
%   >> %(we need to know what dimensions we need to make back_region)
%   >> filename = fullfile('20170405','Cool100d100d80PGCZ4.4_1_raw.ascii');
%   >> image_in = load_image(filename);
%   >> 
%   >> %Select a background region
%   >> row_min=40; row_max=60; col_min=50; col_max=80;
%   >> back_region = make_back_region(image_in,row_min,row_max,col_min,col_max);
%   >> 
%   >> %Make the basis
%   >> max_vectors = 20; %20 is typically a good number for this
%   >> ls_pattern = fullfile('20170405','*_back.ascii');
%   >> file_list = get_file_list(ls_pattern);
%   >> [basis_eig,mean_back] = make_basis_eig(file_list,back_region,max_vectors);
%   >>
%   >> %Plot the three most dominant basis images
%   >> image_size = size(image_in);
%   >> plot_basis_vectors(basis_eig,1:3,image_size);

%Set default image size
if nargin<3
    image_size=[121,121];
end

%Now plot the basis vectors
for index=indices
    figure();
    vector=basis(:,index);
    %There's probably a smarter way to do the next line, but I'm not sure
    %exactly what that way is right now
    image_in=reshape(vector,image_size);
    image(image_in,'CDataMapping','scaled');
    plot_name=sprintf('Basis Vector #%d',index);
    title(plot_name,'Interpreter','none');
    colorbar();
end
end