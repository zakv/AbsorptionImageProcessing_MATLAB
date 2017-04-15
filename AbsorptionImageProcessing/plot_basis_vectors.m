function plot_basis_vectors( basis, indices, image_size)
%Plots image_in with a colorbar
%   basis should be a 2D array of basis vectors from get_basis()
%
%   indices should be a 1D array of indices of the basis vectors to plot
%   (i.e. 1:3 will plot the first three basis vectors)
%
%   image_size is an optional argument.  When given should be a vector
%   giving the dimensions of the image.  If it is not given, the iamges are
%   assumed to be 121x121.
%
%   Example Usage:
%   >> plot_basis_vectors(basis, 1:10, [121,121]);

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