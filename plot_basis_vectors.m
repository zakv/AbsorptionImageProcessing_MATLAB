function plot_basis_vectors( basis, indices )
%Plots image_in with a colorbar
%   basis should be a 2D array of basis vectors from get_basis()
%   indices should be a 1D array of indices
for index=indices
    figure();
    vector=basis(:,index);
    %There's probably a smarter way to do the next line, but I'm not sure
    %exactly what that way is right now
    image_in=reshape(vector,121,121);
    image(image_in,'CDataMapping','scaled');
    plot_name=sprintf('Basis Vector #%d',index);
    title(plot_name,'Interpreter','none');
    colorbar();
end
end