function x = temporal_finite_difference_tenc(im)
    
    % reorder so all frames are consequitive
    x = reshape(im(:,:,:,:,end:-1:1), size(im,1), size(im,2), size(im,3), []);
    
    %finite difference + reorder back
    if size(x,4) > 1
        temp = -1*circshift(x,-1,4) + 2*x - 1*circshift(x,1,4);
        x = cat(4, x(:,:,:,1,:) - x(:,:,:,2,:), temp(:,:,:,2:end-1,:), x(:,:,:,end,:)-x(:,:,:,end-1,:));
        x = reshape(x,size(im));
        x = x(:,:,:,:,end:-1:1);
    else
        x = 0;

    end
 
    
    
end