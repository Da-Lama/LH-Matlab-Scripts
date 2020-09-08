function [h_mask] = create_mask(image,rect_points)

%rect_points = 1, for creating a polygon, rect_points =2 for creating a
%rectangle
if rect_points==1
    ct_im=image;
    figure
    imagesc(ct_im)
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    message = sprintf('Left click to create a polygon.\nSimply lift the mouse button to finish');
    uiwait(msgbox(message));
    hFH = impoly();
    h_mask = hFH.createMask();
    
else
    
    ct_im=image;
    figure
    imagesc(ct_im)
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    message = sprintf('click to create a rectangle');
    uiwait(msgbox(message));
    hFH = imrect();
    h_mask=hFH.createMask();
    
end

end