function [C_position]=HCT_tracker_v2(obj,start_tr,end_tr,bg_s,bg_e,visualization,st_pt,threshold,poly_med,mask)
% start and end frames
start_tr;
end_tr;

% scaler
scaler= 3.35;
% fish start point
% ct_s=fish_start_pt;
st_pt;
% background image and threshold
bg_s; % background start image
bg_e;  % background end image
manual_threshold=threshold; % threshold after image substraction

% Usage of mask
use_mask=2; % if you want to use mask, if yes then '1'(for home compartment and circular mask), '2'(only for rectangular mask) else '0',

% polyval or median for finding the midline of the fish
blob_med=poly_med; % 1: curve fitting, 2: median , none or only centroid: 0

% visualization each iteration
visu=visualization; % if want to visualize the Head, tail and center= '1', for just centroid ='2', none='0'

% calculate angle or not
cal_angle=0;

%% creating mask to remove the home compartment: need to hand draw
%
if use_mask==1
    im_mask=read(obj,start_tr);
    figure
    imagesc(im_mask)
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    message = sprintf('Left click and hold to begin drawing around the home compartment.\nSimply lift the mouse button to finish');
    uiwait(msgbox(message));
    hFH = imfreehand();
    mask = hFH.createMask();
    load('h_mask.mat')
elseif use_mask==2
    % load circular mask
    h_mask=mask;
end


%% initializing starting position of fish: finding the starting centroid point on the fish
if st_pt==1
    ct_im=read(obj,end_tr);%start_tr
    figure
    imagesc(ct_im)
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    title('Click on the approximate center of the fish')
    ct_s=ginput(1);
    ct_roid(1,:)=ct_s;
else
    ct_roid(1,:)=ct_s;
end
%% Background calculation
Bg_I=read(obj,[bg_s bg_e]);
Bg_I=double(mean(Bg_I,3));
Bg_I=mean(Bg_I,4);

%%
%start_tr+1:end_tr
i=2;
for tt=end_tr-1:-1:start_tr
    
    % Reading and loading in current frame
    Cur_I=read(obj,tt);
    Cur_I=double(mean(Cur_I,3));
    
    % Background substraction
    Use_I=Cur_I-Bg_I;
    Use_I=abs(Use_I);
    
    % Applying mask
    if use_mask==1
        Use_I(mask==1)=1;% applying mask to remove home compartment
        Use_I(h_mask==0)=0;% applying circular mask
    elseif use_mask==2
        Use_I(h_mask==0)=0;% applying circular mask
    end
    %%
    
    % Thresholding and image filling
    
    J3 = Use_I > manual_threshold;
    BW = imfill(J3,'holes');
    
    %%
    % finding centroid
    s=regionprops(BW,Use_I,'WeightedCentroid'); %Use_I
    centroids = cat(1, s.WeightedCentroid);
    [r,c]=size(centroids);
    
    % finding Area
    s1=regionprops(BW,'Area');
    Ar = cat(1, s1.Area);
    
    % finding index number of blob with biggest area
    idx=find(Ar==max(Ar));
    
    % if it is more than one, in cases where there are two blobs with same
    % area
    blob_nos=length(idx);
    
    % angle
    
    s2=regionprops(BW,'Orientation');
    Or = cat(1, s2.Orientation);
    
    %%  finding appropriate fish position
    
    
    % making sure the blob in the current frame is not too far from the
    % previous frame.
    if blob_nos>0 % if no blob, then assigns the previous position (usually happens when the fish is in the shelter)
        if blob_nos<2 % if only one blob, then stores the value, and also calculates the distance.
            ct_roid(i,:)=centroids(idx,:); % storing centroid position
            orient(i,:)=Or(idx,:); % orientation
            distCent=sqrt((ct_roid(i-1,1)-centroids(idx,1))^2+(ct_roid(i-1,2)-centroids(idx,2))^2);
            dist(i)=distCent;
        else % more than one blob with same area, finds the nearest blob from the previous frame.
            for id=1:blob_nos
                distCent1(id)=sqrt((ct_roid(i-1,1)-centroids(idx(id),1))^2+(ct_roid(i-1,2)-centroids(idx(id),2))^2);
            end
            min_dist=find(distCent1==min(distCent1));
            ct_roid(i,:)=centroids(idx(min_dist(1)),:); % storing centroid position
            orient(i,:)=Or(idx(min_dist(1)),:);  % orientation
            dist(i)=distCent1(min_dist(1));
            clear distCent1
        end
    else
        ct_roid(i,:)=ct_roid(i-1,:);
        orient(i,:)=orient(i-1,:);
        dist(i)=dist(i-1);
    end
    
    % if the current distance is too far, then assign the previous position
    if dist(i)>50
        ct_roid(i,:)=ct_roid(i-1,:);
        orient(i,:)=orient(i-1,:);
    else
    end
    
    
    %% fitting a line to the fish related pixels
    if blob_med==1
        
        
        s7=regionprops(BW,'PixelList');
        Im7 = cat(1, s7(idx,1).PixelList);
        
        %%
        
        if abs(orient(i))<70 % this condition flips the x and y axis, based on the orientation, so as to overcome the limitaiton of the polyfit and polyval function.
            p=polyfit(Im7(:,1),Im7(:,2),3);
            y1=polyval(p,Im7(:,1));
        else
            p=polyfit(Im7(:,2),Im7(:,1),3);
            y1=polyval(p,Im7(:,2));
        end
        
        
        %%
        % finding the distance between centroid and first and last points on
        % the fitted line.
        
        if blob_nos>0
            if abs(orient(i))<70
                a_dist=sqrt((ct_roid(i,1)-Im7(1,1))^2+(ct_roid(i,2)-y1(1))^2);
                b_dist=sqrt((ct_roid(i,1)-Im7(end,1))^2+(ct_roid(i,2)-y1(end))^2);
            else
                for ik=1:length(y1)
                    all_dist(ik)=sqrt((ct_roid(i,1)-Im7(ik,2))^2+(ct_roid(i,2)-y1(ik))^2);
                end
                
                T_ind=find(all_dist==max(all_dist));
                H_ind=find(all_dist==min(all_dist));
                
            end
            % assingning head and tail position: head closer from centroid than
            % tail
            if abs(orient(i))<70
                if a_dist<b_dist
                    H_position(i,:)=[Im7(1,1),y1(1)];
                    T_position(i,:)=[Im7(end,1),y1(end)];
                elseif b_dist<a_dist
                    H_position(i,:)=[Im7(end,1),y1(end)];
                    T_position(i,:)=[Im7(1,1),y1(1)];
                end
            else
                
                
                H_position(i,:)=[y1(H_ind(1)),Im7(H_ind(1),2)];
                T_position(i,:)=[y1(T_ind(1)),Im7(T_ind(1),2)];

            end
            clear H_ind T_ind all_dist
            
        end
        if abs(orient(i))<70
            a(i)=a_dist;
            b(i)=b_dist;
        end
        % fede update
        posicionH(i,:)=H_position(i,:)/scaler;
        posicionT(i,:)=T_position(i,:)/scaler;
        pcdata(i,:)=polyfit(Im7(:,1)/scaler,Im7(:,2)/scaler,3);
        %% finding median or a line going through the fish
        
    elseif blob_med==2
        
        if blob_nos>0
            % getting the pixel cordinates
            s7=regionprops(BW,'PixelList');
            Im7 = cat(1, s7(idx,1).PixelList);
            
            % getting the median
            xa_s=unique(Im7(:,1));
            for ii=1:length(xa_s)
                id=find(Im7(:,1)==(xa_s(ii)));
                ya_s(ii)=median(Im7(id,2));
            end
            
            yb_s=unique(Im7(:,2));
            for ii=1:length(yb_s)
                id=find(Im7(:,2)==(yb_s(ii)));
                xb_s(ii)=median(Im7(id,1));
            end
            
            
            % assigning proper medians and interpolating so its easy to
            % find the closest point to the centroid
            if length(xa_s)>length(xb_s)
                pt = interparc(200,xa_s,ya_s,'spline');
                d_x=pt(:,1);
                d_y=pt(:,2);
            else
                pt = interparc(200,xb_s,yb_s,'spline');
                d_x=pt(:,1);
                d_y=pt(:,2);
            end
            
            
            % distance of each point in the medial line from centroild
            for id=1:length(d_x)
                meddist(id)=sqrt((ct_roid(i,1)-d_x(id))^2+(ct_roid(i,2)-d_y(id))^2);
            end
            
            % indx of closesest pixel in to centroid
            cent_ind=find(meddist==min(meddist));
            cent_indx=cent_ind(1);
            
            % finding arc lengths
            [arclen_a,seglen_a]=arclength(d_x(1:cent_indx),d_y(1:cent_indx)); % beginning to centroid
            [arclen_b,seglen_b]=arclength(d_x(cent_indx:end),d_y(cent_indx:end)); % centroid to end
            
            % Assigning Head and tail positions, depending which arc length
            % is shorter or longer.
            if arclen_a<arclen_b
                H_position(i,:)=[d_x(1),d_y(1)];
                T_position(i,:)=[d_x(end),d_y(end)];
            elseif arclen_b<arclen_a
                H_position(i,:)=[d_x(end),d_y(end)];
                T_position(i,:)=[d_x(1),d_y(1)];
            elseif arclen_a==arclen_b
                H_position(i,:)=[d_x(end),d_y(end)];
                T_position(i,:)=[d_x(1),d_y(1)];
            end
            
            ct_roid1(i,:)=[d_x(cent_indx),d_y(cent_indx)];
            
            % taking the median line of the fish and fitting it with a
            % polynomial
            
            p1=polyfit(pt(:,1),pt(:,2),3);
            y2=polyval(p1,pt(:,1));
            % fede variables
            posicionH(i,:)=H_position(i,:)/scaler;
            posicionT(i,:)=T_position(i,:)/scaler;
            pcdata(i,:)=polyfit(pt(:,1)/scaler,pt(:,2)/scaler,3);
        else
            % if no blobs
            ct_roid1(i,:)=ct_roid1(i-1,:);
            H_position(i,:)=H_position(i-1,:);
            T_position(i,:)=T_position(i-1,:);
        end
    end
    %% visualization
    
    if visu==1
        % for blob image
        imshow(BW)
        hold on
        plot(H_position(i,1),H_position(i,2),'ro','MarkerFaceColor','r')
        plot(T_position(i,1),T_position(i,2),'bo','MarkerFaceColor','b')
        
        if blob_med==1
            plot(ct_roid(i,1),ct_roid(i,2),'go','MarkerFaceColor','g')
            if abs(orient(i))<70
                plot(Im7(:,1),y1,'r')
            else
                plot(y1,Im7(:,2),'r')
            end
            %pause(0.2)
        elseif blob_med==2
            plot(ct_roid1(i,1),ct_roid1(i,2),'go','MarkerFaceColor','g')
            plot(pt(:,1),pt(:,2),'r')
            plot(pt(:,1),y2,'ro')
        end
        title('Head: Red, Centroid: Green, red-fitted line')
    elseif visu==2
        imshow(BW)
        hold on
        plot(ct_roid(i,1),ct_roid(i,2),'go','MarkerFaceColor','g')
    end
    
    clear xa_s ya_s xb_s yb_s d_x d_y pt
    i=i+1;
end
C_position=ct_roid;

end

