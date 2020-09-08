%% plotting 

eod_Pos=data.eod_Pos;
c_Position=data.c_Position;
C_Position=data.C_Position;
ct_im=read(obj,S_frames(1,1));

for i = 1: length(S_frames)
    
    ana= S_frames(i,1:3);
    fr_no(i,1)=ana(1,fr)+ex-ana(1,1);
    fr_nu(i,1)=length(find(eod_Pos(:,i)>0));
    
end

figure
imshow(ct_im)
hold on
plot(C_Position(:,1),C_Position(:,2),'go','MarkerFaceColor','g')

%% plotting 
figure
imshow(ct_im)
hold on
al=0;
for i=1:length(S_frames)
    
    if i==1
        if fr_no(i,1)>fr_nu(i,1)
            scatter(c_Position(1:fr_nu(i,1),1),c_Position(1:fr_nu(i,1),2),20,eod_Pos(1:fr_nu(i,1),1),'filled')
            position=[c_Position(1:fr_nu(i,1),1),c_Position(1:fr_nu(i,1),2)];
        elseif fr_no(i,1)<=fr_nu(i,1)
            
            scatter(c_Position(1:fr_no(i,1),1),c_Position(1:fr_no(i,1),2),20,eod_Pos(1:fr_no(i,1),1),'filled')
            position=[c_Position(1:fr_nu(i,1),1),c_Position(1:fr_nu(i,1),2)];
        end
        
        al=al+2;
    else
        
        if fr_no(i,1)>fr_nu(i,1)
            scatter(c_Position(1:fr_nu(i,1),al+1),c_Position(1:fr_nu(i,1),al+2),20,eod_Pos(1:fr_nu(i,1),1),'filled')
            
        elseif fr_no(i,1)<=fr_nu(i,1)
            
            scatter(c_Position(1:fr_no(i,1),al+1),c_Position(1:fr_no(i,1),al+2),20,eod_Pos(1:fr_no(i,1),1),'filled')
            
        end
        al=al+2;
        
        
        
    end
end
