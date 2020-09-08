% tracks HCT positions

%% initializing
clear all; close all;
clc;

global SR FR
SR=10000;     % Sampling frequency of audio
FR=30;        % frame rate

%% loading file names

[f_names]=filenames_Test; % loads all the file names

for i= 2 % enter the file number
    
    
    %% choosing the folder
    
    % Specify the folder where the files live.
    myFolder = uigetdir('C:\Users\localadmin\Desktop\30-08-2020\',f_names(i,1:16)) ;
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isdir(myFolder)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
        return;
    end
    
    % Get a list of all files in the folder with the desired file name pattern.
    
    % Gets and loads the video file
    filePattern = fullfile(myFolder, '*.avi'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    
    baseFileName = theFiles.name;
    fullFileName = fullfile(myFolder, baseFileName);
    obj = VideoReader(fullFileName);
    
    %% EOD
    % Gets and loads the video file
    filePattern1 = fullfile(myFolder, '*.daq'); % Change to whatever pattern you need.
    theFiles1 = dir(filePattern1);
    
    baseFileName1 = theFiles1.name;
    fullFileName1 = fullfile(myFolder, baseFileName1);
    [eod2_ttl1,time]=daqread(fullFileName1);
    
    eod=eod2_ttl1(:,2);
    ttl=eod2_ttl1(:,1);
    
    
    
    %%
    
    % Gets and loads the matfile
    load(f_names(i,1:20))
    data=eval(f_names(i,1:16));
    
    % Get the frame nos of approach trajectories
    BG=data.BG;
    S_frames=data.S_frames;
    % creating a mask
    image=read(obj,S_frames(1,1));
    [h_mask] = create_mask(image,2);
    %[r_mask] = load('r_mask.mat');
    %mask = r_mask.r_mask;
    mask=h_mask;
    
    %% initializing variables
    C_position=zeros(150,2);
    H_position=zeros(150,2);
    T_position=zeros(150,2);
    st_pt=1;
    k=1;
    e=1;
    poly_med=0;% 1: curve fitting, 2: median , none or only centroid: 0
    len_l=size(S_frames);
    C_Position=[];
    c_Position=[];
    EOD_pos=[];
    eod_Pos=[];
    
    %% things to play around with
    fr=5; % last frame you want to pick. If your trial end frame is in coloum 5, then 5
    ex=0;
    threshold=20;% play with this for better tracking
    visualization=0;% 2_centroid, 0_none, tail and center= '1'
    
    
    for jj=1:length(S_frames)%1:length(S_frames); Frame_selector(1)1:len_l(1)
        try
            % this function calculates the head, centroid and tail positions
            
            [C_position]=HCT_tracker_v2(obj,S_frames(jj,1),S_frames(jj,fr)+ex,BG(jj,1),BG(jj,2),visualization,st_pt,threshold,poly_med,mask);
            % gets the frame numbers
            ana= S_frames(jj,1:fr);
            %gets the times
            tet1=(SR*ana(1,1)/FR);
            tet2=(SR*(ana(1,fr)+ex)/FR);
            % this finds eods and camera ttls within the above time slots
            [pks_eod,locs_eod,pks_ttl,locs_ttl]=eod_ttl_data(eod(tet1:tet2),time(tet1:tet2),ttl(tet1:tet2));
            % gets the eods, calculates eod rate
            [EOD,Spike,EODR]=soundAnalysis2(-eod(tet1:tet2));
            %[H_position,C_position,T_position]=HCT_tracker_v3(obj,S_frames(jj,1),S_frames(jj,3),BG(jj,1),BG(jj,2),visualization,st_pt,threshold,poly_med,mask);
            
            close all
            %stores centroid position, and eod rates
            C_Position=[C_Position;C_position];
            eod_pos =EODR(locs_ttl);% finds the eod rate at camera ttl times. 
            EOD_pos=[EOD_pos;eod_pos];
            
            if jj==1
                c_Position=C_position;
                eod_Pos=eod_pos;
                
            else
                aa=size(c_Position);
                c_Position(1:length(C_position),aa(2)+1)=C_position(:,1);
                c_Position(1:length(C_position),aa(2)+2)=C_position(:,2);
                eod_Pos(1:length(eod_pos),jj)=eod_pos;
                
            end
            clear H_position C_position T_position eod_pos
            
        catch
            warning('Problem using function');
            ex_1(e)=jj
            e=e+1;
        end
        %.............
        k=k+1
        
    end
    
    % storing data
    data.eod_Pos=eod_Pos;
    data.EODR=EODR;
    data.EOD_pos=EOD_pos;
    data.C_Position=C_Position;
    data.c_Position=c_Position;
    data.mask=mask;
    eval([f_names(i,1:16),'=data;'])
    save(f_names(i,1:20),f_names(i,1:16))
    
    
end


