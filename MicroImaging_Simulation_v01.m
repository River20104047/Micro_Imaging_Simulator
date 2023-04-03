%% This is used to simulate the processing of micro-imaging based on non-ideal conditions

% 2023/02/21 v01 Created by Zijiang Yang

%% Prepare workspace
clc, clear, close all
tic

% rng(777);



%% Input data
% Show lattice on overlaid plot
mp.lattice        = 1; % if 0, then no lattice displayed; if 1, then displace lattice
mp.plotting       = 1; % if 0, then no plot displayed

% mp = measurement point, or micro-plasitic
mp.colorindex_r   = 255 - [255, 124, 128; 255, 180, 128; 248, 243, 141; 66, 214, 164;8, 202, 209; 89, 173, 246; 157, 148, 255; 199, 128, 232];
mp.colorindex     = [255, 124, 128; 255, 180, 128; 248, 243, 141; 66, 214, 164;8, 202, 209; 89, 173, 246; 157, 148, 255; 199, 128, 232];

% Measurement points and lattice
mp.color = '#91F086';  % color of laser spot (we used green); if red, may use this color: #FFCCCC
mp.d     = 5;          % =5 as default. radius of laser spot [μm] (estimated from https://www.jasco-global.com/principle/principles-of-raman-spectroscopy-3-raman-spectroscopy-measurements/#:~:text=Confocal%20optics,typically%20less%20than%201%20%CE%BCm. for x20 objective)
mp.nw    = 100;        % number of measurement points along a-axis
mp.nh    = 100;        % number of measurement points along y-axis
mp.m     = 100;        % measurement interval [μm]

% Microplastic particles
mp.dp    = [100 150 300  100  100];   % diameter of particles [μm]
mp.ct    = [20   30  5    0    0];   % count of particles
mp.ty    = [1    2   3    1    1];   % type of MPs 1-8
mp.ov    = 0               ;    % allowed overlapped area of particles, maybe 0.1

% Other parameters
pr.dd    = mp.m/2          ;   % interpolation interval - using half of m   
pr.it    = 'nearest'       ;   % interpolation method 'nearest','v4','cubic'

% Derived parameters
mp.W     = (mp.nw - 1) * mp.m;        % Width of the lattice [μm]
mp.H     = (mp.nh - 1) * mp.m;        % Height of the lattice [μm]
mp.Np    = sum(mp.ct);                % total number of particles to generate
mp.Dp    = repelem(mp.dp, mp.ct)';    % diameter matrix of particles
mp.Tp    = repelem(mp.ty, mp.ct)';    % type matrix of particles
mp.uniqueType = unique(mp.ty);        % Unique elements (types)
mp.uniqueTyCt = numel(mp.uniqueType); % Number of unique elements

for i = 1:1:mp.Np
    mp.Cl(i,1:3) = mp.colorindex(mp.Tp(i),:); % color (type) matrix of particles
end

%% Generate particles (and lattice)
% gp = generated particles

% Overlayed particles


% Organize input paramaters for generating particles (may remove this)
gp.Count  = sum(mp.ct);
gp.W      = mp.W;
gp.H      = mp.H;
gp.da     = mp.Dp;
gp.rgb    = mp.Cl;

[gp.ct,f1] = Fun01_PrintParticle_v01(gp.Count, gp.W, gp.H, gp.da, gp.rgb);
% set(gca, 'xtick', [], 'ytick', []); 
if mp.plotting == 0
    close
end
% Generate measurement lattice
% (x0,y0) = (0,0), i.e., always in relative coodinates
% Coordinates for center of laser spots
mp.x = repmat((0:(mp.nw-1))*mp.m, mp.nh, 1);
mp.y = repmat((0:(mp.nh-1))*mp.m, mp.nw, 1)';
mp.xy = [mp.x(:), mp.y(:)];

if mp.lattice == 1
    % plot the circles using the rectangle function
    for i = 1:height(mp.xy)
        % calculate the radius of the circle
        r = mp.d/2;
        
        % calculate the position of the lower-left corner of the rectangle
        x_rect = mp.xy(i,1) - r;
        y_rect = mp.xy(i,2) - r;
        
        % plot the circle using the rectangle function
        rectangle('Position',[x_rect y_rect mp.d mp.d],'Curvature',[1 1], 'FaceColor',mp.color,'EdgeColor','black');
        hold on;
    end
    axis equal;
    xlim([-mp.m,mp.nw*mp.m])
    ylim([-mp.m,mp.nh*mp.m])
    xlabel('X [μm]')
    ylabel('Y [μm]')
    title(append('Measurement lattice: ',num2str(mp.nw),'×',num2str(mp.nh),'×',num2str(mp.m),' μm'))
end

%% Individual type of particles
% Generate particle matrix: gp.particle = [(x,y) d Tp (RGB)] = [centers, diameter, type, R,G,B]
gp.particle = [gp.ct, mp.Dp, mp.Tp,mp.Cl];

for i = 1:1:mp.uniqueTyCt  % for the i-th type

    ith.index = find(mp.Tp == mp.uniqueType(i));
    ith.ct    = gp.ct(ith.index,:);
    ith.Dp    = mp.Dp(ith.index,:);

    ith.jmax  = height(mp.xy);    % j refers to index of measurement points
    ith.kmax  = height(ith.ct);   % k refers to index of particles

    % Fine overlapped circles
    overlapping = zeros([ith.jmax,ith.kmax]);

    for j = 1:1:ith.jmax
        for k = 1:1:ith.kmax

            xj  = mp.xy(j,1);
            yj  = mp.xy(j,2);
            xk  = ith.ct(k,1);
            yk  = ith.ct(k,2);

            dp  = ith.Dp(k);
            distance = sqrt((xk - xj)^2 + (yk - yj)^2);

            if distance < 1/2*(dp + mp.d)
                overlapping(j,k) = 1;
            end
        end

        % total overlapped
        ith.ov(:,i) = sum(overlapping,2); % i-th overlapping data
    end
end

% Spatial interpolation
map = [1 1 1
       0.7 0.7 0.7];

for i = 1:1:mp.uniqueTyCt

    [xc, yc]  = meshgrid(min(mp.xy(:,1)):pr.dd:max(mp.xy(:,1)),min(mp.xy(:,2)):pr.dd:max(mp.xy(:,2)));
    zc(:,:,i) = griddata(mp.xy(:,1),mp.xy(:,2),ith.ov(:,i),xc,yc,pr.it);
    figure
    fig = figure;
    [M,c] = contourf(xc,yc,zc(:,:,i),1);
    axis equal
    set(c,'LineColor','none')
    colormap gray   
    xlabel('X(μm)')
    ylabel('Y(μm)')
    axis off
    axis tight 
    % title(append('Contour map of type ',num2str(i)));

    % Convert figure to binary image
    exportgraphics(fig,string(append('Contour map of type ',num2str(i),'.png')),'Resolution',300); % Export figure, it is related to pixels and scaled space
    close
    I         = imread(string(append('Contour map of type ',num2str(i),'.png')));
    BW        = im2bw(I);
    [hti,wdi] = size(BW);
    r    = 1/2*(mp.H/(mp.m + hti) + mp.W /(wdi + mp.m)); % average ratio between pixels and real size. h/hi ≈ w/wi, but there are some small differences due to white space
        


     % extract information from image
     Area = regionprops(BW,'area');
     Cent = regionprops(BW,'centroid');
     ALmj = regionprops(BW,'MajorAxisLength');
     ALmn = regionprops(BW,'MinorAxisLength');

     Area = cat(1, Area.Area);
     ALmj = cat(1, ALmj.MajorAxisLength) * r; % *100/15 = convertion factor
     ALmn = cat(1, ALmn.MinorAxisLength) * r;
     Cent = cat(1, Cent.Centroid);
        
     index = find(ALmj == max(ALmj));

     Area(index) = [];
     ALmj(index) = [];
     ALmn(index) = [];
     Cent(index,:) = [];

     ALmj0(1:length(ALmj),i) = ALmj;
     ALmn0(1:length(ALmn),i) = ALmn;    

     % Plot with identified particle and data
     if mp.plotting == 1
     imshow(BW,'Border','tight','Colormap',map)
     imshow(BW,map)
     hold on
     text(Cent(:,1),Cent(:,2),strcat('(',strcat(num2str(round(ALmj,0))),', ',num2str(round(ALmn,0)),')'),'Color','red')
     % title([string(append('Contour map of type ',num2str(i))),num2str(length(Area))])
     title(append('Type #',num2str(i),' | Major and minor axis length (μm), n = ',num2str(length(Area))))
     xlabel('X-axis')
     ylabel('Y-axis')
     axis off
     end
  
end

%% Overlaid image

if mp.plotting == 1
for i = 1:1:mp.uniqueTyCt    
    if i == 1
        A   = imread(string(append('Contour map of type ',num2str(i),'.png')));
        BW  = im2bw(imread(string(append('Contour map of type ',num2str(i),'.png'))));% Export figure, it is related to pixels and scaled space
        B   = imoverlay(A,BW,mp.colorindex(i,:)/255);         
    elseif i > 1
        A   = imread(string(append('Overlaied ',num2str(i-1),'.png')));
        BW  = im2bw(imread(string(append('Contour map of type ',num2str(i),'.png'))));% Export figure, it is related to pixels and scaled space

        % re-size image
        max_height = max(size(A, 1), size(BW, 1));
        max_width = max(size(A, 2), size(BW, 2));

        A_resized = imresize(A, [max_height, max_height]);
        BW_resized = imresize(BW, [max_height, max_height]);

        B   = imoverlay(A_resized,BW_resized,mp.colorindex(mp.uniqueType(i),:)/255);     
    end

    figOV = figure;
    imshow(B);
    exportgraphics(figOV,string(append('Overlaied ',num2str(i),'.png')),'Resolution',300); % Export figure, it is related to pixels and scaled space

    if i < mp.uniqueTyCt
        close
    end
end
end

%% Summary statistics
% Axis length
T.LengthMajorAxis = array2table(ALmj0,'VariableNames',string(mp.uniqueType));
T.LengthMinorAxis = array2table(ALmn0,'VariableNames',string(mp.uniqueType));

% nominal detection rate (norDR) and average de (for one type with multiple sizes)
for i = 1:1:mp.uniqueTyCt  
    
    % nominal dectection rate
    typeindex = mp.uniqueType(i);
    totalct   = sum(mp.Tp == typeindex);
    Ndetected = nnz(ALmj0(:,i));
    norDR(i)  = Ndetected / totalct;

    % de (equlvalent diameter)
    idx           = find(mp.Tp == mp.uniqueType(i));
    de_1type      = mp.Dp(idx);

    SM.de_mean(i) = mean(de_1type);
    SM.de_std(i)  = std(de_1type);

    SM.ct(i)      = length(de_1type);
end
T.nominalDR   = array2table(norDR,'VariableNames',string(mp.uniqueType)); 

% Advanced summary
AS.de_mean    = SM.de_mean;
AS.de_std     = SM.de_std;
AS.count      = SM.ct;
AS.norDR      = norDR;
for i = 1:1:mp.uniqueTyCt
    AS.lmj_mean(i)  = mean(ALmj0(ALmj0(:,i) ~=0,i));
    AS.lmj_std(i)   = std(ALmj0(ALmj0(:,i) ~=0,i));
    AS.lmj_Q50(i)   = median(ALmj0(ALmj0(:,i) ~=0,i));
    AS.lmj_Q025(i)  = prctile((ALmj0(ALmj0(:,i) ~=0,i)),2.5);
    AS.lmj_Q975(i)  = prctile((ALmj0(ALmj0(:,i) ~=0,i)),97.5);

    AS.lmn_mean(i)  = mean(ALmn0(ALmn0(:,i) ~=0,i));
    AS.lmn_std(i)   = std(ALmn0(ALmn0(:,i) ~=0,i));
    AS.lmn_Q50(i)   = median(ALmn0(ALmn0(:,i) ~=0,i));
    AS.lmn_Q025(i)  = prctile((ALmn0(ALmn0(:,i) ~=0,i)),2.5);
    AS.lmn_Q975(i)  = prctile((ALmn0(ALmn0(:,i) ~=0,i)),97.5);  

end

AS.summary =   [AS.de_mean
                AS.de_std
                AS.count
                AS.norDR
                AS.lmj_mean
                AS.lmj_std
                AS.lmj_Q50
                AS.lmj_Q025
                AS.lmj_Q975
                AS.lmn_mean
                AS.lmn_std
                AS.lmn_Q50
                AS.lmn_Q025
                AS.lmn_Q975];

T.Summary = array2table(AS.summary,'VariableNames',string(mp.uniqueType),'RowNames',{'mean de','std de','N','nDR', 'lmj mean','lmj std','lmj Q50','lmj Q2.5','lmj Q97.5','lmn mean','lmn std','lmn Q50','lmn Q2.5','lmn 97.5'}'); 

disp(T.Summary)

if mp.plotting == 0
    close all
end
%% May add some output for further analysis




















toc