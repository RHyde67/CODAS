% Copywrite R Hyde 2015
% Wrapper to test CODAS clustering
% Suggested defaults are pre-loaded into the dialog box.
% Note that due to the randomised data order cluster colours may vary
% between runs. On some occasion clusters may appear the same colour,
% usually comparing the micro-cluster plot with the data plot will highlight
% the differences.

clear
clear functions
%% Load data from file: un-rem the groups of statements to load data and suitable parameters to pass to the GUI
% DataIn=csvread('DS1.csv',1,0);
% DataIn=DataIn(:,1:2);
% Defaults = {'11', '8', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
% DataIn=csvread('DS2.csv',1,0);
% DataIn=DataIn(:,1:2);
% Defaults = {'13', '8', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
% DataIn=csvread('DS3.csv',1,0);
% DataIn=DataIn(:,1:2);
% Defaults = {'11', '4', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
DataIn=csvread('ChainLink3DNoise.csv',1,0);
DataIn=DataIn(:,1:3);
Defaults = {'0.2', '5', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
%
% DataIn=csvread('ChainLink3D.csv',1,0);
% DataIn=DataIn(:,1:3);
% Defaults = {'0.2', '3', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
% DataIn=csvread('gaussian5000.csv',1,0);
% DataIn=DataIn(:,1:2);
% Defaults = {'0.022', '8', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
% DataIn=csvread('SpiralData2_very_noisy.csv',1,0);
% DataIn=DataIn(:,1:2);
% Defaults = {'0.45', '6', '0'}; % values for Intial Radius, Min Cluster Size, Verbose
% 
% DataIn=csvread('SpiralData2_92.csv',1,0);
% DataIn=DataIn(:,1:end-1);
% Defaults = {'0.3', '2', '0'}; % values for Intial Radius, Min Cluster Size, Verbose


% DataIn=DataIn(:,1:2);

%% Plot Data
figure(1)
clf
if size(DataIn,2)==2 | size(DataIn,2)>3
    scatter(DataIn(:,1),DataIn(:,2),5,'o')
elseif size(DataIn,2)==3
    scatter3(DataIn(:,1),DataIn(:,2),DataIn(:,3),5,'o')
end
title('Raw Data')

%% Randomise data order
DataIn=DataIn(randperm(size(DataIn,1)),:);

%% Normalise DataDataIn=double(DataIn(:,1:2));
% for idx=1:size(DataIn,2)
%     DataIn(:,idx) = (DataIn(:,idx)-min(DataIn(:,idx))) / (max(DataIn(:,idx))-min(DataIn(:,idx)));
% end

%% Initialise
prompt = {'Initial radius:', 'Minimum cluster size:', 'Verbose:'};
dlg_title = 'Input parameters';
num_lines= 1;
Inputs = inputdlg(prompt,dlg_title,num_lines, Defaults);
ClusterParameters=[str2double(cell2mat(Inputs(1,1))),str2double(cell2mat(Inputs(2,1)))];
Verbose=str2double(cell2mat(Inputs(3,1)));

%% Set plot information
NumColours=20;
Colours=distinguishable_colors(NumColours);
if size(DataIn,2)==2
   AXSz=[min(DataIn(:,1)) max(DataIn(:,1)) min(DataIn(:,2)) max(DataIn(:,2))]; 
elseif size(DataIn,2)==3
   AXSz=[min(DataIn(:,1)) max(DataIn(:,1)) min(DataIn(:,2)) max(DataIn(:,2)) min(DataIn(:,3)) max(DataIn(:,3))]; 
else
    Verbose=0;
    sprintf('Too many dimensions to display')
end

%% ### Run CODAS ###
tic
for idx1 = 1 : size(DataIn,1)
    NewSample=DataIn(idx1,:);
    
    [Clusters]=CODAS2_ver02(ClusterParameters,NewSample);
    
    if Verbose==1
       figure(1)
       hold off
       if size(DataIn,2)==2
%            scatter(DataIn(1:idx1,1),DataIn(1:idx1,2),5,'o')
           for idx2=1:size(Clusters.Centre,1)
               ClrNum=rem(Clusters.global(idx2),NumColours)+1;
               circles(Clusters.Centre(idx2,1),Clusters.Centre(idx2,2),Clusters.Radius(idx2,1),...
                   'points',20,'facecolor',Colours(ClrNum,:),'facealpha',0.2);
           end
       elseif size(DataIn,2)==3
         for idx2=1:size(Clusters.Centre,1)
            CCol=rem(Clusters.global(idx2),20)+1;
            if Clusters.Count(idx2)<MinClusSize % if outlier

            else % if cluster
                [x, y, z]=sphere;
                Clusters.Radius(idx2);
                hmesh=mesh((x*Clusters.Radius(idx2))+Clusters.Centre(idx2,1),...
                    (y*Clusters.Radius(idx2))+Clusters.Centre(idx2,2),...
                    (z*Clusters.Radius(idx2))+Clusters.Centre(idx2,3));
                set(hmesh,'FaceColor',Clr(CCol,:),'FaceAlpha',0.2,'LineStyle','none');
            end
            drawnow
         end
       end
       title('Received Data with Micro-Cluster Overlay Coloured by Global Cluster')
       axis(AXSz)
       drawnow
    end
    
end   
%% ### END CODAS ###
t2=toc;
sprintf('Run time: %2f s, Time per Sample: %4f s',t2,t2/size(DataIn,1))

%% Plot Clusters
figure(2)
clf
if size(DataIn,2)==2 | size(DataIn,2)>3
%    scatter(DataIn(1:idx1,1),DataIn(1:idx1,2),5,'o')
   for idx2=1:size(Clusters.Centre,1)
       ClrNum=rem(Clusters.global(idx2),NumColours)+1;
       if Clusters.Count(idx2)>ClusterParameters(1,2)
       circles(Clusters.Centre(idx2,1),Clusters.Centre(idx2,2),Clusters.Radius(idx2,1),...
           'points',20,'facecolor',Colours(ClrNum,:),'facealpha',0.2);
       end
       title(sprintf('End of Data\n Micro-Clusters Coloured by Global Cluster\n in First 2 Dimensions'))
   end
elseif size(DataIn,2)==3
   for idx2=1:size(Clusters.Centre,1)
       
        ClrNum=rem(Clusters.global(idx2),NumColours)+1;
        if Clusters.Count(idx2)>ClusterParameters(1,2)
            [x, y, z]=sphere;
            Clusters.Radius(idx2);
            hmesh=mesh((x*Clusters.Radius(idx2))+Clusters.Centre(idx2,1),...
                (y*Clusters.Radius(idx2))+Clusters.Centre(idx2,2),...
                (z*Clusters.Radius(idx2))+Clusters.Centre(idx2,3));
            set(hmesh,'FaceColor',Colours(ClrNum,:),'FaceAlpha',0.2,'LineStyle','none');
        end
        hold on
        title(sprintf('End of Data\n Micro-Clusters Coloured by Global Cluster'))
   end
   drawnow
end

% axis(AXSz)

%% Plot Data coloured by cluster
D2C=pdist2(DataIn,Clusters.Centre); % distances to centres
[D2C,NC]=min(D2C,[],2); % min distances & nearest cluster
GC=Clusters.global(NC); % global cluster of nearest centre
GC(Clusters.Count(NC)<ClusterParameters(1,2))=0; % set small cluster to zero

% Renumber clusters for colouring
[C, ia, ic] = unique(GC);
C=[1:size(C,1)];
renum=C(ic).';
GC=renum;
PlotClrs=Colours(rem(GC,NumColours)+1,:);
if size(DataIn,2)==2 | size(DataIn,2)>3
    figure(3)
    clf
    scatter(DataIn(GC~=1,1),DataIn(GC~=1,2),5,PlotClrs(GC~=1,:))
    hold on
    scatter(DataIn(GC==1,1),DataIn(GC==1,2),5,'xk')
    title(sprintf('End of Data\n Data Placed into Cluster Results\n in First 2 Dimensions'))
elseif size(DataIn,2)==3
    figure(3)
    clf
    scatter3(DataIn(GC~=1,1),DataIn(GC~=1,2),DataIn(GC~=1,3),5,PlotClrs(GC~=1,:))
    hold on
    scatter3(DataIn(GC==1,1),DataIn(GC==1,2),DataIn(GC==1,3),5,'xk')
    title(sprintf('End of Data\n Data Placed into Cluster Results'))
else
    sprintf('Too many dimensions to display')
end

% axis(AXSz)