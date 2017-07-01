%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Traction Force Microscopy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If used please cite:
% David Podorefsky, Peyton Laboratory at University of Massachusetts Amherst
% &
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2
clear; close all; clc
%% Workspace set up (directory of Stack folder of conditions)
folder_name = uigetdir('','Select Experiment Folder');
cd(folder_name)
masterfol = dir;
masfolel = numel(masterfol)-2;
%% Image Construction %%
% Converts a time series of image stacks from .tiff files (RGB) to three dimensional
% arrays and stores each time point in a folder as a .mat file for the DAPI
% and EGFP channels
% NOTE: Names of cell folders, conditions, and cell names must start with a number

c = 2; %number of channels
z = 96; %slices in stack

for p = 1:masfolel
    cd(masterfol(p+2).name) %Navagate to condition within experiment folder
    cellsfol = dir;
    cellfolel = numel(cellsfol)-2;
    for s=1:cellfolel
        Cell_Volume = cellsfol(s+2).name
        cd(Cell_Volume) %Navagate to cell        
        mkdir('Volumes') %Create folder to store     
        addpath(genpath(pwd))  
        scene = dir; %Create directory variable
        el = numel(scene)-3; %Determine number of elements in scene  
        t = el/c/z; %timepoints in set
        for j = 1:t %time point
            vol = cell(1,c);	
            for k=1:c   %channel
                a = 2+k; %starting index number **use 2+k if name starts with number
                if k == 1
                   RGB = 3; %DAPI
                elseif k == 2 
                   RGB = 2; %EGFP
                elseif k == 3
                   RGB = 1; %red
                end			 
                image = scene(a+(z*c)*(j-1)).name;
                stack = imread(image); %reads first slice of time point         
                stack = stack(:,:,RGB);     
                for i=2:z %z slice
                    index = a+(z*c)*(j-1)+c*(i-1);
                    image = scene(index).name;
                    temp = imread(image);
                    temp = temp(:,:,RGB);  %Comment if not (~,~,3)
                    stack = cat(3, stack, temp);
                end
                vol{k} = stack;
            end
            if j < 10
               volsave = [pwd,'\Volumes\vol0',num2str(j),'.mat'];
             else
               volsave = [pwd,'\Volumes\vol',num2str(j),'.mat'];
            end
            save(volsave,'vol')    
        end
        cd ..    
    end
    cd ..
end

%% Hilbert Space to Store Computations (TFM.mat)
masterfol = dir;
masfolel = numel(masterfol)-2;
TFM = cell(1,masfolel);
for p = 1:length(TFM)
    cd(masterfol(p+2).name) %Navagate to condition
    cellsfol = dir;
    cellfolel = numel(cellsfol)-2;
    TFM{p} = cell(1,cellfolel);
    for s = 1:cellfolel
        cd(cellsfol(s+2).name) %Navagate to cell        
        scene = dir; %Create directory variable
        el = numel(scene)-3; %Determine number of elements in scene  
        tp = el/c/z - 1; %timepoints in set (dt)
        for t=1:tp
            if t == 1
               TFM{p}{s} = cell(tp,1);
            end
            TFM{p}{s}{t} = cell(1,4);
        end
        cd ..
    end
    cd ..
end
save('TFM.mat','TFM')

%% Displacement Field Calculation (FIDVC) %%
% Estimates displacements between beads in successive time pointns via FIDVC algorithm
% NOTE: Image stack should be at least 1.5 times subset size (128x1.5=192,64x1.5=96)
% The default subset size is 128x128x64, so we recommend that the minimum input volume size should be 192x192x96.
% The size of the input image stack should be divisible by 0.5 times the size of the subset.
% 512/(128/2)=8, 96/(64/2)=3, 128/(64/2)=4 (valid sizes)
% 

sSize = [128 128 64]; %subset size
incORcum = 'incremental';
filename{1} = 'vol*.mat'; %name of vol
filename{2} = 1; %channel to run IDVC (DAPI)

%If you must stop FIDVC, set "a" to condition and "b" to cell
i = 0; %Leave turned off unless you need to use it, i = 1

if i == 1
    a = 1; %i.e. start from 7th cell in 1st condition
    b = 7;
 else
    a = 1; %normal operation
    b = 1;
end

masterfol = dir;
for p = a:length(TFM)
    cd(masterfol(p+2).name) %Navagate to condition
    cellsfol = dir;
    cellfolel = numel(cellsfol)-2;	
    for s = b:cellfolel    
        FIDVC_name = cellsfol(s+2).name
        cd(cellsfol(s+2).name)
        cd('Volumes')
        [u, cc, dm] = funIDVC(filename, sSize, incORcum); 
        save('FIDVC.mat','u','cc','dm');
        cd ('..\..')
    end
    cd ..
end

%Assemble displacements into TFM.mat file
masterfol = dir;
for p=1:length(TFM)
    cd(masterfol(p+2).name) %Navagate to condition
    cellsfol = dir;
    cellfolel = numel(cellsfol)-2;
    for s = 1:cellfolel
        cd(cellsfol(s+2).name)
        cd('Volumes')
        load('FIDVC.mat')
        for t = 1:length(u)
            for i = 1:4
                TFM{p}{s}{t}{1,i} = u{t}{1,i};
            end
        end
        cd ('..\..')
    end
    cd ..
end
save('TFM.mat','TFM')

%% Traction Force Computation (LD-3D-TFM) %%
% FUNCTIONRUNTFM calculates 3D Traction Force Microscopy metrics
% [Fij, Sij, Eij, Uhat, ti, tiPN] = functionRunTFM(u,dm, surface, normals,
% materialProps, materialModel) calculates Deformation gradient Fij,
% Cauchy stress Sij, Lagrangian strain Eij, Strain energy density Uhat,
% surface traction vector ti, and in-plane and out-of-plane
% components of surface tractions
clear;
masterfol = dir;
load TFM.mat

%Spacing (dm)
cd(masterfol(3).name)
cellsfol = dir;
cd(cellsfol(3).name)
cd Volumes
load vol01.mat
dm=size(vol{1},1)/(size(TFM{1}{1}{1}{1},1)-1); %meshgrid spacing (8 by default)
cd('..\..\..')

E = 2000; %Pa,  Young's Modulus
nu = 0.45; %Poisson's Ratio

for p = 1:length(TFM)
    cd(masterfol(p+2).name) %Navagate to condition
    cellsfol = dir;
    for s = 1:length(TFM{p})
        cd(cellsfol(s+2).name) %Navagate to cell
        cd('Volumes')
        u = TFM{p}{s};
        sizeI = (size(u{1}{4})-1)*dm; %Size of 3D image

        %Create surface
        [surface{1}{1},surface{1}{2}] = meshgrid(1:dm:sizeI(1)+1,1:dm:sizeI(2)+1); 
        surface{1}{3} = ceil(sizeI(3))/2*ones(size(surface{1}{1}))+1;

        %Normals
        normals{1}{1} = zeros(size(surface{1}{1}));
        normals{1}{2} = zeros(size(surface{1}{1}));
        normals{1}{3} = ones(size(surface{1}{1}));

        %Define material
        materialModel = 'neohookean';
        materialProps = [E, nu];  % [Young's Modulus, Poisson's Ratio]
 
        %Calculate Surface based on bead pattern
        filePrefix = 'vol*';
        F = dir(filePrefix);
        I0 = importdata(F(1).name);
        dataChannel = 1;
        [s1,n1,angle001] = findSurface(I0{dataChannel}, dataChannel, dm, 12,16,0);
        [surface, normals] = calculateSurfaceUi(s1, n1, u);
    
        %Calculate TFM metrics
        [Fij, Sij, Eij, U, ti, tiPN] = fun3DTFM(u, dm, surface, normals, materialModel, materialProps);
        disp([p,s])
        
        for i = 1:9
        gradU = sqrt((Fij{1}{1,1}-1).^2 + (Fij{1}{1,2}).^2 + (Fij{1}{1,3}).^2 + ...
            (Fij{1}{2,1}).^2 + (Fij{1}{2,2}-1).^2 + (Fij{1}{2,3}).^2 + ...
            (Fij{1}{3,1}).^2 + (Fij{1}{3,2}).^2 + (Fij{1}{3,3}-1).^2);
        end

        %Store in TFM file
        for t = 1:length(ti)
            for i = 1:4
                TFM{p}{s}{t}{2,i} = ti{t}{1,i};
            end
        end
        
        cd ('..\..')
    end
    cd ..
end
save('TFM.mat','TFM'); %Contains displacement and traction force information

%% Visualization %%
mkdir('Figures') %Figures folder for storage
addpath('Figures')
%% 2D boundary plot
load TFM.mat
masterfol = dir;
masfolel = numel(masterfol)-4;

sc = 6.2; %[px/um]
dz = 0.34; %[um/z]
tf = 30; %[min] time frame
fr = 3; %frames per second
gif = 1; %gif creation on (1) or off (0)

for p = 1:masfolel
    cd(masterfol(p+2).name)
    cellsfol = dir;
    cellfolel = numel(cellsfol)-2;
    mkdir('../Figures',masterfol(p+2).name) %Create condition folder in figures
    for s = 1:cellfolel
        str = ['../Figures/',masterfol(p+2).name,'/'];
        mkdir(str,cellsfol(s+2).name)
        cd(cellsfol(s+2).name)
        cd('Volumes')
        V = dir('vol*'); %directory of volumetric matrices
        tmax = size(V,1); %maximum time points        
        
        %Determine cross section to view
        ZIdx = zeros(1,tmax); 
        for t = 1:tmax %Determine most intense slice at each time point            
            load(V(t).name)
            I = vol{2}; %EGFP channel
            zmax = size(I,3);
            Zint = zeros(zmax,1);
            for z = 1:zmax
                Zint(z) = sum(sum(I(:,:,z))); %intensity of each slice
            end
            [~,Idx] = max(Zint); %max intensity slice for time point
            ZIdx(t) = Idx;
        end
        ZMax = ceil(mean(ZIdx));        
        
        %Determine displacement slice to view
        
		%1.Find max displacement slice
		% UIdx = zeros(length(TFM{p}{s}),1);     
        % for t = 1:length(TFM{p}{s})
            % zmax = size(TFM{p}{s}{t}{1,4},3);
            % zmean = zeros(zmax,1);
            % for z = 1:zmax
                % zmean(z) = mean2(TFM{p}{s}{t}{1,4}(:,:,z));
            % end
            % [~,I]=max(zmean);
            % UIdx(t)=I;
        % end
        % UMax=ceil(mean(UIdx));
        
		%2. Use slice of max cell intensity
		UMax=ceil(ZMax/dm); %
        
        %Construct Boundary
        BS = cell(tmax,1);
        for t = 1:tmax
            load(V(t).name)
            I = vol{2}(:,:,ZMax);
            % figure
            % imshow(I)          
            BW = im2bw(I,0.3);
            % figure 
            % imshow(BW)            
            [B,~] = bwboundaries(BW,'noholes'); %Make boundaries
            %Find max (the true) boundary
            if ~isempty(B) %If boundaries exist, which boundary to use
                BIdx = zeros(size(B,1),1);
                for b = 1:size(B,1)
                    BIdx(b) = size(B{b},1);
                end
                [~,Idx] = max(BIdx);
                BS{t} = B{Idx}; %Boundary for time point
              else
                disp('Uh oh');
            end  
        end

        %The Plot
        sizeI = size(vol{1});
        dm=size(vol{1},1)/(size(TFM{1}{1}{1}{1},1)-1); %meshgrid spacing (8 by default)
        
        plotIdx = cell(1,3);
        for i = 1:2, plotIdx{i} = 1:dm:sizeI(i)+1; plotIdx{i} = plotIdx{i}./sc; end
		plotIdx{3} = 1:dm:sizeI(3)+1; plotIdx{3} = plotIdx{3}.*dz;
        
        %Determine color map bounds, displacement and traction        
            %Max
            for t = 2:length(TFM{p}{s})
                cmax1(t) = max(max(TFM{p}{s}{t}{1,4}(:,:,UMax)));
                cmax2(t) = max(max(TFM{p}{s}{t}{2,4}));
            end
            cmax1=ceil(max(cmax1))./(sc*tf).*60; %microns/hr
            cmax2=ceil(max(cmax2))*sc^2./(tf*60); %Pa/s

            %Min
            for t = 2:length(TFM{p}{s})
                cmin1(t) = min(min(TFM{p}{s}{t}{1,4}(:,:,UMax)));
                cmin2(t) = min(min(TFM{p}{s}{t}{2,4}));
            end
            cmin1=floor(min(cmin1))./(sc*tf).*60; %microns/hr
            cmin2=floor(min(cmin2))*sc^2./(tf*60); %Pa/s
        
        for t = 1:length(TFM{p}{s})
            %Plot displacement map
            if t == 1 && s==1 && p==1
               figure
            end 
            figure(1);
            ax1 = gca;
            [~,h] = contourf(ax1,plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{1,4}(:,:,UMax)./(sc*tf).*60,50);
            c = colorbar;
            caxis([cmin1,cmax1])
            set(h,'linestyle','none'); axis image
            
            %Plot boundary
            hold on
            plot(BS{t+1}(:,2)./sc,BS{t+1}(:,1)./sc,'w','LineWidth',2)
            hold off

            %Mask cell projection onto plot
            M = sizeI(1:2);
            green = cat(3, zeros(M), ones(M), zeros(M));
            load(V(t+1).name)
            I = vol{2}(:,:,ZMax);
            %INITPSF = ones(size(S)); %blind deconvolution
            %[J,PSF] = deconvblind(S,INITPSF);
            hold on
            x0 = [plotIdx{1}(1) plotIdx{1}(end)];
            y0 = [plotIdx{2}(1) plotIdx{2}(end)];
            g=image(x0,y0,green);       
            set(g, 'AlphaData',I)
            hold off 
            
            %Arrows
            [u,v] = gradient(TFM{p}{s}{t}{1,4}(:,:,UMax)./(sc*tf).*60);
			
            %Crop arrows and plot index
            k = 2; %croping passes
                A = u;
                for i = 1:k
                    A = A(2:2:size(A,1)-1,:);
                    A = A(:,2:2:size(A,2)-1);
                end
                u = A;

                A = v;
                for i = 1:k
                    A = A(2:2:size(A,1)-1,:);
                    A = A(:,2:2:size(A,2)-1);
                end
                v = A;

                B = plotIdx{1};
                for i = 1:k
                    B = B(2:2:end-1);
                end
                x = B;

                B = plotIdx{2};
                for i = 1:k
                    B = B(2:2:end-1);
                end
                y = B;         
            [x,y] = meshgrid(x,y);
            
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            colormap(ax1,jet) %contour colormap
            set(gca,'Ydir','reverse')           
            title('Displacement Field')
            ylabel(c,'Bead displacement (\mum/hr)')
            xlabel('\mum'); ylabel('\mum')
            axis([1 sizeI(2) 1 sizeI(1)]./sc)
            box off
            set(gcf, 'color', [1 1 1]);
            drawnow
            
            %Save gif and .tiff
            outputdir = ['../../../Figures/',masterfol(p+2).name,'/',cellsfol(s+2).name,'/'];            
            
            if gif == 1
                %Capture frame
                frame = getframe(1);
                im1{t} = frame2im(frame);
                [G,map] = rgb2ind(im1{t},256);
                
                filename = [outputdir,'dis.gif'];
                if t == 1;
                   imwrite(G,map,filename,'gif', 'LoopCount',inf,'DelayTime',1/fr);
                  else
                   imwrite(G,map,filename,'gif','WriteMode','append','DelayTime',1/fr);
                end
            end
  
            dname = [outputdir,'dis0',num2str(t),'.tif'];
            saveas(gcf,dname)
            
            %Plot Traction
            if t == 1 && s==1 && p==1
               figure
            end
            figure(2);
            ax2 = gca;
            [~,h] = contourf(plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{2,4}*sc^2./(tf*60),50); 
            c = colorbar;
            caxis([cmin2,cmax2])
            set(h,'linestyle','none'); axis image
            
            %Plot boundary
            hold on
            plot(BS{t+1}(:,2)./sc,BS{t+1}(:,1)./sc,'w','LineWidth',3)
            hold off
            
            %Plot cell
            hold on
            g=image(x0,y0,green);       
            set(g, 'AlphaData',I)
            hold off 
            
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            colormap(ax2,jet)
            set(gca,'Ydir','reverse') 
            title('Traction Field')
            ylabel(c,'Traction Force (Pa/s)')
            xlabel('\mum'); ylabel('\mum')
            axis([1 sizeI(2) 1 sizeI(1)]./sc)
            box off
            set(gcf, 'color', [1 1 1]);
            drawnow
            
            if gif == 1
                %Capture frame             
                frame = getframe(2);
                im2{t} = frame2im(frame);
                [G,map] = rgb2ind(im2{t},256);
                
                %Save gif and .tiff
                filename = [outputdir,'trac.gif'];            
                if t == 1;
                    imwrite(G,map,filename,'gif', 'LoopCount',inf,'DelayTime',1/fr);
                  else
                    imwrite(G,map,filename,'gif','WriteMode','append','DelayTime',1/fr);
                end    
            end
            
            tname = [outputdir,'trac0',num2str(t),'.tif'];
            saveas(gcf,tname)
        end
        cd ('..\..')
    end
    cd('..')
end

%% Box plot 
load TFM.mat
p = length(TFM); %number of conditions
US = cell(1,p);	 %displacements of cells
FS = cell(1,p);  %traction of cells 

for i = 1:p %condition number
	s = length(TFM{i});	 %number of cells
	
    %Determine max time points for storage
    for j = 1:s  %cell number
        t(j) = length(TFM{i}{j}); %time points for each cell
    end 
    tmax = max(t); %max time points
    
    %Compute average displacements/force per cell per time frame and store in vectors
    for j = 1:s		
		t = length(TFM{i}{j});
		U = zeros(tmax,1);
        F = zeros(tmax,1);
        for k = 1:t						
			U(k) = mean2(TFM{i}{j}{k}{1,4}); %mean displacement in time frame in 3D space
			F(k) = mean2(TFM{i}{j}{k}{2,4}); %traction force								
        end
        
        %Store vectors for each cell in a cell
        if j == 1
		   US{i} = U;
		   FS{i} = F;
		 else 			
			US{i} = horzcat(US{i},U);
			FS{i} = horzcat(FS{i},F);
        end
    end
    
    %If values are zero, replace w/ nan (needed for box plot)
	US{i}(US{i}==0)=nan;
	FS{i}(FS{i}==0)=nan;	
end

%Group into population averages
for i = 1:p
    U_P{i} = US{i};
    U_P{i} = reshape(U_P{i},[size(U_P{i},1)*size(U_P{i},2),1]);
    U_P{i} = U_P{i}(~isnan(U_P{i}));
    
    F_P{i} = FS{i};
    F_P{i} = reshape(F_P{i},[size(F_P{i},1)*size(F_P{i},2),1]);
    F_P{i} = F_P{i}(~isnan(F_P{i}));
end

%Convert cells of cell populations to a matrix w/ nan for space
for i = 1:p
    d(i) = length(U_P{i});
end
dmax = max(d); %maximum spaces
UP = zeros(dmax,p); %fill with empty space
FP = zeros(dmax,p);
for i = 1:p
	for	j = 1:length(U_P{i})
		UP(j,i) = U_P{i}(j); %Fill space with displacements
		FP(j,i) = F_P{i}(j);
	end
end

UP(UP==0)=nan;
FP(FP==0)=nan;

%Condition index
for i = 1:p
    N(i) = masterfol(i+2).name;
end

%Populations' Displacements
figure
x = boxplot(UP./(sc*tf)*60,'Labels',{'Control','VEGF'});
ylim([0 6.5])
set(gca,'YTick',0:1:7)
ylabel('Speed (um/hr)')
title('Bead displacement')

%Populations' Tractions
figure
x = boxplot(FP.*sc^2./tf.*60,'Labels',{'Control','6uM CytoD'});
ylim([0 1300])
set(gca,'YTick',0:200:1300)
ylabel('Force (Pa/hr)')
title('Traction forces')

%For each population
for i=1:p
	figure
    boxplot(US{i}./(sc*tf)*60)
    xlabel('Cell')
    ylabel('Speed (\mum/hr)')
    str=['Bead displacements (P',num2str(i),')'];
    title(str)    
    
    figure
	boxplot(FS{i}.*sc^2./tf.*60/10^2)
    xlabel('Cell')
    ylabel('Force (Pa/hr)')
    str=['Cellular traction force (P',num2str(i),')'];
    title(str)
end

%% 3D Recontruction of Cell w/ displacement vectors & contour map of traction  #still in dev
clear U
load TFM.mat
for p = 1:length(TFM)
    for l = 1:length(TFM{1,p})
        if p==1
          cell_name = num2str(masterscene(l+2).name)
        else
          cell_name = num2str(masterscene(l+2+length(TFM{1,1})).name)
        end
        
        tp = size(TFM{1,p}{1,l},1)+1; %Determine time points in cell 4D image 
        for t = 1:tp
            load([pwd,'\',cell_name,'\Vol\vol',num2str(t),'.mat'])
            I1 = vol{2}; %reads image
            I2 = zeros(size(I1));
            Zint = zeros(size(I1,3),1);
            for i=1:size(I1,3)   
                Im=I1(:,:,i);
                BW = im2bw(Im,0.3);
                [B,L] = bwboundaries(BW,'noholes');    
                Zint(i) = sum(sum((double(BW)))); %intensity to determine projection    
                if ~isempty(B) %If boundaries exist       
                   maxB=size(B{1},1); %Size of first
                   max_index=1;
                   for j=2:numel(B)    
                       if size(B{j},1)>maxB         
                          max_index=j;
                          maxB=size(B{j},1);
                       end
                   end
                   Bb = B{max_index}; 
                   I3 = zeros(size(Im));
                   for k = 1:size(Bb,1)      
                       I3(Bb(k,2),Bb(k,1))=1;
                   end
                   I2(:,:,i)=I3;   
                end   
            end
            [~,Idx] = max(Zint);

            [~,n,zz] = size(I1);   %create 3D grid to plot
            [X,Y,Z] = meshgrid(1:n,1:n,1:zz);

            figure
            rec = patch(isosurface(X,Y,Z,I2,0.01)); %reconstruct
            rec.FaceColor = 'g';
            rec.EdgeColor = 'none';

            hold on

            %Plot displacement vectors
            if t==1    %Determine spacing of original displacent matrix
               [~,nn,zz]=size(TFM{1,p}{1,l}{1}{1});
            else
               [~,nn,zz]=size(TFM{1,p}{1,l}{t-1}{1});
            end

            %Crop initiation
            cn = 4; %crop number
            for kk = 1:cn 
                for ii=1:3
                    if t==1 && kk==1
                       A = TFM{1,p}{1,l}{1}{1,ii};
                    elseif kk==1
                       A = TFM{1,p}{1,l}{t-1}{1,ii};
                    else
                       A = R{ii};               
                    end          

                    %Remove excess arrows%%
                    [r,c,pg]=size(A); %determine size of displacement matrix      

                    %Odd slices only in A cropped volume on first pass
                    if kk==1
                        for jjj=1:pg					
                            if mod(jjj,2)~=0 
                                if jjj==1
                                    Ac=A(:,:,jjj);
                                else 
                                    Ac=cat(3,Ac,A(:,:,jjj));
                                end
                            end
                        end
                    pgc=size(Ac,3);
                    clear A
                    A=Ac;
                    end		

                    for jj=1:pgc
                        Aa=A(:,:,jj);
                        %Mark rows for deletion
                        for i = 1:r
                            if mod(i,2)==0
                               Aa(i,:) = 0;
                            end 
                        end 
                        %Mark columns for deletion
                        for j = 1:c
                            if mod(j,2)==0
                               Aa(:,j) = 0;
                            end 
                        end
                        %Delete null rows and columns 
                        Aa(~any(Aa,2),:) = [];  %rows
                        Aa(:,~any(Aa,1)) = [];  %columns
                        %Create 3D matrix of resized displacement
                        if jj==1
                           Bb = Aa;
                        else
                           Bb = cat(3,Bb,Aa);
                        end	
                    end
                    if ii==1 %Assign results to matrices
                        U=Bb;
                    elseif ii==2
                        V=Bb;
                    else
                        W=Bb;
                    end  
                end
                R = {U V W};		
            end

            [~,cc,pgg]=size(U);
            sx=ceil(nn/cc);
            sy=ceil(nn/cc);
            sz=ceil(zz/pgg);     
            [X,Y,Z] = meshgrid(1:dm*sx:n+1,1:dm*sy:n+1,1:dm*sz:p+1); %Image spacing to displacement spacing     
            quiver3(X,Y,Z,U,V,W,'color',[1,0,0],'linewidth',1)

            %Plot traction force contour map
            X = 1:dm:n+1;
            Y = 1:dm:n+1;
            Z = ceil(Idx/dm+1);
            if t==1
                contourf(X,Y,TFM{1,p}{1,l}{1}{2,4},50,'linestyle','none')
            else
                contourf(X,Y,TFM{1,p}{1,l}{t-1}{2,4},50,'linestyle','none')
            end
            h = colorbar;
            ylabel(h, 'pN/m^2')
            %Adjust colors to be same for each in a time point
            colormap(jet)

            %Adjust plot
            view(3);
            axis([1 512 1 512 0 96])
            daspect([6.2 6.2 2.94])
            camlight
            lighting gouraud

            hold off
        end
    end
end