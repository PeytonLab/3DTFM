%% 3D Traction Force Microscopy
% Authors,
% David Podorefsky, Peyton Laboratory at University of Massachusetts Amherst
% &
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

%% Directory Setup
clear; clc;
directory = uigetdir('','Select Output Folder');
cd(directory);
dir_list = dir;
dir_el = numel(dir_list)-2;

%% Save image stacks as vol*.mat t-stacks
c = 2;
for p = 1:dir_el
    cd(dir_list(p+2).name);
    cells_list = dir;
    cells_el = numel(cells_list)-2;
    for s = 1:cells_el     
        cellfol = cells_list(s+2).name;
        cd(cellfol);     
        cell_list = dir;
        cell_name = cell_list(3).name;     
        info = imfinfo(cell_name);

        str = strrep(cell_name,'.tif',''); %t    
        idx1 = strfind(str,'t=')+1;   
        idx2 = length(str);
        for i = 1:idx2-idx1
            t_str(i) = str(idx1+i);
        end
        t = str2double(t_str);
      
        el = length(info); %z
        z = el/t/c;        
        
        mkdir('Volumes');
        addpath(genpath(pwd));   
        
        vol = cell(1,c);        
        for i = 1:t    
            for j = 1:c                  
                for k = 1:z         
                    if k == 1
                       stack = imread(cell_name,j+(z*c)*(i-1));
                      else
                       slice = imread(cell_name,j+(z*c)*(i-1)+c*(k-1));
                       stack = cat(3, stack, slice);
                    end
                end
                vol{j} = stack;
            end
            if i < 10
                volsave = [pwd,'\Volumes\vol0',num2str(i),'.mat'];
              else
                volsave = [pwd,'\Volumes\vol',num2str(i),'.mat'];
            end
            save(volsave,'vol');
        end      
        cd('..');    
    end
    cd('..');
end

%% Displacement Field Calculation (FIDVC)
sSize = [128 128 64]; %subset size
incORcum = 'incremental';
filename{1} = 'vol*.mat'; %name of vol
filename{2} = 1; %channel to run IDVC (DAPI)

%If you must stop FIDVC, set "a" to population and "b" to cell
i = 0; %Leave turned off unless you need to use it, i = 1
if i == 1
   a = 4; %i.e. start from 7th cell in 1st condition
   b = 10;
 else
   a = 1; %normal operation
   b = 1;
end

%Run FIDVC
for p = a:dir_el
    cd(dir_list(p+2).name);
    cells_list = dir;
    cells_el = numel(cells_list)-2;	
    for s = b:cells_el    
        FIDVC_name = cells_list(s+2).name;
        cd([cells_list(s+2).name,'\Volumes']);
        [u, cc, dm] = funIDVC(filename, sSize, incORcum); 
        save('FIDVC.mat','u','cc','dm');
        cd('..\..');
    end
    cd('..');
end

%% Computation Save Space (TFM.mat)
TFM = cell(1,dir_el);
for p = 1:length(TFM)
    cd(dir_list(p+2).name) %Navagate to condition
    cells_list = dir;
    cells_el = numel(cells_list)-2;
    TFM{p} = cell(1,cells_el);
    for s = 1:cells_el
        cd([cells_list(s+2).name,'\Volumes']);    
        vol_list = dir('vol*'); %Directory of volumetric matrices
        vol_el = numel(vol_list);
        TFM{p}{s} = cell(vol_el-1,1);
        for t = 1:vol_el-1
            TFM{p}{s}{t} = cell(1,4);
            load('FIDVC.mat') 
            for i = 1:4
                TFM{p}{s}{t}{1,i} = u{t}{1,i};
            end
        end      
        cd('..\..');
    end
    cd('..');
end
save('TFM.mat','TFM');

%% Displacement Field, Median Drift Registration
for p = 1:length(TFM)
    for s = 1:length(TFM{p})                    
        for t = 1:length(TFM{p}{s})         
            for i = 1:4
                u = TFM{p}{s}{t}{i}
                u_r = u - median(median(median(u)));
                TFM{p}{s}{t}{i} = u_r;
            end
        end
    end
end
save('TFM.mat','TFM')

%% Traction Force Computation (LD-3D-TFM) %%
dm = 8;
for p = 1:length(TFM)
    
    %Hydrogel Material Properties
    if p == 1
        E = 2500; %Pa,  Young's Modulus
    else
        E = 500; %Pa
    end 
    nu = 0.45; %Poisson's Ratio

    cd(dir_list(p+2).name);
    cells_list = dir;
    for s = 1:length(TFM{p})
        cd([cells_list(s+2).name,'\Volumes']);
        
        sizeI = (size(TFM{p}{s}{1}{4})-1)*dm; %Size of 3D image
        
        %Create surface
        [surface{1}{1},surface{1}{2}] = meshgrid(1:dm:sizeI(1)+1,1:dm:sizeI(2)+1); 
        surface{1}{3} = ceil(sizeI(3))/2*ones(size(surface{1}{1}))+1;
        
        %Normals
        normals{1}{1} = zeros(size(surface{1}{1}));
        normals{1}{2} = zeros(size(surface{1}{1}));
        normals{1}{3} = ones(size(surface{1}{1}));
        
        %Define material
        materialModel = 'neohookean';
        materialProps = [E, nu];
        
        %Calculate Surface based on bead pattern
        F = dir('vol*');
        I0 = importdata(F(1).name);
        dataChannel = 1;
        [s1,n1,angle001] = findSurface(I0{dataChannel}, dataChannel, dm, 12,16,0);
        [surface, normals] = calculateSurfaceUi(s1, n1, u);   
        
        %Calculate TFM metrics
        [Fij, Sij, Eij, U, ti, tiPN] = fun3DTFM(u, dm, surface, normals, materialModel, materialProps);
        
        disp([p,s]);
        
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
        cd('..\..');
    end
    cd('..');
end
save('TFM.mat','TFM'); %Contains displacement and traction force information

%% Visualization Setup
mkdir('Figures');
addpath('Figures');
dir_list = dir;
dir_el = numel(dir)-4;

%Parameters
sc = 6.2; %[px/um]
dz = 0.34; %[um/z]
tf = 15; %time frame [min]
dm = 8;

%% Cross-Sectional Animation
fr = 3; %frames per second
gif = 1; %gif creation on (1) or off (0)
for p = 1:dir_el
    mkdir('..\Figures',dir_list(p+2).name) %Create condition folder in figures
    cd(dir_list(p+2).name)
    cells_list = dir;
    cells_el = numel(cells_list)-2;
    for s = 1:cells_el
        mkdir(['..\Figures\',dir_list(p+2).name,'\'],cells_list(s+2).name)
        cd([cells_list(s+2).name,'\Volumes']);
        vol_list = dir('vol*'); %directory of volumetric matrices
        tmax = size(vol_list,1); %maximum time points            
        
        %Determine cross section to view
        ZIdx = zeros(1,tmax);  
        for t = 1:tmax %Determine most intense slice at each time point            
            load(vol_list(t).name)
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
		UMax = ceil(ZMax/dm); %Use slice of max cell intensity       
        
        %Construct cell boundary
        B_cell = cell(tmax,1);
        for t = 1:tmax
            load(vol_list(t).name)
            I = vol{2}(:,:,ZMax);
            for j = 1:3
                I = medfilt2(I); %Median filter, thrice
            end  
            BW = im2bw(I,graythresh(I)); %Binarize        
            [B,L,N,~] = bwboundaries(BW,'noholes'); %Boundaries
            if N > 0 %Cell boundary
                B_list = zeros(N,1);
                for b = 1:N
                    B_list(b) = size(B{b},1);
                end
                [~,Idx] = max(B_list);
                B_cell{t} = B{Idx}; %Boundary for time point
            else
                disp(['No boundary detected for s=',num2str(s),' at t=',num2str(t)]);
            end 
        end
        
        %Plot index
        sizeI = size(vol{1});   
        plotIdx = cell(1,3);      
        for i = 1:2, plotIdx{i} = 1:dm:sizeI(i)+1;
            plotIdx{i} = plotIdx{i}./sc; 
        end
        plotIdx{3} = 1:dm:sizeI(3)+1; 
        plotIdx{3} = plotIdx{3}.*dz;       
        
        %Color map bounds          
        for t = 1:length(TFM{p}{s})
            cmin1(t) = min(min(TFM{p}{s}{t}{1,4}(:,:,UMax)));  %Min
            cmax1(t) = max(max(TFM{p}{s}{t}{1,4}(:,:,UMax)));  %Max
            cmin2(t) = min(min(TFM{p}{s}{t}{2,4}));
            cmax2(t) = max(max(TFM{p}{s}{t}{2,4}));
        end
        cmin1 = floor(min(cmin1))./(sc*tf).*60; %um/hr
        cmax1 = ceil(max(cmax1))./(sc*tf).*60; %um/hr
        cmin2 = floor(min(cmin2))*sc^2./(tf*60); %Pa/s      
        cmax2 = ceil(max(cmax2))*sc^2./(tf*60); %Pa/s
        
        for t = 1:length(TFM{p}{s})             
            %% Displacements
            if t == 1 && s == 1 && p == 1
               figure
            end 
            figure(1);
            ax1 = gca;
            [~,h] = contourf(ax1, plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{1,4}(:,:,UMax)./(sc*tf).*60,50);
            c = colorbar;
            caxis([cmin1,cmax1]);
            set(h,'linestyle','none'); axis image;
            
            %Plot Boundary
            if ~isempty(B_cell{t+1})    
                hold on
                plot(B_cell{t+1}(:,2)./sc,B_cell{t+1}(:,1)./sc,'Color','White','LineWidth',2)
                hold off              
                
                %Mask cell projection onto plot
                M = sizeI(1:2);
                green = cat(3, zeros(M), ones(M), zeros(M));
                load(vol_list(t+1).name)
                I = vol{2}(:,:,ZMax);
                hold on
                x0 = [plotIdx{1}(1) plotIdx{1}(end)];
                y0 = [plotIdx{2}(1) plotIdx{2}(end)];
                g = image(x0,y0,green);       
                set(g, 'AlphaData',I)
                hold off
            end
            
            %Displacement arrows
            W = [TFM{p}{s}{t}{1,1}(:,:,UMax),TFM{p}{s}{t}{1,2}(:,:,UMax)]./(sc*tf).*60;
            X = zeros(1,2);
            k = 2; %croping passes
            for i = 1:2              
                A = W(i);
                B = plotIdx{i};
                for j = 1:k
                    A = A(2:2:size(A,1)-1,:); %Arrows
                    A = A(:,2:2:size(A,2)-1);
                    B = B(2:2:end-1); %Plot index
                end
                W(i) = A;
                X(i) = B;   
            end            
            u = W(1); v = W(2);
            x = X(1); y = X(2); 
            [x,y] = meshgrid(x,y);           
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            colormap(ax1,jet); %contour colormap
            set(gca,'Ydir','reverse','FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold');
            title('Displacement Field','fontsize',13);
            ylabel(c,'Bead displacement (um/hr)','fontsize',12);
            xlabel('um'); ylabel('um');
            axis([1 sizeI(2) 1 sizeI(1)]./sc);
            box off;
            set(gcf,'color','white');
            drawnow
            
            %Save .gif and .tif
            outputdir = ['../../../Figures/',dir_list(p+2).name,'/',cells_list(s+2).name,'/'];            
            if gif == 1
                %Capture frame
                frame = getframe(1);
                im1{t} = frame2im(frame);
                [G,map] = rgb2ind(im1{t},256);
                
                filename = [outputdir,'dis.gif'];
                if t == 1
                   imwrite(G,map,filename,'gif', 'LoopCount',inf,'DelayTime',1/fr);
                  else
                   imwrite(G,map,filename,'gif','WriteMode','append','DelayTime',1/fr);
                end
            end
            name1 = [outputdir,'dis0',num2str(t),'.tif'];
            saveas(gcf,name1)
                        
            %% Traction Forces
            if t == 1 && s == 1 && p == 1
               figure
            end
            figure(2);
            ax2 = gca;
            [~,h] = contourf(ax2, plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{2,4}*sc^2./(tf*60),50); 
            c = colorbar;
            caxis([cmin2,cmax2]);
            set(h,'linestyle','none'); axis image;
            
            if ~isempty(B_cell{t+1})
                %Plot boundary
                hold on
                plot(B_cell{t+1}(:,2)./sc,B_cell{t+1}(:,1)./sc,'color','white','LineWidth',3);
                hold off

                %Plot cell
                hold on
                g = image(x0,y0,green);       
                set(g,'AlphaData',I)
                hold off
            end
            
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            colormap(ax2,jet)
            set(gca,'Ydir','reverse','FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold')
            title('Traction Field','fontsize',13)
            ylabel(c,'Traction force (Pa/s)','fontsize',12)
            xlabel('um'); ylabel('um')
            axis([1 sizeI(2) 1 sizeI(1)]./sc)
            box off
            set(gcf,'color', [1 1 1]);
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
            
            name2 = [outputdir,'trac0',num2str(t),'.tif'];
            saveas(gcf,name2)
        end
        cd ('..\..')
    end
    cd('..')
end

%% Box plot 
load TFM.mat
p = length(TFM); %number of conditions
CellDisplacement = cell_name(1,p);	 %displacements of cells
CellTraction = cell_name(1,p);  %traction of cells 

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
			U(k) = abs(mean2(TFM{i}{j}{k}{1,4})); %mean displacement in time frame in 3D space
			F(k) = abs(mean2(TFM{i}{j}{k}{2,4})); %traction force								
        end
        
        %Store vectors for each cell in a cell
        if j == 1
		   CellDisplacement{i} = U;
		   CellTraction{i} = F;
        else 			
           CellDisplacement{i} = horzcat(CellDisplacement{i},U);
		   CellTraction{i} = horzcat(CellTraction{i},F);
        end
    end
    
    %If values are zero, replace w/ nan (needed for box plot)
	CellDisplacement{i}(CellDisplacement{i}==0)=nan;
	CellTraction{i}(CellTraction{i}==0)=nan;	
end

%Group into population averages
for i = 1:p
    U_P{i} = CellDisplacement{i};
    U_P{i} = reshape(U_P{i},[size(U_P{i},1)*size(U_P{i},2),1]);
    U_P{i} = U_P{i}(~isnan(U_P{i}));
    
    F_P{i} = CellTraction{i};
    F_P{i} = reshape(F_P{i},[size(F_P{i},1)*size(F_P{i},2),1]);
    F_P{i} = F_P{i}(~isnan(F_P{i}));
end

%Convert cells of cell populations to a matrix w/ nan for space
for i = 1:p
    d(i) = length(U_P{i});
end
dmax = max(d); %maximum spaces
PhaseDisplacement = NaN(dmax,p); %fill with empty space
PhaseTraction = NaN(dmax,p);

for i = 1:p
    for	j = 1:length(U_P{i})
		PhaseDisplacement(j,i) = U_P{i}(j); %Fill space with displacements	
    end
    for k = 1:length(F_P{i})
        PhaseTraction(k,i) = F_P{i}(k);
    end
end

%Plot
for i=1:p %For each population
	figure
    x = boxplot(CellDisplacement{i}./(sc*tf)*60);
    xlabel('Cell')
    ylabel('Speed (um/hr)')
    str=['Bead Displacements (P',num2str(i),')'];
    title(str)
    set(x,{'linew'},{1})
    set(gcf,'color','white')
    box off

    figure
	x = boxplot(CellTraction{i}.*sc^2./(tf*60));
    xlabel('Cell')
    ylabel('Force (Pa/s)')
    str=['Cellular traction force (P',num2str(i),')'];
    title(str)
    set(x,{'linew'},{1})
    set(gcf,'color','white')
    box off
end

%Population displacements
figure
x = boxplot(PhaseDisplacement./(sc*tf)*60,'Labels',{'20wt,0.5mM','20wt,2mM','5wt,0.5mM','5wt,2mM'});
set(x(7,:),'Visible','off')
ylim([0 2])
set(gca,'YTick',0:1:8,'FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold')
ylabel('Speed (um/hr)','fontsize',12)
title('Bead Displacement','fontsize',13)
set(x,{'linew'},{1})
set(gcf,'color','white')
box off

%Population tractions
figure
x = boxplot(PhaseTraction.*sc^2./(tf*60),'Labels',{'20wt,0.5mM','20wt,2mM','5wt,0.5mM','5wt,2mM'});
set(x(7,:),'Visible','off')
ylim([0 70])
set(gca,'YTick',0:10:100,'FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold')
ylabel('Force (Pa/s)','fontsize',12)
title('Cellular Traction Forces','fontsize',13)
set(x,{'linew'},{1})
set(gcf,'color','white')
box off

%DJP, 1/11/2018