%% 3D Traction Force Microscopy
% Authors,
% David Podorefsky, Peyton Laboratory at University of Massachusetts Amherst
% &
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2
%
% DJP, 5/27/2018

%% Directory Setup
clear; clc;
directory = uigetdir('','Select Output Folder');
cd(directory);
dir_list = dir;
dir_el = numel(dir_list)-2;

%% Save image stacks as vol.mat t-stacks
c = 2; %or 3
for p = 1:dir_el
    cd(dir_list(p+2).name);
    cells_list = dir;
    cells_el = numel(cells_list)-2;
    for s = 1:cells_el     
        cellfol = cells_list(s+2).name;
        cd(cellfol);     
        cell_list = dir;
        cell_name = cell_list(3).name     
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

%Run FIDVC
for p = 1:dir_el
    cd(dir_list(p+2).name);
    cells_list = dir;
    cells_el = numel(cells_list)-2;	
    for s = 1:cells_el    
        FIDVC_name = cells_list(s+2).name
        cd([cells_list(s+2).name,'\Volumes']);
        [u, cc, ~] = funIDVC(filename, sSize, incORcum); 
        save('FIDVC.mat','u','cc');
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

%% Displacement Field Registration

% %CytoD space
% dt1 = [35 35 35 35 35 35 30 30 30 30 30 30 30 35 35]./60; %hr
% dt2 = [35 35 35]./60; %hr
% dt3 = [30 30 30 30 30 30]./60; %hr
% dt4 = [35]./60; %hr
% dt = {dt1,dt2,dt3,dt4};
% 
% %DOE Phase space
% dt = [30 30 15 15 15 30 15 30 15 15 15 19]./60;  %/hr

%MMP14(TR) Space
dt = 34/60; %/hr

for p = 1:length(TFM)
    for s = 1:length(TFM{p})                    
        for t = 1:length(TFM{p}{s})         
            for i = 1:4
                u = TFM{p}{s}{t}{i};
                
                %square crop                
                %outsize = [65 65 13]-[2 2 2]; 
                outsize = [129 129 13] - [4 4 2];
                
                sizeI = size(u);
                diffsize = sizeI-outsize;               
                padL = ceil(diffsize./2);
                padR = floor(diffsize./2);    
                u = u(1+padL(1):sizeI(1)-padR(1),1+padL(2):sizeI(2)-padR(2),1+padL(3):sizeI(3)-padR(3));                             
                [a,b,~] = size(u);
                if a > b
                   u = u(1+a-b:end,:,:);
                  elseif b > a
                     u = u(:,1+b-a:end,:);
                end
                            
                u_r = u - median(median(median(u))); %subtract median           
                
                %u_r = u_r./dt{p}(s); %px to px/hr                        
                %u_r = u_r./dt(p); %px to px/hr
                u_r = u_r./dt; %px to px/hr
                                
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
    Evec = [2 2 2 2 3.75 3.75 3.75 3.75 3.75 3.75 3.75 1.25 1.25 1.25 1.25]; %Young's Modulus, kPa
    
    E = Evec(p)*1000; %Pa
    E = 2000;
    
    nu = 0.45; %Poisson's Ratio

    cd(dir_list(p+2).name);
    cells_list = dir;
    for s = 1:length(TFM{p})
        cd([cells_list(s+2).name,'\Volumes']);
        u = TFM{p}{s};
        sizeI = (size(u{1}{1,4})-1)*dm; %Size of 3D image
        
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
        
        %Match size with bead image
        outsize = size(I0{dataChannel});
        diffsize = outsize - sizeI;      
        padL = ceil(diffsize./2);
        padR = floor(diffsize./2);        
        I0{dataChannel} = I0{dataChannel}(1+padL(1):outsize(1)-padR(1),1+padL(2):outsize(2)-padR(2),1+padL(3):outsize(3)-padR(3));
                
        [s1,n1,angle001] = findSurface(I0{dataChannel}, dataChannel, dm, 12,16,0);
        [surface, normals] = calculateSurfaceUi(s1, n1, u);   
        
        %Calculate TFM metrics
        [Fij, Sij, Eij, U, ti, tiPN] = fun3DTFM(u, dm, surface, normals, materialModel, materialProps);
        
        disp(['Calculating metrics of phase ',num2str(p),', cell ',num2str(s)]);
        
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