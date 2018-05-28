%% 3) 3D Cell Reconstruction in Displacement Quiver & Traction Map Projection (+ TR MAP)

%DJP 5/27/18 :')

load TFM.mat
for p = 11%:length(TFM) %Phase
    mkdir('Figures',dir_list(p+2).name) %Create condition folder in figures
    
	cd(dir_list(p+2).name)
    cells_list = dir;
    cells_el = numel(cells_list)-2;    
    
    for s = 9%:cells_el %Cell
        mkdir(['..\Figures\',dir_list(p+2).name,'\'],cells_list(s+2).name)
        
		cd([cells_list(s+2).name,'\Volumes']);
        vol_list = dir('vol*'); %directory of volumetric matrices
        tmax = size(vol_list,1)-1; %maximum time points     
        
        for t = 1%:length(tmax)       
            
            load(vol_list(t).name)
            
            %Reconstruct cell
            dim = size(vol{2});           
            v = vol{2};                           
            endsize = (size(TFM{p}{s}{t}{1,4})-1)*dm; %Desired size
            outsize = dim; %Size of image 
            diffsize = outsize - endsize;      
            padL = ceil(diffsize./2);
            padR = floor(diffsize./2);      
            v = v(1+padL(1):outsize(1)-padR(1),1+padL(2):outsize(2)-padR(2),1+padL(3):outsize(3)-padR(3));       
            
            for z = 1:size(v,3)
                I = v(:,:,z);
                for j = 1:3
                    I = medfilt2(I); %Median filter, thrice
                end                
                BW = im2bw(I,0.2); %Binarize
                C(:,:,z) = BW;               
            end
            dim = size(C);
            
            %% Plot index 
            plotIdx = cell(1,3); %displacement dimensions     
            for i = 1:2, plotIdx{i} = 1:dm:dim(i)+1;
                plotIdx{i} = plotIdx{i}./sc; 
            end
            plotIdx{3} = 1:dm:dim(3)+1;
            plotIdx{3} = plotIdx{3}.*dz;     
            
            %% Cell Reconstruction
            [X,Y,Z] = meshgrid((1:dim(1))./sc,(1:dim(2))./sc,(1:dim(3)).*dz); %cell dimensions
            [F,V] = MarchingCubes(X,Y,Z,C,0.8);
            FV = patch('Faces',F,'Vertices',V);
            FV2 = smoothpatch(FV,0,5,1,1);
            FV2.FaceColor = 'r';
            FV2.EdgeAlpha = 0;
            FV2.FaceLighting = 'gouraud';
            axis on
            daspect([1 1 1]) %um,um,um
            set(gcf,'color','c')  
            
            %% Displacement Arrows
            %Displacement arrow crop
            clear u v w
            for z = 1:size(plotIdx{3},2)
                uvw = cell(1,3);
                xyz = cell(1,3);
                k = 3; %croping passes
                for i = 1:3              
                    if i == 3
                       uvw{i} = TFM{p}{s}{t}{1,i}(:,:,z).*dz;
                    else
                       uvw{i} = TFM{p}{s}{t}{1,i}(:,:,z)./sc;
                    end
                    A = uvw{i};
                    B = plotIdx{i};
                    for j = 1:k
                        A = A(2:2:size(A,1)-1,:); %Arrows
                        A = A(:,2:2:size(A,2)-1);
                        B = B(2:2:end-1); %Plot index
                    end
                    uvw{i} = A;
                    xyz{i} = B;   
                end            
                u(:,:,z) = double(uvw{1}); v(:,:,z) = double(uvw{2}); w(:,:,z) = double(uvw{3});
            end
            
            x = double(xyz{1}); y = double(xyz{2}); z = double(plotIdx{3});
            [px,py,pz] = meshgrid(x,y,z);

            hold on
            
            q = quiver3(px,py,pz,u,v,w);            
            q.Color = [0 0.5 1];
            q.LineWidth = 1;
            q.AutoScaleFactor = 3;
            
            %% Traction Map
            [~,h] = contourf(plotIdx{1},plotIdx{2},TFM{p}{s}{t}{2,4}.*sc^2/3600,50);
            set(h,'linestyle','none');
            colormap(jet)
            c = colorbar;
            
            axis([10 70 10 70 0 30])

        end
    end
end