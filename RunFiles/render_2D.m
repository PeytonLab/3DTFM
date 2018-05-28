%% Visualization
mkdir('Figures'); %Create Figures folder
addpath('Figures');

dir_list = dir; %Direction setup
dir_el = numel(dir)-4;

%Scaling constants
sc = 6.2; %[px/um]
dz = 0.34; %[um/z]
dm = 8; %spacing of displacement grid

load('TFM.mat');

%% 1) Cross-Sectional Displacement and Traction Field (Animation)
gif = 1; %gif creation on (1) or off (0)
fr = 3; %frames per second

for p = 1%:length(TFM) %Phase
    mkdir('Figures',dir_list(p+2).name) %Create condition folder in figures
    cd(dir_list(p+2).name)
    cells_list = dir;
    cells_el = numel(cells_list)-2;   
    for s = 3%:cells_el %Cell
        mkdir(['..\Figures\',dir_list(p+2).name,'\'],cells_list(s+2).name)
        cd([cells_list(s+2).name,'\Volumes']);
        vol_list = dir('vol*'); %directory of volumetric matrices
        tmax = size(vol_list,1)-1; %maximum time points       
        %% Determine cross section to view
        ZIdx = zeros(1,tmax);  
        UIdx = zeros(1,tmax);
        for t = 1:tmax %Determine most intense slice at each time point            
            load(vol_list(t+1).name)
            I = vol{2}; %EGFP channel
            
            %Size matching       
            endsize = (size(TFM{p}{s}{t}{1,4})-1)*dm; %Desired size
            outsize = size(I); %Size of image 
            diffsize = outsize - endsize;      
            padL = ceil(diffsize./2);
            padR = floor(diffsize./2);      
            I = I(1+padL(1):outsize(1)-padR(1),1+padL(2):outsize(2)-padR(2),1+padL(3):outsize(3)-padR(3));
            
            %Find cross section at each time point
            zmax = size(I,3);
            umax = size(TFM{p}{s}{t}{1,4},3);          
            Zint = zeros(zmax,1);
            Uint = zeros(umax,1);            
            for z = 1:zmax
                Zint(z) = sum(sum(I(:,:,z))); %intensity of each slice      
            end
            for z = 1:umax
                Uint(z) = sum(sum(TFM{p}{s}{t}{1,4}(:,:,z)));
            end
            [~,Idx1] = max(Zint); %max intensity slice for time point
            [~,Idx2] = max(Uint);        
            ZIdx(t) = Idx1;
            UIdx(t) = Idx2;
        end   
        UMax = ceil((mean(ZIdx)/dm+1 + mean(UIdx))/2); %average of max intensity and displacement slice, in DVC coordinates
        ZMax = ceil((UMax-1)*dm); %Z slice coordinates
        
        %% Construct cell boundary
        B_cell = cell(tmax,1); %Boundary for each time point
        IB = cell(1,tmax); %Cell slice image for each time point
        for t = 1:tmax
            load(vol_list(t+1).name)
            I = vol{2}; %EGFP channel
            
            %Crop
            endsize = (size(TFM{p}{s}{t}{1,4})-1)*dm; %Desired size
            outsize = size(I); %Size of image 
            diffsize = outsize - endsize;      
            padL = ceil(diffsize./2);
            padR = floor(diffsize./2);      
            I = I(1+padL(1):outsize(1)-padR(1),1+padL(2):outsize(2)-padR(2),1+padL(3):outsize(3)-padR(3));
                        
            IB{t} = I(:,:,ZMax); %Cell slice
            
            for j = 1:3
                IB{t} = medfilt2(IB{t}); %Median filter, thrice
            end
            
            BW = im2bw(IB{t},graythresh(IB{t})); %Binarize
            IB{t} = BW;
            
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
        
        %% Plot index
        sizeI = size(I);   
        plotIdx = cell(1,3);      
        
        for i = 1:2, plotIdx{i} = 1:dm:sizeI(i)+1;
            plotIdx{i} = plotIdx{i}./sc; 
        end
        
        plotIdx{3} = 1:dm:sizeI(3)+1;
        plotIdx{3} = plotIdx{3}.*dz;       
        
        %% Color map bounds          
        for t = 1:length(TFM{p}{s})
            cmin1(t) = min(min(TFM{p}{s}{t}{1,4}(:,:,UMax)));  %Min
            cmax1(t) = max(max(TFM{p}{s}{t}{1,4}(:,:,UMax)));  %Max
            cmin2(t) = min(min(TFM{p}{s}{t}{2,4}));
            cmax2(t) = max(max(TFM{p}{s}{t}{2,4}));
        end
        cmin1 = floor(min(cmin1))/sc; %um/hr
        cmax1 = ceil(max(cmax1))/sc; %um/hr
        cmin2 = floor(min(cmin2))*sc^2/3600; %Pa/s      
        cmax2 = ceil(max(cmax2))*sc^2/3600; %Pa/s
                
        for t = 4%1:length(TFM{p}{s})            
            
            %% Displacement field cross section
            if t == 1 && s == 1 && p == 1
               figure
            end
            
            figure(1);
            ax1 = gca;
            [~,h] = contourf(ax1, plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{1,4}(:,:,UMax)./sc, 50);
            c = colorbar;
            caxis([cmin1,cmax1]);
            set(h,'linestyle','none'); axis image;
            
			%Plot Boundary
            if ~isempty(B_cell{t})                   
                hold on
                plot(B_cell{t}(:,2)./sc,B_cell{t}(:,1)./sc,'Color',[0 1 0],'LineWidth',1.75) %Green boundary
                hold off              
                
                %Mask cell projection onto plot
                M = sizeI(1:2); %Green colormap
                green = cat(3, zeros(M), ones(M), zeros(M));                
                
                hold on
                x0 = [plotIdx{1}(1) plotIdx{1}(end)];
                y0 = [plotIdx{2}(1) plotIdx{2}(end)];
                g = image(x0,y0,green); %Cast green     
                set(g, 'AlphaData',IB{t}*0.3) %Mask cell image
                hold off
            end
            
            %Displacement arrows
            W = cell(1,2);
            X = cell(1,2);
            k = 2; %cropping passes
            for i = 1:2              
                W{i} = TFM{p}{s}{t}{1,i}(:,:,UMax)./sc;
                A = W{i};
                B = plotIdx{i};
                for j = 1:k
                    A = A(2:2:size(A,1)-1,:); %Arrows
                    A = A(:,2:2:size(A,2)-1);
                    B = B(2:2:end-1); %Plot index
                end
                W{i} = A;
                X{i} = B;   
            end            
            u = W{1}; v = W{2};
            x = X{1}; y = X{2}; 
            [x,y] = meshgrid(x,y);           
            
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            
            %Colormap
            cmap1 = [0.0416666679084301;0.0401785708963871;0.0386904776096344;0.0372023805975914;0.0357142873108387;0.0342261902987957;0.0327380970120430;0.0312500000000000;0.0297619048506022;0.0282738097012043;0.0267857145518065;0.0252976194024086;0.0238095242530108;0.0223214291036129;0.0208333339542151;0.0193452388048172;0.0178571436554194;0.0163690485060215;0.0148809524253011;0.0133928572759032;0.0119047621265054;0.0104166669771075;0.00892857182770968;0.00744047621265054;0.00595238106325269;0.00446428591385484;0.00297619053162634;0.00148809526581317;0;0.0694444403052330;0.138888880610466;0.208333328366280;0.277777761220932;0.347222208976746;0.416666656732559;0.486111104488373;0.555555522441864;0.625000000000000;0.694444417953491;0.763888895511627;0.833333313465118;0.840579688549042;0.847826063632965;0.855072438716888;0.862318813800812;0.869565188884735;0.876811563968658;0.884057939052582;0.891304314136505;0.898550689220429;0.905797064304352;0.913043439388275;0.920289874076843;0.927536249160767;0.934782624244690;0.942028999328613;0.949275374412537;0.956521749496460;0.963768124580383;0.971014499664307;0.978260874748230;0.985507249832153;0.992753624916077;1];
            cmap2 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0434782616794109;0.0869565233588219;0.130434781312943;0.173913046717644;0.217391297221184;0.260869562625885;0.304347813129425;0.347826093435288;0.391304343938828;0.434782594442368;0.478260874748230;0.521739125251770;0.565217375755310;0.608695626258850;0.652173936367035;0.695652186870575;0.739130437374115;0.782608687877655;0.826086938381195;0.869565188884735;0.913043498992920;0.956521749496460;1];
            cmap3 = [0;0.0357142873108387;0.0714285746216774;0.107142858207226;0.142857149243355;0.178571432828903;0.214285716414452;0.250000000000000;0.285714298486710;0.321428567171097;0.357142865657806;0.392857134342194;0.428571432828903;0.464285701513290;0.500000000000000;0.535714268684387;0.571428596973419;0.607142865657806;0.642857134342194;0.678571403026581;0.714285731315613;0.750000000000000;0.785714268684387;0.821428596973419;0.857142865657806;0.892857134342194;0.928571403026581;0.964285731315613;1;0.930555582046509;0.861111104488373;0.791666686534882;0.722222208976746;0.652777791023254;0.583333313465118;0.513888895511627;0.444444447755814;0.375000000000000;0.305555552244186;0.236111119389534;0.166666671633720;0.159420296549797;0.152173921465874;0.144927546381950;0.137681156396866;0.130434781312943;0.123188406229019;0.115942031145096;0.108695656061172;0.101449280977249;0.0942028984427452;0.0869565233588219;0.0797101482748985;0.0724637731909752;0.0652173906564713;0.0579710155725479;0.0507246404886246;0.0434782616794109;0.0362318865954876;0.0289855077862740;0.0217391308397055;0.0144927538931370;0.00724637694656849;0];           
            cmap = [cmap1,cmap2,cmap3];
            
            colormap(ax1,cmap); %set colormap
            set(gca,'Ydir','reverse','FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold');
            title('Velocity Field','fontsize',13);
            ylabel(c,'Bead velocity (\mum/hr)','interpreter','Tex','fontsize',12);
            xlabel('\mum'); ylabel('\mum');
            axis([1 sizeI(2) 1 sizeI(1)]./sc);
            set(gcf,'color','white');
            box off
                                  
            %% Traction Forces cross section
            if t == 1 && s == 1 && p == 1
               figure
            end
            figure(2);
            ax2 = gca;
            [~,h] = contourf(ax2, plotIdx{1}, plotIdx{2}, TFM{p}{s}{t}{2,4}*sc^2/3600,50); 
            c = colorbar;
            caxis([cmin2,cmax2]);
            set(h,'linestyle','none'); axis image;
            
            if ~isempty(B_cell{t})
                %Plot boundary
                hold on
                plot(B_cell{t}(:,2)./sc,B_cell{t}(:,1)./sc,'color',[0 1 0],'LineWidth',1.75);
                hold off

                %Plot cell
                hold on
                g = image(x0,y0,green);       
                set(g,'AlphaData',IB{t}*0.3)
                hold off
            end
            
            %Plot arrows
            hold on
            quiverc(x,y,u,v)
            hold off
            
            %Image parameters
            colormap(ax2,cmap)
            set(gca,'Ydir','reverse','FontName','Arial','linewidth',0.75,'fontsize',10,'fontweight','bold')
            title('Traction Field','fontsize',13)
            ylabel(c,'Traction force (Pa/s)','fontsize',12)
            xlabel('\mum'); ylabel('\mum')
            axis([1 sizeI(2) 1 sizeI(1)]./sc)
            box off
            set(gcf,'color', [1 1 1]);
            
            
            %% TR map
            load(vol_list(t+1).name)
            R = cell(1,3);
            for rr = 1:3
                I = vol{rr}; %EGFP channel            
                
                %Crop
                endsize = (size(TFM{p}{s}{t}{1,4})-1)*dm; %Desired size
                outsize = size(I); %Size of image 
                diffsize = outsize - endsize;      
                padL = ceil(diffsize./2);
                padR = floor(diffsize./2);      
                I = I(1+padL(1):outsize(1)-padR(1),1+padL(2):outsize(2)-padR(2),1+padL(3):outsize(3)-padR(3));
                
                R{rr} = I(:,:,ZMax);
            end
            
            Q = R{3}-0.8*R{2}-0.2*R{1};  %TR signal equalizer      
            
                       
            figure
            D = smooth2a(double(Q),100,100);
%           D = smooth2a(double(R{1}),200,200);
            colormap(jet)
            imagesc(D)
             
            for i = 1:2
                TRidx{i} = 1:size(D,i);
            end
            
            [DX,DY] = gradient(D,30,30);
            DD = {DX, DY};
            
            %Displacement arrows
            W = cell(1,2);
            X = cell(1,2);
            k = 5; %croping passes
            for i = 1:2              
                W{i} = DD{i};
                A = W{i};
                B = TRidx{i};
                for j = 1:k
                    A = A(2:2:size(A,1)-1,:); %Arrows
                    A = A(:,2:2:size(A,2)-1);
                    B = B(2:2:end-1); %Plot index
                end
                W{i} = A;
                X{i} = B;   
            end            
            u = W{1}; v = W{2};
            x = X{1}; y = X{2}; 
            [x,y] = meshgrid(x,y);           
            
            %Plot arrows
            hold on
            quiver(x,y,u,v)
            hold off
            axis image
            box off
            set(gcf,'color', [1 1 1]);
                                   
            %Save .gif and .tif
            
			name = [outputdir,'TR_0',num2str(t),'.tif'];
            saveas(gcf,name)  
            close(3)
			
            outputdir = ['../../../Figures/',dir_list(p+2).name,'/',cells_list(s+2).name,'/'];            
            
			for f = 1:2               
                figure(f)
                %Save gif
                if gif == 1
                    %Capture frame
                    frame = getframe(gcf);
                    im{t} = frame2im(frame);
                    [G,map] = rgb2ind(im{t},256);         
                    
                    if f == 1
                       filename = [outputdir,'dis.gif'];
                      else
                       filename = [outputdir,'trac.gif'];
                    end
                    
                    if t == 1
                       imwrite(G,map,filename,'gif', 'LoopCount',inf,'DelayTime',1/fr);
                      else
                       imwrite(G,map,filename,'gif','WriteMode','append','DelayTime',1/fr);
                    end
                end
                
                %Save tif
                if f == 1
                   name = [outputdir,'dis0',num2str(t),'.tif'];
                else
                   name = [outputdir,'trac0',num2str(t),'.tif'];        
                end
                saveas(gcf,name)                    
           
            end   
        end
        
    cd ('..\..')
    end
cd('..')
end