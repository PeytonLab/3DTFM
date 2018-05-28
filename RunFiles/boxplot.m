%% 2) Box plot of population tractions and displacements
load TFM.mat
p = length(TFM); %number of conditions

S_U = cell(1,p); %displacements of cells
S_F = cell(1,p); %traction of cells 

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
			U(k) = mean(abs(TFM{i}{j}{k}{1,4}(:)./sc)); %mean displacement, um/hr
            F(k) = mean(abs(TFM{i}{j}{k}{2,4}(:).*(sc^2/3600))); %traction force, Pa/s		
        end
        
        %Store vectors for each cell in a cell
        if j == 1
		   S_U{i} = U;
		   S_F{i} = F;
        else 			
           S_U{i} = horzcat(S_U{i},U);
		   S_F{i} = horzcat(S_F{i},F);
        end
    end
    
    %If values are zero, replace w/ nan (needed for box plot)
	S_U{i}(S_U{i}==0)=nan;
	S_F{i}(S_F{i}==0)=nan;	
end

%Group into population averages
for i = 1:p
    P_U{i} = nanmean(S_U{i})';  
    P_F{i} = nanmean(S_F{i})';
end

%Convert cells of cell populations to a matrix w/ nan for space
for i = 1:p
    d(i) = length(P_U{i});
end
dmax = max(d); %maximum spaces
PhaseDisplacement = NaN(dmax,p); %fill with empty space
PhaseTraction = NaN(dmax,p);

for i = 1:p
    for	j = 1:length(P_U{i})
		PhaseDisplacement(j,i) = P_U{i}(j); %Fill space with displacements	
    end
    for k = 1:length(P_F{i})
        PhaseTraction(k,i) = P_F{i}(k);
    end
end

%Plot
for i=7%1:p %For each population
	figure
    x = boxplot(S_U{i});
    xlabel('Cell')
    ylabel('Speed (um/hr)')
    str=['Bead Displacements (P',num2str(i),')'];
    title(str)
    set(x,{'linew'},{1})
    set(gcf,'color','white')
    box off

    figure
	x = boxplot(S_F{i});
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
x = boxplot(PhaseDisplacement,'Labels',{'10wt_1mMRGD_25MMP','20wt_0.5mMRGD_0MMP_v1','20wt_0.5mMRGD_0MMP_v2','20wt_0.5mMRGD_50MMP', ...
    '20wt_2mMRGD_0MMP','20wt_2mMRGD_50MMP','20wt_2mMRGD_50MMP_v2','20wt_2mMRGD_50MMP_v3', ...
    '5wt_0.5mMRGD_0MMP','5wt_0.5mMRGD_50MMP','5wt_2mMRGD_0MMP'});
set(x(7,:),'Visible','off')
ylim([0 2])
set(gca,'YTick',0:1:8,'FontName','Arial','linewidth',0.75,'fontsize',8,'fontweight','bold')
ylabel('Speed (um/hr)','fontsize',12)
title('Bead Displacement','fontsize',13)
set(x,{'linew'},{1})
set(gcf,'color','white')
box off

%Population tractions
figure
x = boxplot(PhaseTraction,'Labels',{'10wt_1mMRGD_25MMP','20wt_0.5mMRGD_0MMP_v1','20wt_0.5mMRGD_0MMP_v2','20wt_0.5mMRGD_50MMP', ...
    '20wt_2mMRGD_0MMP','20wt_2mMRGD_50MMP','20wt_2mMRGD_50MMP_v2','20wt_2mMRGD_50MMP_v3', ...
    '5wt_0.5mMRGD_0MMP','5wt_0.5mMRGD_50MMP','5wt_2mMRGD_0MMP'});
set(x(7,:),'Visible','off')
ylim([0 70])
set(gca,'YTick',0:10:100,'FontName','Arial','linewidth',0.75,'fontsize',8,'fontweight','bold')
ylabel('Force (Pa/s)','fontsize',12)
title('Cellular Traction Forces','fontsize',13)
set(x,{'linew'},{1})
set(gcf,'color','white')
box off