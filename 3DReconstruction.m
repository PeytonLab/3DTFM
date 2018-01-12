%% Cell Surface Rendering

%Cell Reconstruction
dim = size(vol{2});
V = cell(dim(3),1);
for z = 1:size(vol{2},3)   
    I = vol{2}(:,:,z);
    for j = 1:3
        I = medfilt2(I); %Median filter, thrice
    end    
    BW = im2bw(I,0.15); %Binarize
    [B,L,N,~] = bwboundaries(BW,'noholes'); %Boundaries 
    M = zeros(dim(1),dim(2));    
    if N > 0 %Cell boundary
        B_list = zeros(N,1);
        for b = 1:N
            B_list(b) = size(B{b},1);
        end
        [~,Idx] = max(B_list);      
        for x = 1:size(B{Idx},1)
            M(B{Idx}(x,1),B{Idx}(x,2)) = 1;
        end
        M = double(L == Idx);
    else
        disp('No boundary')
    end
    if z == 1
       V = M;
    else
       V = cat(3,V,M);
    end
end

% figure; %Image Scan
% for c = 2
%     for z = 1:size(vol{c},3)
%         imagesc(vol{c}(:,:,z));
%         pause(0.1);
%     end
% end
% 
% figure %Boundary Scan
% for z = 1:size(vol{c},3)
%     imagesc(V(:,:,z));
%     pause(0.05);
% end

figure
isosurface(V,0.8); 