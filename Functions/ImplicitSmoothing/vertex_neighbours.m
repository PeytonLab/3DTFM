function Ne=vertex_neighbours(FV)
% This function VERTEX_NEIGHBOURS will search in a face list for all 
% the neigbours of each vertex.
%
% Ne=vertex_neighbours(FV)
%

Ne=vertex_neighbours_double(FV.Faces(:,1),FV.Faces(:,2),FV.Faces(:,3),FV.Vertices(:,1),FV.Vertices(:,2),FV.Vertices(:,3));