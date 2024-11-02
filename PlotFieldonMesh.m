function PlotFieldonMesh(coordinates,nodes, component )
%--------------------------------------------------------------------------
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com                                         |
%          http://sites.google.com/site/kolukulasivasrinivas/             |    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To plot the profile of a component on mesh
% Synopsis :
%           ProfileonMesh(coordinates,nodes,component)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y ] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]      
%           component - The components whose profile to be plotted
%           -----> components = a column vector in the order of node
%                               numbers
%--------------------------------------------------------------------------


nel = length(nodes) ;                  % number of elements
nnode = length(coordinates) ;          % total number of nodes in system
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
Z = zeros(nnel,nel) ;
profile = zeros(nnel,nel) ;
%
for iel=1:nel   
     for i=1:nnel
     nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=coordinates(nd(i),1);    % extract x value of the node
     Y(i,iel)=coordinates(nd(i),2);    % extract y value of the node
     end   
     profile(:,iel) = component(nd') ;         % extract component value of the node 
end
X([3,4],:) = X([4,3],:);
Y([3,4],:) = Y([4,3],:);
profile([3,4],:) = profile([4,3],:);    
% Plotting the FEM mesh and profile of the given component
     % f3 = figure ;
     % f3.Position(3) = f3.Position(4) * Lpara(1) / Lpara(2);
     % set(f3,'name','Postprocessing','numbertitle','off') ;
     % plot(X,Y,'k')
     patch(X,Y,profile, 'LineStyle', 'none', 'FaceAlpha', 0.5)
     % patch(X,Y,profile)
     axis off ;
     axis equal
     % caxis(clim)
     % colorbar
     % F = getframe(gca);
     % writeVideo(v,F);
     % close
     % Colorbar Setting
     % SetColorbar

     % f3 = figure ;
     % f3.Position(3) = f3.Position(4) * Lpara(1) / Lpara(2);
     % set(f3,'name','Postprocessing','numbertitle','off') ;
     % pdeplot(mesh, XYData = gather( component) );
     % axis off ;
     % caxis(clim)
     % colorbar
     % F = getframe(gca);
     % writeVideo(v,F);
     % close





 end

              
         
 
   
     
       
       

