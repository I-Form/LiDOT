function [V,E,NormalV,S_E] = struts2beams(N,S,varargin)
% [V,E,Normal Vertices,Strut-element array] = struts2beams(N,S,ElementsPerStrut,node_el_lengths)
% Converts a lattice defined using lattice nodes and struts into elements
% with normals for the abaqus beam section definition
% Element type: 1 for Linear (B31), 2 for Quadratic (B32)
% Node_el_lenghts is the length of each element at the intersections -
% GenerateLattice for example input

% Process inputs and set defaults
switch nargin
    case 2
        ElementsPerStrut=6;
        elementtype=1;
        node_el_lengths=[];
    case 3
        ElementsPerStrut=varargin{1};
        elementtype=1;
        node_el_lengths=[];
    case 4
        ElementsPerStrut=varargin{1};
        elementtype=varargin{2};
        node_el_lengths=[];
    case 5
       ElementsPerStrut=varargin{1};
       elementtype=varargin{2};
       node_el_lengths=varargin{3};

end

%% Subdivide Struts for multiple elements per strut
strut_vectors=N(S(:,2),:) - N(S(:,1),:);  %gets vector of direction of each strut
strut_lengths=sqrt( (strut_vectors(:,1)).^2 + (strut_vectors(:,2)).^2 + (strut_vectors(:,3)).^2);   %Get length of each strut

% Determine locations of new nodes along each strut 
newnodes = 0 : 1/ElementsPerStrut : 1 ;
node_locations = newnodes(2:end-1);
NodesperStrut=numel(node_locations);     %Nodes to be added per strut

%Cycle through nodes to be added
if isempty(node_el_lengths)==0  % If increased diameter at nodes - elements created with corresponding length
    % Vector addition to nodes for first and last element in each strut
    % First element in strut
    node_el_1= N(S(:,1),:) + ((node_el_lengths(:,1)./strut_lengths).*strut_vectors);   
    % Last element in strut 
    node_el_2= N(S(:,1),:) + ((1 - node_el_lengths(:,2)./strut_lengths).*strut_vectors);  

    % New strut vector to account for new length of middle section
    mid_vectors=node_el_2 - node_el_1;  %gets vector of direction of each strut

    mid_nodes=zeros(size(S,1)*(NodesperStrut-2),3);      %Empy mid node array
    % Get mid-nodes through vector addition. Sub-dividing mid_vectors
    for i = 1:NodesperStrut-2
        mid_nodes((1+(i-1)*size(S,1)):(i)*size(S,1),:)= node_el_1 + (i/(ElementsPerStrut-2))*mid_vectors;   
    end
    mid_nodes = [node_el_1; mid_nodes; node_el_2];

else
    mid_nodes=zeros(size(S,1)*(NodesperStrut),3);      %Empy mid node array
    % Get mid-nodes through vector addition. Sub-dividing strut_vectors
    for i = 1:NodesperStrut
        mid_nodes((1+(i-1)*size(S,1)):(i)*size(S,1),:)= N(S(:,1),:) + (i/(ElementsPerStrut))*strut_vectors;   
    end
end

% Create element IDs and arrays for elements in each strut and their
% vertices.
MidNodeIDs = length(N) + reshape(1:length(mid_nodes),size(S,1),NodesperStrut);  % IDs of new nodes
S=[S(:,1) MidNodeIDs S(:,2)];      %add to strut array
V=[N ; mid_nodes];                  %add new nodes to node array

%% Get element array from strut array 
% List all nodes
node_list=S';  
node_list=node_list(:); 

% Order adjacent nodes into elements
E=[node_list(1:end-1) node_list(2:end)];    

% Remove elements produced from adjacent nodes not within the same strut
E_true=logical(ones(length(E),1));                      
E_true(size(S,2):size(S,2):end)=false; 
Es=E(E_true,:);  %list of elements

%get strut array indexing each element in each strut
S_E= (reshape(1:size(Es,1),ElementsPerStrut,size(S,1)))';

%% Get middle nodes if quadratic elements used
if elementtype==2    %Gets mid-nodes if quadratic elements selected  
    qMidNodes=0.5*( V(Es(:,1),:)  +  V(Es(:,2),:) ); %new edge mid-points
    qMidNodeIDs= ( length(V) + (1:length(qMidNodes)) )';
    V = [V; qMidNodes]; %Join point sets
    Es= [Es(:,1) qMidNodeIDs Es(:,2)]   ;    %new quadratic element set
end

%% Find normals for Abaqus beam section definition
[a,~]=vectorOrthogonalPair(strut_vectors);

NormalV=a+N(S(:,1),:);   %Adding vector to nodes to get a normal node
NormalInd=zeros(size(Es,1),1);

for j=1:ElementsPerStrut
    NormalInd(j:ElementsPerStrut:end)= length(V) + (1:length(NormalV));  %Indices
end
E = [Es NormalInd];     %Adds indices of normal nodes to element array      


%Plot lattice and normal nodes
% figure
% plot3(NormalV(:,1),NormalV(:,2),NormalV(:,3),'*y')
% axisGeom(gca,15);
% camlight headlight;
% drawnow;
% hold on
%     Xstruts= [V(E(:,1),1) V(E(:,2),1)]';
%     Ystruts= [V(E(:,1),2) V(E(:,2),2)]';
%     Zstruts= [V(E(:,1),3) V(E(:,2),3)]';
%     plot3(Xstruts,Ystruts,Zstruts,'-k','LineWidth',3,'Marker','o')
% hold off

end

%%
% Lattice Inverse Design & Optimisation Tool 
% 06/12/2023 - Brian McDonnell - University of Galway
% GNU AFFERO GENERAL PUBLIC LICENSE - See LICENSE file details