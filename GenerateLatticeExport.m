function Geometry = GenerateLatticeExport(lattice_settings)
% Generate lattices in a beam element format for .3MF (Diameters defined at nodes)
% Input is a structure containing lattice parameters - input empty structure for example/defaults
% Output is a structure inlcuding lattice node, element, diameter
% definitions

%% Defaults
% Define base mesh & unit cell dimensions for lattice  
Opts.Shape='Square';      % 'Square' or 'Cylinder';
Opts.ZDir='Axial';        % If cylinder - orientation of vertical struts ["Radial","Hoop","Axial"] 
Opts.Dim=[1 1 3];         % Number of unit cells in [x,y,z] . Cylinders - [radial,axial,inner-diameter dimension]
Opts.Size=[5 5 5];     % Dimensions of unit cell [x,y,z] or [radial hoop(approx) axial]

% Z Graded aspect ratio - define via one:
Opts.AR_gradient=[];     % Aspect ratio of each layer of unit cells
Opts.layer_heights=[ ];   % heights of each layer

% Define mesh structure for custom base mesh - check output of GIBBON
% functions e.g hexMeshBox - output type 2.
Opts.meshStruct = [];

% Unit cell type - "BCC","BCCZ","FCC","FCCZ","Octet" or "Edge","Dual","tetDual","tetEdge"
Opts.lattice_type="BCC";    % Unit cell type - one of above
Opts.el_per_strut=4;        % elements per strut 

% Diameter definition - including graded & taper lattices
Opts.nominal_diam = 1;  % Strut diameter - nominal diam if taper defined
Opts.taper=[1 0.7 1];    % Strut taper definition - e.g [1 0.7 1] for 0.7 taper
Opts.gradient =[1 0.5];       % Z-Graded strut thicknes (Makes nominal_diam redundant). E.g [1 0.5] 

% Material Transition Plane - (Not relevant but empty variable needed to avoid errors on ViewLattice)
Opts.MaterialTransition = []; % [Normal ; Point In plane]
%% Change parameters accoring to Input
Fields=fieldnames(lattice_settings);
for i=1:numel(Fields)
    Opts.(Fields{i})=lattice_settings.(Fields{i});
end

%% Create Lattice - Wirelattice - struts only
% Aspect ratio gradient - can be defined via aspect ratio or height of each
% layer of unit cells. - Opts.AR_gradient OR Opts.layer_heights

if isempty(Opts.AR_gradient)==0 ||  isempty(Opts.layer_heights)==0 %Angle Gradient
    % First, Generate base mesh.    Calculate heights for each layer
    if isempty(Opts.AR_gradient)==0         %Based on aspect ratio inputs
        Opts.AR_gradient_factors=interp1(1:numel(Opts.AR_gradient),Opts.AR_gradient,...
                                linspace(1,numel(Opts.AR_gradient),Opts.Dim(3)));  
        layer_heights=Opts.AR_gradient_factors*Opts.Size(3);
    else                                    % Based on layer height input
        if length(Opts.layer_heights) == Opts.Dim(3)
            layer_heights = Opts.layer_heights;
        else    % Interpolate if not enough values inputted
            layer_heights=interp1(1:numel(Opts.layer_heights),Opts.layer_heights,...
                                linspace(1,numel(Opts.layer_heights),Opts.Dim(3)));
        end
    end

    for i=1:Opts.Dim(3)     %Create lattice for this layer
        layer = hexMeshBox([Opts.Size(1:2).*Opts.Dim(1:2),layer_heights(i)/1],[Opts.Dim(1:2) 1]);
        %translate to correct height
        if i==2  
            layer.V(:,3)=layer.V(:,3) + layer_heights(1)/2 + layer_heights(i)/2;
        elseif i>2
            layer.V(:,3)=layer.V(:,3) + layer_heights(1)/2 + layer_heights(i)/2 + sum(layer_heights(2:(i-1)));
        end
        %Correct indices in element arrays
        if i>1
            layer.E=layer.E+size(V,1);
            V(size(V,1)+1 : size(V,1)+size(layer.V,1),:)=layer.V;
            E(size(E,1)+1 : size(E,1)+size(layer.E,1),:)=layer.E;
        elseif i==1
            V=layer.V;
            E=layer.E;
        end
    end

    % Merge Duplicate nodes
    [E,V]=mergeVertices(E,V);    
    meshStruct.elements = E;
    meshStruct.nodes = V;

    % Generate lattice
    [N,S,fill_vol] = WireLattice(Opts.lattice_type,'InputMesh',meshStruct,'ZDir',Opts.ZDir);

elseif isempty(Opts.meshStruct) == 0

    [N,S,fill_vol] = WireLattice(Opts.lattice_type,'InputMesh',Opts.meshStruct);

else
    [N,S,fill_vol] = WireLattice(Opts.lattice_type,Opts.Dim,Opts.Size,Opts.Shape,'ZDir',Opts.ZDir);
end
%% Convert to Elements
[V,E,NormalV,S_E] = struts2beams(N,S,Opts.el_per_strut);%Opts.el_per_strut,Opts.element_type);

%% Calculate Strut Diameters
% For export to 3MF, diameter defined at each node instead of for each
% element

% Get diameters based on gradient/diameter definition
if isempty(Opts.gradient)==0        %If Opts.gradient is enabled
    % get mid points of each element
    El_endNodes=[V(E(:,1),:)  ;  V(E(:,2),:) ]; 

    % Get unique heights 
    [unique_heights,i_unique,i_original] = unique(El_endNodes(:,3));

    % Heights for each term in input gradient definition
    gradient_heights = linspace(min(unique_heights),max(unique_heights),numel(Opts.gradient));
   
    % Diameters for elements at each height
    unique_diameters = interp1(gradient_heights,Opts.gradient,unique_heights,'makima');
    
    % Create full diameter array based on original position in element array
    element_diams = unique_diameters(i_original);

    % Reshape so there is two diam values on each row
    element_diams = reshape(element_diams,size(E,1),2);

else
    element_diams=ones(size(E,1),2)*Opts.nominal_diam;
    unique_diameters = [];
end

% Get diameter multipliers based on Taper factor
if isempty(Opts.taper)==0       %If Opts.taper is enabled
    %Interpolate taper definition to get diameter multiplier for every element in a strut
    taper_factors = interp1(1:numel(Opts.taper),Opts.taper,linspace(1,numel(Opts.taper),Opts.el_per_strut+1)); 
    
    % structure so there are two diameters for each element
    taper_factors_2 = [taper_factors(1:end-1)' taper_factors(2:end)'];

    %Repeat for every strut in lattice
    taper_multipliers=repmat(taper_factors_2,[size(S,1) 1]);
else
    taper_multipliers=ones(size(E,1),2);
end

% Combine gradient and taper
element_diams=element_diams.*taper_multipliers;


%% Store Data to re-create lattice
Geometry.lattice_params=Opts;
Geometry.Diameters=element_diams;
if isempty(Opts.gradient)==0 %If Gradient exists
    Geometry.gradient_factors=[unique_diameters unique_heights];
end
Geometry.V=V;
Geometry.element_normals=NormalV;
Geometry.E=E;
Geometry.N=N;
Geometry.S=S;
Geometry.S_E=S_E;
Geometry.material_ID=[];

end

%%
% Lattice Inverse Design & Optimisation Tool 
% 06/12/2023 - Brian McDonnell - University of Galway
% GNU AFFERO GENERAL PUBLIC LICENSE - See LICENSE file details