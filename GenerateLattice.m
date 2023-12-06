function Geometry = GenerateLattice(lattice_settings)
% Geometry = GenerateLattice(lattice_settings)
% Generate lattices in a beam element format for FEA (Diameter defined per element)
% Input is a structure containing lattice parameters - input empty structure for example/defaults
% Output is a structure containing beam element geometry of the lattice

%%  Default Settings
% Define base mesh & unit cell dimensions for lattice  
Opts.Shape='Square';      % 'Square' or 'Cylinder';
Opts.ZDir='Axial';        % If cylinder - orientation of vertical struts ["Radial","Hoop","Axial"] 
Opts.Dim=[1 1 1];         % Number of unit cells in [x,y,z] . Cylinders - [radial,axial,inner-diameter dimension]
Opts.Size=[10 10 10];     % Dimensions of unit cell [x,y,z] or [radial hoop(approx) axial]

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
Opts.taper=[1 1 1];    % Strut taper definition - e.g [1 0.7 1] for 0.7 taper
Opts.gradient =[ ];       % Z-Graded strut thicknes (Makes nominal_diam redundant). E.g [1 0.5] 

% FE-related
Opts.element_type=1;
Opts.node_el_multiplier = [];
% Material Transition Plane
Opts.MaterialTransition = []; % [Normal ; Point In plane]
%% Change parameters accoring to Input
Fields=fieldnames(lattice_settings);
for i=1:numel(Fields)
    Opts.(Fields{i})=lattice_settings.(Fields{i});
end

%% Create Lattice - Wirelattice - struts only
if isempty(Opts.AR_gradient)==0 ||  isempty(Opts.layer_heights)==0 
    % If gradient of unit cell aspect ratio is defined
    % First, Generate base mesh.    Calculate heights for each layer
    if isempty(Opts.AR_gradient)==0         %Based on aspect ratio inputs
        Opts.AR_gradient_factors=interp1(1:numel(Opts.AR_gradient),Opts.AR_gradient,...
                                linspace(1,numel(Opts.AR_gradient),Opts.Dim(3)));
    
    
        layer_heights=Opts.AR_gradient_factors*Opts.Size(3);
    else                                    % Based on direct input
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
    % If lattice is based on input mesh
    [N,S,fill_vol] = WireLattice(Opts.lattice_type,'InputMesh',Opts.meshStruct);

else
    % No unit cell aspect ratio or input mesh
    [N,S,fill_vol] = WireLattice(Opts.lattice_type,Opts.Dim,Opts.Size,Opts.Shape,'ZDir',Opts.ZDir);
end

%% Convert to Elements

[V,E,NormalV,S_E] = struts2beams(N,S,Opts.el_per_strut,Opts.element_type);%Opts.el_per_strut,Opts.element_type);

%% Calculate Strut Diameters
% Get diameters based on gradient/diameter definition
if isempty(Opts.gradient)==0        %If Opts.gradient is enabled
    % get mid points of each element
    El_MidNodes=0.5*( V(E(:,1),:)  +  V(E(:,2),:) ); 

    % Get unique heights 
    [unique_heights,i_unique,i_original] = unique(El_MidNodes(:,3));

    % Heights for each term in input gradient definition
    gradient_heights = linspace(min(unique_heights),max(unique_heights),numel(Opts.gradient));
   
    % Diameters for elements at each height
    unique_diameters = interp1(gradient_heights,Opts.gradient,unique_heights,'makima');
    
    % Create full diameter array based on original position in element array
    element_diams = unique_diameters(i_original);

else
    element_diams=ones(size(E,1),1)*Opts.nominal_diam;
    unique_diameters = [];
end

% Get diameter multipliers based on Taper factor
if isempty(Opts.taper)==0       %If Opts.taper is enabled
    %Interpolate taper definition to get diameter multiplier for every element in a strut
    taper_factors = interp1(1:numel(Opts.taper),Opts.taper,linspace(1,numel(Opts.taper),Opts.el_per_strut)); 
    
    %Repeat for every strut in lattice
    taper_multipliers=repmat(taper_factors',[size(S,1) 1]);
else
    taper_multipliers=ones(size(E,1),1);
end

% Combine gradient and taper
element_diams=element_diams.*taper_multipliers;


%% Re-mesh after assigning diameters - (Only if node_el_multiplier specified)
% length of elements for node_el_multiplier are dependant on diameter, need
% to be re-meshed after diameter is calculcated.
if isempty(Opts.node_el_multiplier) == 0  % recreate lattice if node lengths on
    %Node Element lengths
    node_el_lengths = [element_diams(S_E(:,1),:) element_diams(S_E(:,end),:)]/2;
    
    % Convert to Elements
    [V,E,NormalV,S_E] = struts2beams(N,S,Opts.el_per_strut,Opts.element_type,node_el_lengths);%Opts.el_per_strut,Opts.element_type);
    
    % Create new diameters with increase at nodes 
    if isempty(Opts.gradient)==0  %If Opts.gradient is enabled
        % get mid points of each element
        El_MidNodes=0.5*( V(E(:,1),:)  +  V(E(:,2),:) ); 
    
        % Get unique heights 
        [unique_heights,i_unique,i_original] = unique(El_MidNodes(:,3));
    
        % Heights for each term in gradient definition
        gradient_heights = linspace(min(unique_heights),max(unique_heights),numel(Opts.gradient));

        % Diameter for each unique height
        unique_diameters = interp1(gradient_heights,Opts.gradient,unique_heights,'makima');
        
        % Diameter for each element
        element_diams = unique_diameters(i_original);
    else
        element_diams=ones(size(E,1),1)*Opts.nominal_diam;
        unique_diameters = [];
    end

    if isempty(Opts.taper)==1      % If no taper defined, create definition to apply node increase
        Opts.taper = [1 1 ];
    end
    %Interpolate values for every element in a strut
    taper_factors = interp1(1:numel(Opts.taper),Opts.taper,linspace(1,numel(Opts.taper),Opts.el_per_strut)); 

    % Applied increase in diam at nodes through taper
    taper_factors(1) = taper_factors(1)*Opts.node_el_multiplier;
    taper_factors(end) = taper_factors(end)*Opts.node_el_multiplier;

    %Repeat for every strut in lattice
    taper_multipliers=repmat(taper_factors',[size(S,1) 1]);
        
    % Combine gradient and taper
    element_diams=element_diams.*taper_multipliers;
end

%% Calculate clearence required at top and bottom of lattice to offset compression platens
% Helps ensure lattice and plates are in contact at beginning of step,
% however does not account for ABAQUS reducing effective contact size of elements

% gets vector of direction of each element
el_vectors=V(E(:,end-1),:) - V(E(:,1),:);  
% Get length of each element (for volume calculation)
el_lengths=sqrt( (el_vectors(:,1)).^2 + (el_vectors(:,2)).^2 + (el_vectors(:,3)).^2);   

% Angle of a strut to xy plane
el_angle = atan2d(el_vectors(1,3),sqrt(el_vectors(1,2)^2 + el_vectors(1,1)^2)) ;

% get max diameter at top and bottom of lattice
if isempty(Opts.gradient)==0 %If Gradient exists
    [~,i_bottom] = min(unique_heights);
    [~,i_top] = max(unique_heights);
    max_d_bot = element_diams(i_unique(i_bottom)) ;
    max_d_top = element_diams(i_unique(i_top)) ;

elseif isempty(Opts.taper)==0 %If only Taper
    max_d_bot = max(Opts.taper) * Opts.nominal_diam ;
    max_d_top =  max_d_bot ;
else % if neither
    max_d_bot = Opts.nominal_diam ;
    max_d_top = Opts.nominal_diam ;
end

% Clearence calculated from the angle of the elements and the diameter
if isempty(Opts.node_el_multiplier) == 0
    clearence = [(max_d_bot/(2*Opts.node_el_multiplier))*cosd(el_angle) (max_d_top/(2*Opts.node_el_multiplier))*cosd(el_angle)];
else
    clearence = [(max_d_bot/2)*cosd(el_angle) (max_d_top/2)*cosd(el_angle)];
end

%% Assign Material ID (if multi material)
if isempty(Opts.MaterialTransition) == 0        %If material tranaitionis enabled
    % initialise variables
    dist_plane = zeros(size(E,1),1);
    material_ID = zeros(size(E,1),1);

    for i = 1:size(Opts.MaterialTransition,1)
        for j = 1 : size(E,1)
            % Material assigned based on normal to the plane and element centre, 
            % if angle is less than 90 deg -> dot product is positive -> element is on 'other' side of plane 
            dist_plane(j) = dot(El_MidNodes(j,:) - Opts.MaterialTransition(i,4:6), Opts.MaterialTransition(i,1:3));
        end
       
        if i == 1
            material_ID( dist_plane <= 0 ) = Opts.MaterialSectionIDs(1);
        end
        
        material_ID( dist_plane > 0 ) = Opts.MaterialSectionIDs(i+1);
    end
    Geometry.material_ID=material_ID;
else
    Geometry.material_ID=ones(size(E,1),1);
end

%% Get volume and relative density
element_vol=sum((el_lengths*pi()/4).*(element_diams.^2),'all') ;
relative_density=element_vol/fill_vol;

%% Store Data to re-create lattice
Geometry.lattice_params=Opts;

Geometry.Diameters=element_diams;
Geometry.relative_density=relative_density;
Geometry.total_element_volume=element_vol;
Geometry.clearence=clearence;
if isempty(Opts.gradient)==0 %If Gradient exists
    Geometry.gradient_factors=[unique_diameters unique_heights];
end
Geometry.V=V;
Geometry.element_normals=NormalV;
Geometry.E=E;
Geometry.N=N;
Geometry.S=S;
Geometry.S_E=S_E;

   

end