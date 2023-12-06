function [StressStrain,AnalysisTime,Energies]=ExplicitBeamCompressionSP(Geometry,FEA_Opts)
% Wites ABAQUS input file and runs job for compression of beam element lattice structure. 
%
% SP - Sanwich panels which are modelled with hex elements and attched to/share nodes with lattice
%
% INPUTS
% Geometry - output of GenerateLattice function
% FEA_Opts - analysis options, see below or OptRun for details.
%
% ---- DEMO ---- 
% run with empty input strtuctures for unit cell demo/default settings
% e.g, run "ExplicitBeamCompressionSP(GenerateLattice(struct),struct)"
% 

Plots=0;

abaqusPath='abq2019';%'/usr/bin/abaqus'; %Abaqus excute command or path
%% DEFAULT Analysis Options
%abaqusInpFileNamePart= ['opt_' num2str(Taper(2)*10) ];   % Defining Abaqus file names
Opts.FileName= 'OptRun';   % Defining Abaqus file names
Opts.appliedStrain=-0.75;   %Applied strain to top plate (in Y direction)

% Step Settings 
Opts.Time=1;    %Step time length
Opts.MassScalingOn=1;
Opts.MSIncSize=1e-6;   %Mass scaling - defined stable increment size
Opts.Num_intervals=200;   %numbr of output intervals
% Contact and BCs
Opts.StrutContact=1;          %Strut contact on (1) / off (0) / 2 defining surfaces - not needed for explicit
Opts.BCsetting=5;            % 4 - Sandwich panels with hex elements instead of rigid plates,  5 - Sandwich panels with symmetry BCs

Opts.PanelMesh_ZSize = [2 2];      % [Panel thickness in mm , No. of Elements through thickness]
% Number of Element sin x-y direction per unit cell - must be multiple of 2 if SquareSymm is used for lattice
Opts.PanelMesh_XYSizeFactor = 2;
%% DEFAULT Material Parameters
Opts.Material.PlasticOn= 1;   %0 for off, 1 for on
Opts.Material.DamageOn=0;

% 17-4 PH (B) No Heat treatment - Data Sheet
Opts.Material.Name='Datasheet 174PH';
Opts.Material.Elastic=[200000 0.3];
Opts.Material.Plastic=[620 0; 1100 0.16];
Opts.Material.Density=7.8e-9;

Opts.FrictionCoeff=0;

%% Change parameters accoring to "Options Input"
Fields=fieldnames(FEA_Opts);
for i=1:numel(Fields)
    Opts.(Fields{i})=FEA_Opts.(Fields{i});
end

%% Define variables for lattice Geometry
V = Geometry.V;
E = Geometry.E;
element_normals = Geometry.element_normals;
S_E = Geometry.S_E;
Diameters = Geometry.Diameters;
clearence = Geometry.clearence;

% Define some variables based on inputs
el_per_strut=Geometry.lattice_params.el_per_strut;  
elementtype=size(E,2) - 2;      
ConExcl=round((el_per_strut/4));           %Number of elements excluded from contact at nodes


%% Path names Control parameters
savePath =fullfile(fileparts(mfilename('fullpath')),'\abq_out');
%savePath=fullfile(defaultFolder,'data','temp');

abaqusInpFileName=fullfile(savePath,[Opts.FileName,'.inp']); %INP file name

%% Get Size of Part
GeomSize(1,:) = max(V,[],1); GeomSize(2,:) = min(V,[],1); 
GeomSize(3,:)=abs(diff(GeomSize,[],1));
BotPHeight = GeomSize(2,3) - clearence(1) ; 
TopPHeight = GeomSize(1,3) + clearence(2) ;
% TopPHeight=GeomSize(1,3)+0.5*(max(Diameters(ETop)));  %Set bottom plate to minumum Y
% BotPHeight=GeomSize(2,3)-0.5*(max(Diameters(EBotttom)));  %Set bottom plate to max Y
displacementMagnitude=GeomSize(3,3)*Opts.appliedStrain; 

% Top and bottom of lattice
Vind_lat=1:size(V,1);
VBottom_lat=Vind_lat(V(:,3)==min(V(:,3)));
VTop_lat=Vind_lat(V(:,3)==max(V(:,3))); 

% Adds normal nodes (used for beam section definition) to node array
V_lat = [V; element_normals];  

normal_IDs = unique(E(:,end));

%E_Vectors=V(E(:,end),:) - V(E(:,1),:);  %gets vector of direction of each element
%E_midpoints= (V(E(:,end),:) + V(E(:,1),:))./2;  %gets midpoint of each element
%E_Lengths=sqrt( (E_Vectors(:,1)).^2 + (E_Vectors(:,2)).^2 + (E_Vectors(:,3)).^2);   %Get length of each element


%% Create Sandwich Panels
boxDim=[Geometry.lattice_params.Size(1:2).*Geometry.lattice_params.Dim(1:2)...
        Opts.PanelMesh_ZSize(1)];
boxEl=[Geometry.lattice_params.Dim(1:2).*Opts.PanelMesh_XYSizeFactor...
        Opts.PanelMesh_ZSize(2)];
[meshStruct]=hexMeshBox(boxDim,boxEl,2);

E_p=meshStruct.elements;
V_p=meshStruct.nodes;

% Create copies of panels top and bottom of lattice
% X-Y displacement based on min x and y positions
p_disp =  min(V(:,1:2)) - min(V_p(:,1:2));

% Z displacement based on distance from plate surfaces to lattice ends
pT_disp =  [p_disp ( max(V(VTop_lat,3)) - min(V_p(:,3)) )];
pB_disp =  [p_disp ( max(V(VBottom_lat,3)) - max(V_p(:,3)) )];

V_pT = V_p + pT_disp;
V_pB = V_p + pB_disp;

% join sets
[Ep,Vp] = joinElementSets({E_p,E_p},{V_pT,V_pB}) ;

[E_j,V_join] = joinElementSets({E,Ep},{V_lat,Vp});

[E_j,V_join,~,i_old]  = mergeVertices(E_j,V_join);

E_beam = E_j{1};
E_hex = E_j{2};

normal_IDs = i_old(normal_IDs);
%% Top/Bottom Surfaces, and Corner Nodes for BCs
% % Find top/bottom and sides
Vind=1:size(V_join,1);
V_BCcheck = V_join;
V_BCcheck(normal_IDs,:) = repmat(mean(V_join),numel(normal_IDs),1); % Make sure normal points for beam section definitions don't effevt

VBottom=Vind(V_BCcheck(:,3)==min(V_BCcheck(:,3)));
VTop=Vind(V_BCcheck(:,3)==max(V_BCcheck(:,3))); 

bcXFace= Vind(V_BCcheck(:,1)==min(V_BCcheck(:,1))); 
bcYFace= Vind(V_BCcheck(:,2)==min(V_BCcheck(:,2))); 

%% Define Element & Node IDs
beam_element_IDs = (1:size(E_beam,1))';
hex_element_IDs = ( (size(E_beam,1)+1) : (size(E_beam,1)+size(E_hex,1)) )';
element_IDs = ( 1 : (size(E_beam,1)+size(E_hex,1)) )';
nodeIds=(1:1:size(V_join,1))';

num_materials = numel(unique(Geometry.material_ID));

%% Defining the abaqus input structure
% See also |abaqusStructTemplate| and |abaqusStruct2inp| and the abaqus user
% manual

%%--> Heading
abaqus_spec.Heading.COMMENT{1}='Job name: Lattice Compression';
abaqus_spec.Heading.COMMENT{2}='Generated by: GIBBON';

%%--> Preprint
abaqus_spec.Preprint.ATTR.echo='NO';
abaqus_spec.Preprint.ATTR.model='NO';
abaqus_spec.Preprint.ATTR.history='NO';
abaqus_spec.Preprint.ATTR.contact='YES';

%--> Part

% Node
abaqus_spec.Part{1}.COMMENT='This section defines the part geometry in terms of nodes and elements';
abaqus_spec.Part{1}.ATTR.name='Cube';
abaqus_spec.Part{1}.Node={nodeIds,V_join};


if elementtype==2   %Quadratic Elements
    abaqus_spec.Part{1}.Element{1}.ATTR.type='B32';
elseif elementtype==1   %Linear Elements
    abaqus_spec.Part{1}.Element{1}.ATTR.type='B31';
end
abaqus_spec.Part{1}.Element{1}.VAL={beam_element_IDs,E_beam};

abaqus_spec.Part{1}.Element{2}.ATTR.type='C3D8';
abaqus_spec.Part{1}.Element{2}.ATTR.type='C3D8';
abaqus_spec.Part{1}.Element{2}.VAL={hex_element_IDs,E_hex};


% Element sets 
% All Elements
abaqus_spec.Part{1}.Elset{1}.ATTR.elset='AllEl'; 
abaqus_spec.Part{1}.Elset{1}.VAL=element_IDs;

% Elements at lattice nodes (to stiffen nodes)
NodeEls=reshape( S_E(:,[1 end]) , numel(S_E(:,[1 end])),1);
abaqus_spec.Part{1}.Elset{2}.ATTR.elset='NodeEL';
abaqus_spec.Part{1}.Elset{2}.VAL=NodeEls' ;

abaqus_spec.Part{1}.Elset{3}.ATTR.elset='Panel_Elements';
abaqus_spec.Part{1}.Elset{3}.VAL=hex_element_IDs ;



% Set & Section Definitions - For each group of unique diameter & material
for i = 1:num_materials
    mat_i_ID = find(Geometry.material_ID == i);
    mat_i_dia=round(Diameters(mat_i_ID),6);
    [UniqueDiams,~,u_ID]=unique(mat_i_dia,'stable');
    num_dia = numel(UniqueDiams);

    if i == 1
        set_count = 0;
    end

    for j = 1:num_dia
        %diameter_group = find(ub == j);
        diameter_group = mat_i_ID(u_ID == j);
        abaqus_spec.Part{1}.Elset{3 + set_count + j}.ATTR.elset=['Mat_' num2str(i) 'D_' num2str(j)];
        abaqus_spec.Part{1}.Elset{3 + set_count + j}.VAL=diameter_group;

        abaqus_spec.Part{1}.Beam_section{set_count + j}.ATTR.elset=['Mat_' num2str(i) 'D_' num2str(j)];
        abaqus_spec.Part{1}.Beam_section{set_count + j}.ATTR.material=Opts.Material(i).Name;
        abaqus_spec.Part{1}.Beam_section{set_count + j}.ATTR.temperature='GRADIENTS';
        abaqus_spec.Part{1}.Beam_section{set_count + j}.ATTR.section='CIRC';
        abaqus_spec.Part{1}.Beam_section{set_count + j}.VAL{1,1}={UniqueDiams(j)/2};
    end
    
    set_count = set_count + j;

end


abaqus_spec.Part{1}.Solid_section.ATTR.elset='Panel_Elements';
abaqus_spec.Part{1}.Solid_section.ATTR.material=Opts.Material(1).Name;


% %Surfaces - rigid plates
% abaqus_spec.Part{2}.COMMENT='This section defines the part geometry in terms of nodes and elements';
% abaqus_spec.Part{2}.ATTR.name='rigid_plate';
% PlateRefNode=[0 0 0];
% abaqus_spec.Part{2}.Node={1,PlateRefNode};
% 
% abaqus_spec.Part{2}.surface.ATTR.type='cylinder';
% abaqus_spec.Part{2}.surface.ATTR.name='rigid_plate';
% abaqus_spec.Part{2}.surface.VAL{1,1}={{'START'}, [2*GeomSize(3,1) 0]};
% abaqus_spec.Part{2}.surface.VAL{2,1}={{'LINE'},[-(2*GeomSize(3,1)) 0] };
% 
% %Rigid Body
% abaqus_spec.Part{2}.rigid_body.ATTR.analytical_surface='rigid_plate';
% abaqus_spec.Part{2}.rigid_body.ATTR.ref_node=1;

%%--> Assembly
abaqus_spec.Assembly.ATTR.name='Assembly-1';
abaqus_spec.Assembly.Instance{1}.ATTR.name='Cube-assembly';
abaqus_spec.Assembly.Instance{1}.ATTR.part='Cube';

% Reference Point
abaqus_spec.Assembly.Node{1}.VAL={1, [ GeomSize(3,1)/2 GeomSize(3,2)/2 GeomSize(1,3)*2]};

% abaqus_spec.Assembly.Instance{2}.ATTR.name='rigid_plate_1';
% abaqus_spec.Assembly.Instance{2}.ATTR.part='rigid_plate';
% abaqus_spec.Assembly.Instance{2}.VAL{1,1}={[ 0 0 BotPHeight]};
% abaqus_spec.Assembly.Instance{2}.VAL{2,1}={[ 0 0 BotPHeight 1 0 BotPHeight -90]};
% 
% abaqus_spec.Assembly.Instance{3}.ATTR.name='rigid_plate_2';
% abaqus_spec.Assembly.Instance{3}.ATTR.part='rigid_plate';
% abaqus_spec.Assembly.Instance{3}.VAL{1,1}={[ 0 0 TopPHeight]};
% abaqus_spec.Assembly.Instance{3}.VAL{2,1}={[ 0 0 TopPHeight 1 0 TopPHeight 90]};


abaqus_spec.Assembly.Nset{1}.ATTR.nset='Top';
abaqus_spec.Assembly.Nset{1}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Nset{1}.VAL=VTop;

abaqus_spec.Assembly.Nset{2}.ATTR.nset='Bottom';
abaqus_spec.Assembly.Nset{2}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Nset{2}.VAL=VBottom;

abaqus_spec.Assembly.Nset{3}.ATTR.nset='all';
abaqus_spec.Assembly.Nset{3}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Nset{3}.VAL=1:1:size(V_join,1);

abaqus_spec.Assembly.Nset{4}.ATTR.nset='XFace';
abaqus_spec.Assembly.Nset{4}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Nset{4}.VAL=bcXFace';

abaqus_spec.Assembly.Nset{5}.ATTR.nset='YFace';
abaqus_spec.Assembly.Nset{5}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Nset{5}.VAL=bcYFace';


abaqus_spec.Assembly.Nset{6}.ATTR.nset='Ref_point';
abaqus_spec.Assembly.Nset{6}.VAL=1;

% 
% abaqus_spec.Assembly.Nset{6}.ATTR.nset='TopFace';
% abaqus_spec.Assembly.Nset{6}.ATTR.instance='Cube-assembly';
% abaqus_spec.Assembly.Nset{6}.VAL=VTop;
% 
% abaqus_spec.Assembly.Nset{7}.ATTR.nset='BotFace';
% abaqus_spec.Assembly.Nset{7}.ATTR.instance='Cube-assembly';
% abaqus_spec.Assembly.Nset{7}.VAL=VBottom;

abaqus_spec.Assembly.Elset{1}.ATTR.elset='AllElA';
abaqus_spec.Assembly.Elset{1}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Elset{1}.VAL=element_IDs;

abaqus_spec.Assembly.Elset{2}.ATTR.elset='NodeElA';
abaqus_spec.Assembly.Elset{2}.ATTR.instance='Cube-assembly';
abaqus_spec.Assembly.Elset{2}.VAL=NodeEls;

% Define surfaace for coupling of top surface
abaqus_spec.Assembly.surface{1}.ATTR.type='Node';
abaqus_spec.Assembly.surface{1}.ATTR.name='Top_Surf';
abaqus_spec.Assembly.surface{1}.VAL='Top';
 
abaqus_spec.Assembly.Coupling.ATTR.constraint_name='Ref_coupling';
abaqus_spec.Assembly.Coupling.ATTR.ref_node='Ref_point';
abaqus_spec.Assembly.Coupling.ATTR.surface='Top_Surf';
abaqus_spec.Assembly.Coupling.Kinematic = ' ';

abaqus_spec.Surface_Interaction.VAL=1;
abaqus_spec.Surface_Interaction.Friction.VAL=Opts.FrictionCoeff;
abaqus_spec.Surface_Interaction.Time_Points.ATTR.NAME='Start';
abaqus_spec.Surface_Interaction.Time_Points.VAL=0;

% Displacement Amplitude
abaqus_spec.Amplitude.ATTR.name='Amp-1';
abaqus_spec.Amplitude.VAL=[ 0 0 Opts.Time 1];


%%--> Material
for i = 1:size(Opts.Material,2)
    abaqus_spec.Material{i}.ATTR.name=Opts.Material(i).Name;
    abaqus_spec.Material{i}.Elastic=Opts.Material(i).Elastic;
    abaqus_spec.Material{i}.Density=Opts.Material(i).Density;
    if Opts.Material(i).PlasticOn==1
        abaqus_spec.Material{i}.Plastic=Opts.Material(i).Plastic;
    end
end

%%--> Contact
%Defining contact properties
abaqus_spec.Surface_Interaction.ATTR.name='IntProp-1';
abaqus_spec.Surface_Interaction.VAL=1;
abaqus_spec.Surface_Interaction.Friction.VAL=Opts.FrictionCoeff;
abaqus_spec.Surface_Interaction.Time_Points.ATTR.NAME='Start';
abaqus_spec.Surface_Interaction.Time_Points.VAL=0;


% Contact between plates and struts, %% Contact exclusions Surface
% abaqus_spec.Assembly.surface{1}.ATTR.type='ELEMENT';
% abaqus_spec.Assembly.surface{1}.ATTR.name='NoCONTACTSurf';
% abaqus_spec.Assembly.surface{1}.VAL='ContactExclusions';

abaqus_spec.Surface_Interaction.Contact.ATTR.op='NEW';
abaqus_spec.Surface_Interaction.Contact_Inclusions.ATTR.ALL_EXTERIOR='';
abaqus_spec.Surface_Interaction.Contact_Property_Assignment.VAL={'' ,'','IntProp-1'};

% abaqus_spec.Surface_Interaction.Contact_Exclusions='';
% abaqus_spec.Surface_Interaction.Contact_Exclusions.VAL='NoCONTACTSurf';  

% 
% %%--> Boudnary Conditions for initial step
% Constrain bottom panel displacement 
abaqus_spec.Surface_Interaction.Boundary{1}.VAL{1,1}={'Bottom' , [1 3]};
% Contrain top  panel
abaqus_spec.Surface_Interaction.Boundary{1}.VAL{2,1}={'Top' , [1 2]};

% Constrain REF Node in all displacement except Z
abaqus_spec.Surface_Interaction.Boundary{1}.VAL{3,1}={'Ref_point' , [1 2]};

if Opts.BCsetting == 5  % Symmetry
    abaqus_spec.Surface_Interaction.Boundary{1}.VAL{4,1}={'XFace' , 1};
    abaqus_spec.Surface_Interaction.Boundary{1}.VAL{5,1}={'XFace' , [5 6]};
    abaqus_spec.Surface_Interaction.Boundary{1}.VAL{6,1}={'YFace' , 2};
    abaqus_spec.Surface_Interaction.Boundary{1}.VAL{7,1}={'YFace' , 4};
    abaqus_spec.Surface_Interaction.Boundary{1}.VAL{8,1}={'YFace' , 6};
end
% elseif Opts.BCsetting == 0   %Prevent rigid body movements
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{5,1}={'XFace' , 1};
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{6,1}={'YFace' , 2};
% elseif Opts.BCsetting == 2   %Constrained lateral displacement top and bottom
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{5,1}={'BotFace' ,[1 2]};
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{6,1}={'TopFace' ,[1 2]}; 
% elseif Opts.BCsetting == 3   %Constrained lateral displacement & rotation top and bottom
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{5,1}={'BotFace' ,[1 2]};
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{6,1}={'BotFace' ,[4 6]};
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{7,1}={'TopFace' ,[1 2]}; 
%     abaqus_spec.Surface_Interaction.Boundary{1}.VAL{8,1}={'TopFace' ,[4 6]}; 
% end



%%--> Step
abaqus_spec.Step.ATTR.name='Step-1';
abaqus_spec.Step.ATTR.nlgeom='YES';
abaqus_spec.Step.Dynamic.ATTR.Explicit='';
abaqus_spec.Step.Dynamic.VAL={' ' , Opts.Time};
if  Opts.MassScalingOn==1
    abaqus_spec.Step.Dynamic.Fixed_Mass_Scaling.ATTR.dt=Opts.MSIncSize;
    abaqus_spec.Step.Dynamic.Fixed_Mass_Scaling.ATTR.type='below min';
    elseif  Opts.MassScalingOn==2
    abaqus_spec.Step.Dynamic.Variable_Mass_Scaling.ATTR.dt=Opts.MSIncSize;
    abaqus_spec.Step.Dynamic.Variable_Mass_Scaling.ATTR.type='below min';
    abaqus_spec.Step.Dynamic.Variable_Mass_Scaling.ATTR.Number_interval=Opts.MSIntervals;
    
end


% Boundary Conditions for step 1
abaqus_spec.Step.Boundary.ATTR.type='Displacement';
abaqus_spec.Step.Boundary.ATTR.amplitude='Amp-1';
abaqus_spec.Step.Boundary.VAL={'Ref_point',[3 3],displacementMagnitude};


%Output
abaqus_spec.Step.Restart.ATTR.write='';
abaqus_spec.Step.Restart.ATTR.number_interval=20;
abaqus_spec.Step.Restart.ATTR.time_marks='NO';

abaqus_spec.Step.Output{1}.ATTR.field='';
abaqus_spec.Step.Output{1}.ATTR.number_interval=Opts.Num_intervals;
abaqus_spec.Step.Output{1}.ATTR.variable='PRESELECT';
abaqus_spec.Step.Output{2}.ATTR.history='';
abaqus_spec.Step.Output{2}.ATTR.variable='PRESELECT';
abaqus_spec.Step.Output{3}.ATTR.field='';
abaqus_spec.Step.Output{3}.ATTR.number_interval=Opts.Num_intervals;
abaqus_spec.Step.Output{3}.Element_Output.ATTR.Directions={'Yes'};
abaqus_spec.Step.Output{3}.Element_Output.VAL={'STATUS','EVOL'};

%%
if ~exist(savePath, 'dir')
   mkdir(savePath)
end

abaqusStruct2inpBMcD(abaqus_spec,abaqusInpFileName);

%textView(abaqusInpFileName);

%% Run the job using Abaqus
tic

lockFileName=fullfile(savePath,[Opts.FileName,'.lck']);
if exist(lockFileName,'file')
    warning('Lockfile found and deleted')
    delete(lockFileName);
end

%%

oldPath=pwd; %Get current working directory
cd(savePath); %Set new working directory to match save patch

runFlag=system([abaqusPath,' inp=',abaqusInpFileName,' job=',Opts.FileName,' interactive ask_delete=OFF',' double']);

cd(oldPath); %Restore working directory


AnalysisTime=toc;

%% Import and visualize abaqus results
% Run pythin script to create a report contaning reaction force and energy data
addpath('abq_out')

if exist('abq_out/RFReport.rpt', 'file')==2
    delete('abq_out/RFReport.rpt')
end
if exist('abq_out/Energy.rpt', 'file')==2
    delete('abq_out/Energy.rpt')
end

system(['abaqus cae noGUI=ExplicitReportSP.py -- ' savePath [' \' Opts.FileName '.odb']]);
delete('abaqus.rpy');

% Import data reaction force data
RFdata = importdata('RFReport.rpt',' ',5) ;
[timeVec,i_unique,~] = unique(RFdata.data(:,1));
RF = RFdata.data(i_unique,2); %E_Lengths

% Calculate displacement from time and displacment magnitude
raw_displacement = (abs(displacementMagnitude)*timeVec)/Opts.Time ;

% 
% % Adjust data to remove first n rows with 0 reaction force
% data_start_ind = find(RF ~= 0, 1, 'first') - 1;
% adj_displacement = raw_displacement(data_start_ind : end) - raw_displacement(data_start_ind);
% adj_RF = RF(data_start_ind : end);
% 
% % get strain from displacement magnitude and distance between plates
% % account for distance removed in previous step
% OvrStrain= adj_displacement/((TopPHeight - BotPHeight) - raw_displacement(data_start_ind)) ; 
% 
% % Get stress data from initial model dimensions and reaction force
% OvrStress = abs(adj_RF)/(GeomSize(3,1)*GeomSize(3,2));  

OvrStrain= raw_displacement/(TopPHeight - BotPHeight) ; 
OvrStress = abs(RF)/(GeomSize(3,1)*GeomSize(3,2));  

% Import energy data
Energydata = importdata('Energy.rpt') ;
Energies=Energydata.data;
Energies=[Energies (Energies(:,3)./Energies(:,2))];


%% Plot
if Plots ==1
    %% Plot Overall Stress/strain Curve
    figure
    plot(OvrStrain,OvrStress,'-*')
    ylabel('\sigma Compressive Stress (MPa)')
    xlabel('\epsilon Strain')
end

StressStrain=[OvrStrain OvrStress];

end

%%
% Lattice Inverse Design & Optimisation Tool 
% 06/12/2023 - Brian McDonnell - University of Galway
% GNU AFFERO GENERAL PUBLIC LICENSE - See LICENSE file details