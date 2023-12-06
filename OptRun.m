function [OUT] = OptRun(varargin)
% Objective Function for optimisation
% Generates lattice based on input, runs FEA and returns error to target
% curve (as defined below)
global Eval EvalData PopulationData Pop 
%% Process Inputs 
FEA_flag = 0;  %
if nargin == 1
    IN = varargin{1};
elseif nargin == 2
    % If FEA flag == 1, FEA will run even if same variables are in EvalData
    IN = varargin{1};
    FEA_flag = varargin{2}; 
end

IN = [2 70 100 90 80 70];    % Variables Example
var_names = {'Num Cells'; 'Taper'; 'Diam 1'; 'Diam 2'; 'Diam 3'; 'Diam 4'};
%% Function Settings
% Set function mode. % 0 - Standard, run FEA no plots & return objective. 
func_mode=3;         % 1 - visualize lattice, 2 - run FEA and plot diff to target, 3 - check target curve, 4 - print ABQ input file

Objective=1;            % 0 - Curve fit error (vector), 1 - Curve fit error (scalar RMSE), 2 - 1/Strength , 3 - 1/Energy Absorption, 4 - MULTI-OBJ
strength_cutoff=0.1;    % Strength is calculated as max stress at strain <= strength_cutoff

%% Lattice Geometry Settings
if isstruct(IN) == 1
    lat_opts = IN;
else
    % Scale variables down
    IN(2:end)=IN(2:end)/100;  
    
    lat_opts.lattice_type="BCC"; % one of  "BCC" "BCCZ" "FCC" "FCCZ" "OCTET" "DUAL"
    lat_opts.Shape='SquareSym';   % 'Square' or 'Cylinder' . 'SquareSym' Halves inner cells for symmetry BCs
    lat_opts.ZDir='Axial';        %["Radial","Hoop","Axial"] orientation of vertical struts for cylinders
    
    % Nominal number of unit cells in [x,y,z] [radial,axial,inner-diameter] 
    lat_opts.Dim=[3 3 5]; 
    % Nominl unit cell size in mm  [x,y,z] or [radial,hoop,axial]
    cell_size = 4 ; 
    
    % Change size of unit cell depending on number of Z cells in input
    lat_opts.Size=[cell_size cell_size cell_size*(lat_opts.Dim(3)/IN(1))];
    % Change lat_opts.Dim to match input
    lat_opts.Dim(3)=IN(1);
    
    lat_opts.nominal_diam = IN(3); % Nominal diameter (used if no gradient is defined)
    
    % Set variables to empty arrays [] to turn off.
    % Strut taper - local variation in strut diameter. Applied to allstruts
    lat_opts.taper=[]; % [1 IN(2) 1];                          
    % Graded Strut thickness from bottom to top of lattice 
    lat_opts.gradient =[lat_opts.nominal_diam IN(4:end)];    
    lat_opts.AR_gradient=[];                                 % Graded Aspect Ratio
    lat_opts.layer_heights = [];                             % Graded Aspect Ratio via layer heights
    
    % Mesh Settings
    lat_opts.el_per_strut=6;             % Elements per strut (mesh density)
    lat_opts.element_type=1;             % 1 - B31, 2 - B32 
    lat_opts.node_el_multiplier=1.4;   % Increase in diameter at strut intersections

    % Material Transition Plane for bi-material lattices. 
    % Set func_mode=1 to view plane
    % 2nd Material Section in direction of Normal. Extend arrays to include more sections
    lat_opts.MaterialTransition = [0 0 1  7.5 7.5 7.5] ; %  [Normal   Point In plane]
    lat_opts.MaterialSectionIDs = [1 1];   

end


%% FEA Settings
Options.FileName='OptRun';   % ABAQUS file name   

% Boundary Conditons
Options.appliedStrain= -0.7;   %This is an approximate value based on height of lattice

Options.BCsetting=1;  % 0 - prevent rigid body movements, 1 - x and y symmetry,
% 2 - Displacement constraints top and bottom, 3 - Displacement & Rotation Top and bottom
% 4 - Sandwich panels with hex elements instead of rigid plates,  5 - Sandwich panels with symmetry BCs
% Else - No BC on lattice, Friction only

% For sandwich panels :
Options.PanelMesh_ZSize = [3 3];  % [Panel thickness in mm , No. of Elements through thickness]
% Number of Element sin x-y direction per unit cell - must be multiple of 2 if SquareSymm is used for lattice
Options.PanelMesh_XYSizeFactor = 2 ;  

% Contact friction coefficient
Options.FrictionCoeff=0.3;

% ABAQUS Standad or Explicit
FEMethod=2;     % 1 - Standard ,  2 - Explicit (better for large models with more contact)

% Material Options - 17-4PH Stainless Steel - ZB Data
Options.Material(1).Name='AM 17-4PH';
Options.Material(1).Density=7.8e-9;        % Units - tonne/mm^3    = g/cm^3 * 10-9  - Steel = 7.8e-9
Options.Material(1).Elastic=[163726 0.3];
Options.Material(1).PlasticOn=1;
Options.Material(1).Plastic=[443 0;            % Galway 17-4PH (tensile input)
                            537.81 0.00384
                            580.96 0.00752
                            674.89 0.0139
                            705.50 0.0185
                            792.82 0.0272
                            946.13 0.0488
                            1061.61 0.0703
                            1131.49 0.0932
                            1159.43 0.1108];  % UTS

% Options.Material(1).Plastic=[1300 0;      Literature
%                 1380 0.0085783
%                 1453.58 0.013648
%                 1498.05 0.0329];  % UTS

% Damage Settings
Options.Material(1).DamageOn=0;      % 0 - Damage/element deletion off, 1 - On
Options.Material(1).FractureStrain=0.12;   
Options.Material(1).DispAtFailure=0.0015;

% Multi Material Options Material 2 - Normal Direction of plance
Options.Material(2).Name='AM 316L';
Options.Material(2).Density=7.8e-9;    
Options.Material(2).Elastic=[222000 0.3];
Options.Material(2).PlasticOn=1;
Options.Material(2).Plastic=[300 0;          %  Literature 316L
                           327.2  0.005446
                           416.374 0.04287
                           520.714 0.10153
                           653.58 0.20163
                           786.44 0.3314
                           830.23 0.3893];  % UTS

Options.Material(2).DamageOn=0;
Options.Material(2).FractureStrain=0.40;
Options.Material(2).DispAtFailure=0.0015;

Options.MaxDeg = 0.95; %Max stiffness degredation for elements before they are deleted

%Implicit Settings
% [Initial Time increment, Step time period, minimum increment time, max increment time]
Options.StepParameters=[0.1 1 1e-25 0.1];   
Options.StrutContact=0;                     % Strut contact on (1) / off (0)

%Explicit settings
Options.Time=2;                 % Step time length  
Options.MassScalingOn=1;        % 0 off, 1 fixed, 2 variable
Options.MSIncSize=1e-05;        % Mass scaling - defined stable increment size  
Options.MSIntervals = 1000;     % Intervals for variable mass scaling (if enabled)
Options.Num_intervals=100;      % number of output interval

%% Define custom hex mesh to fill with lattice
% Example of using a custom base mesh for lattice.
% (lattice variables controlling cell size/aspect ratio are overwritten)
% optionStruct.sphereRadius=20;
% optionStruct.coreRadius=15;
% optionStruct.numElementsMantel=1;
% optionStruct.numElementsCore=4;
% optionStruct.makeHollow=1;
% optionStruct.outputStructType=2;
% 
% %Creating sphere
% [meshStruct]=hexMeshSphere(optionStruct);
% 
% % Un-comment to include:
% lat_opts.meshStruct = meshStruct;

%% Define Target Curve 
% Target Curve - Interpolated to get target curve
% Flat - Target A
%target_curve=[0 0 ; 0.025 7.5; 0.05 15; 0.25 15; 0.5  15; 0.75 15]; 

% Linear Target B - YS 10 MPa
%target_curve=[0 0; 0.025 5; 0.05 10; (0.05 + 0.7/3) 15; 0.75  25]; 

%  Target C densification
%target_curve=[0 0; 0.025 5; 0.05 10; 0.15 10; 0.3  10; 0.4 10; 0.5 60];

% Linear Target D - YS 15 MPa
% target_curve=[0 0; 0.025 7.5; 0.05 15;0.4 25; 0.75  35]; 

%Non-Linear increase - Target E
target_curve=[0 0; 0.025 5; 0.05 10; 0.3  12; 0.5 16; 0.6 19; 0.75 25];


target_yield = [0.05, 10]; % Used for graphing to check final curve

% Amplify error in specific region of stress-strain curve
error_region = [0.15 0.15];   % [centre of region (Strain value), size of region (+- both sides of centre)
error_amplification_factor = 2;   % Factor to multiply error by

%% Generate full curve and graph
% Smoothly interpolates between values to create full curve
num_points = 100;   %number of points to compare 
strain_target = linspace(target_curve(1,1), target_curve(end,1), num_points);  % Strain values of target curve

% First generates temporary intermediate curve:
% new points either side of defined points to help create smooth curve
temp_strain = 0;
for i=2:size(target_curve,1)
    if i==3    % Yield Point
        e_new=[target_curve(i,1)*0.9 target_curve(i,1)*1.5];
    else        
        e_new=[target_curve(i,1)*0.95 target_curve(i,1)*1.1];
    end
    temp_strain=[temp_strain e_new]; % Append new strain values
end
temp_strain = sort(temp_strain); % Sort points 

% Linearly interpolate from defined curve to get stress values
temp_stress=interp1(target_curve(:,1),target_curve(:,2),temp_strain,'linear','extrap');

% Interpolate with smoothing to get final target curve
target_stress=interp1(temp_strain,temp_stress,strain_target,'makima');

% 2% offset Line for visualisation
offset = [target_curve([1 2],1)+0.002 target_curve([1 2],2);
          (target_curve(2,1)*2)+0.002 target_curve(2,2)*2];

% Amplified error region
error_region_x = [error_region(1)-error_region(2) error_region(1)-error_region(2) ...
        error_region(1)+error_region(2) error_region(1)+error_region(2)];
error_region_y = [0 max(target_curve(:,2))*1.2 ...
        max(target_curve(:,2))*1.2 0];


if func_mode==3  % Plot Target curves
    figure
    subplot(2,1,2)
    plot(temp_strain,temp_stress,'-x')
    xlabel('Strain')
    ylabel('Stress (MPa)')
    ylim([0 max(target_curve(:,2))*1.2])
    title('Before Interpolation')

    subplot(2,1,1)
    amp_error = patch(error_region_x,error_region_y,'yellow','LineStyle','none');
    hold on
    plot(strain_target,target_stress,'-r','LineWidth',2)
    plot(offset(:,1),offset(:,2),'-k')
    plot([target_yield(1) target_yield(1)],[target_curve(1,2) target_curve(end,2)],':k')
    plot([target_curve(1,1) target_curve(end,1)],[target_yield(2) target_yield(2)],':k')
    xlabel('Strain')
    ylabel('Stress (MPa)')
    title('Final Target')
    %xlim([0 0.75])%
    ylim([0 max(target_curve(:,2))*1.2])
    legend({['Amplified Error Region x' num2str(error_amplification_factor)]...
        ,'Target','2% Offset','Specified Yield Point'},'Location','south')
    legend('boxoff')
    set(gca, 'Layer', 'top')
    return
end

%% Status message
% Status Message if optimisation is running
Eval = Eval+1;  
Pop = Pop +1 ;
msgfile = fopen('OptimisationStatus.txt','a+');
timenow = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(msgfile,'Evaluation %d\n',Eval); 
fprintf(msgfile,'%23s\n',timenow);         
fprintf(msgfile,'INPUTS: %s \n',num2str(IN,6));
disp(Eval)
%% Check if same variables have been input before
% if Eval is defined and optimisation is running, check if the same set of
% variables have been inputted before in this problem to reduce uneccesary
% runs of FEA
if exist('Eval', 'var')==1 & Eval > 1   % Optimisation is running
    if FEA_flag==0  % if FEA_flag = 1, FEA is forced regardless of input
        % List of all previous variables
        VarHistory=vertcat(EvalData(:).variable_vals);

        % Check if the same variabls have been inputted before
        VarRepeat = sum(IN==VarHistory,2)==size(VarHistory,2);
        if nnz(VarRepeat) > 0
            % Copy all data from previous entry
            index=find(VarRepeat);
            OUT = EvalData(index(1)).OUT ;
            EvalData(Eval).variable_vals=IN;
            EvalData(Eval).variable_names=var_names;
            EvalData(Eval).target_curve=EvalData(index(1)).target_curve;
            EvalData(Eval).StressStrain=EvalData(index(1)).StressStrain;
            %EvalData(Eval).normResult=EvalData(index(1)).normResult;
            EvalData(Eval).error_points=EvalData(index(1)).error_points;
            EvalData(Eval).Strength=EvalData(index(1)).Strength;
            EvalData(Eval).Energy=EvalData(index(1)).Energy;
            %EvalData(Eval).dist=EvalData(index(1)).dist; 
            EvalData(Eval).Time=EvalData(index(1)).Time;

            if FEMethod==2
                EvalData(Eval).Energies=EvalData(index(1)).Energies;
            end

            EvalData(Eval).Geometry=EvalData(index(1)).Geometry; 
            EvalData(Eval).OUT=EvalData(index(1)).OUT; 
            EvalData(Eval).attr_out=EvalData(index(1)).attr_out;
            %EvalData(Eval).print_out=EvalData(index(1)).print_out; 
            
            PopulationData(Pop) = EvalData(Eval);
            
            fprintf(msgfile,'%s : %9.5f \n', EvalData(Eval).attr_out,...
                                             EvalData(Eval).OUT);    
            fprintf(msgfile,'Final Strain: %9.3f \n\n', EvalData(Eval).StressStrain(end,1));    
            fclose(msgfile);
            return

        end
    end
end

%% Create Lattice

Geometry = GenerateLattice(lat_opts);  % Generate lattice geometry
%    save('EGeom','Geometry')

if func_mode==1       % Visualise in 3D
    disp(['Estimated volume fraction = ' num2str(Geometry.relative_density)]);
    ViewLattice(Geometry);

    if isempty(lat_opts.gradient)==0 
        % Plot the gradient in strut thickness
        heights = abs(min(Geometry.gradient_factors(:,2)))+Geometry.gradient_factors(:,2);
        figure
        plot(Geometry.gradient_factors(:,1),heights)
        hold on
        plot(Geometry.lattice_params.gradient,linspace(min(heights),...
            max(heights),numel(Geometry.lattice_params.gradient)),'xr')
        hold off
        xlim([0 1.2])
        ylabel('Lattice Height')
        xlabel('Strut Thickness')
    end
    return
end
%% RUN FEA
if Options.BCsetting == 4 || Options.BCsetting == 5 % For Hex sandwich panels
    [StressStrain,AnalysisTime,Energies]=ExplicitBeamCompressionSP(Geometry,Options);
elseif FEMethod ==1     % Abaqus Standard
    [StressStrain,AnalysisTime]=StandardBeamCompression(Geometry,Options);
    %[StressStrain,AnalysisTime]=ImplicitBeamFEA(Geometry,Options);
elseif FEMethod==2      % Abaqus Explicit
    [StressStrain,AnalysisTime,Energies]=ExplicitBeamCompression(Geometry,Options);
end


%% Process Results
max_strength = max(StressStrain(:, 2));
strength = max(StressStrain(StressStrain(:,1)<= strength_cutoff , 2));
Energy = trapz(StressStrain(:,1),StressStrain(:,2));

%% Get Error to Target
strain_vals = strain_target(strain_target <= StressStrain(end,1) );
target_stress_resize=target_stress(1:numel(strain_vals));

% Linearly interpolate FEA result to get stress at strain_vals
result_stress=interp1(StressStrain(:,1),StressStrain(:,2),strain_vals,'linear');

% Get difference between FEA stress and target stress
error = target_stress_resize - result_stress;

% Amplify Error around specified region (e.g yield point)
error_region = (strain_vals >= error_region(1) - error_region(2)) & (strain_vals <= error_region(1) + error_region(2));
error(error_region) = error(error_region) * error_amplification_factor;

XDistpoints=[strain_vals; strain_vals];
YDistpoints=[target_stress_resize;  result_stress];

if func_mode==2
    %Plot the stress-strain curve comparison
    figure
    plot(strain_vals,target_stress_resize,'k')
    hold on 
    plot(strain_vals,result_stress,'r','LineWidth', 2)
    % plot([V V]',[zeros(size(V,1),1) ones(size(V,1),1)]','r:','LineWidth',0.2)  %plot bars
    plot(XDistpoints,YDistpoints,'b:.','LineWidth',1)  %distances
    hold off
    legend(["Target" "Result"],'Location','south')
    xlim([-0.01 strain_vals(end)+0.01])
    xlabel('Strain')
    ylabel('Stress (MPa)')
   
    if FEMethod==2
        % If Explicit - plot Kinetic vs Internal energy
        figure
        yyaxis left
        plot(Energies(:,1),Energies(:,2),Energies(:,1),Energies(:,3))
        ylabel('Energy (10^-3 J)')
        hold on
        yyaxis right
        ylabel('Ratio Kinetic to Internal Energy')
        plot(Energies(:,1),Energies(:,4))
        xlabel('Time')
        legend({'Internal Energy', 'Kinetic Energy' ,'Ratio of Kinetic to Internal'})
        hold off
    end
    
end

switch Objective
    case 0    % array of each error reading (lsqnonlin)
        OUT = error;  % - Error - output to optimisaation algorithm
        attr_out = "Root Mean Square Error";
        print_out = sqrt(mean((error.^2)));
    case 1      % RMSE 
        OUT = sqrt(mean((error.^2)));
        attr_out = "Root Mean Square Error";
        print_out = sqrt(mean((error.^2)));
    case 2
        OUT = 1/strength;
        attr_out = "Strength";
        print_out = strength;
    case 3
        OUT = 1/Energy;
        attr_out = "Specific Energy Absorption";
        print_out = Energy;
    case 4      % RMSE and Volume
        OUT = [sqrt(mean((error.^2))), Geometry.total_element_volume];
        attr_out = "Root Mean Square Error, Total Element Volume";
        print_out = [sqrt(mean((error.^2))), Geometry.total_element_volume];
end




%% Update Global variables for optimisation
if exist('Eval', 'var')==1 & Eval > 0 %if Eval is defined
    EvalData(Eval).variable_vals=IN;
    EvalData(Eval).variable_names=var_names;
    EvalData(Eval).StressStrain=StressStrain;
    EvalData(Eval).target_curve=[strain_vals' target_stress_resize'];
    %EvalData(Eval).normalised_stress=normResult;
    EvalData(Eval).error_points={XDistpoints YDistpoints};
    EvalData(Eval).Strength=max_strength;
    EvalData(Eval).Energy=Energy;
    %EvalData(Eval).dist=error; 
    EvalData(Eval).Time=AnalysisTime;
    if FEMethod==2
        EvalData(Eval).Energies=Energies;
    end
    EvalData(Eval).Geometry=Geometry; 
    EvalData(Eval).OUT=OUT; 
    EvalData(Eval).attr_out=attr_out; 
    %EvalData(Eval).print_out=print_out;
    
     if abs(StressStrain(end,1)) < (abs(Options.appliedStrain)*0.9)   
        % Create Copy of obd if it doesn't run to completion (otherwise
        % it will be overwritten)
        % movefile(['abq_out\' Options.FileName '.odb'],['abq_out\' Options.FileName num2str(Eval) '.odb']); 
        % copyfile([Options.FileName '.odb'],[Options.FileName num2str(Eval) '.odb']);

        % Normalise error with respect to the max strain in the simulation
        % -  otherwise optimisation will favour simulations which fail
        OUT(1) =  OUT(1) * ( abs(Options.appliedStrain)/abs(StressStrain(end,1)) );
        EvalData(Eval).OUT=OUT; 
     end

    PopulationData(Pop) = EvalData(Eval);
    fprintf(msgfile,'%s : %9.5f \n',attr_out ,print_out );    
    fprintf(msgfile,'Final Strain: %9.3f \n\n',(StressStrain(end,1)));    
    fclose(msgfile);
end


end