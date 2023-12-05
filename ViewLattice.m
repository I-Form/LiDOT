function ViewLattice(INPUT,varargin)
%% Plots the lattice in 3D to view strut thickness, include ParmHistory to 
% LatticePlotFunc(Geometry or EvalData, Name value pairs)
% Geometry input as given by GenerateLattice
% EvalData as produced by Optimisation run

% --- Name Value Pairs ---
% "n_sides" - number of sides on cylinders for plotting - lower to speed up code
% "Generate"  - ==0 uses nodes/diameters as in Geometry. ==1 Calls GenerateLattice based on lattice params (can be changed etc...)
% "Config" - for EvalData - 1 shows full data, 2 is just stress strain curve and lattice
% "VarNames" - String array of var names for plots

% Plotting lattice geometry only - struct 1 or more geometries
% Tile on to tile each input in subplots (MAX 6)
% "TileOn" - 1 for subplots of multiple geometries
% "TileConfig" - layout of subplots [rows columns]
% "TileTitles" - titles of each subplot

%% Parse inputs - check if Geometry or EvalData is passed in
% sets animation on if more than one geometry inputs.
% plots stress strain data if EvalData is passed in
Fields=fieldnames(INPUT);
if strcmp(Fields{1},'variable_vals')  %If EvalData is passed in
    AnimateOn=1;                     % animate multiple curves/geoms
    EvalData=INPUT;  
    num_iters = size(EvalData,2);
    Geometry=vertcat(EvalData.Geometry);
    StressStrains={EvalData.StressStrain};
    variable_vals=vertcat(EvalData.variable_vals);
    
    GeomOnly=0;         %stress strain and geometry passed in
    if any(strcmp(Fields,'target_curve'))  %if req curve is on
        target_curve=EvalData(1).target_curve;  
       % Distpoints=vertcat(EvalData.Distpoints);
    end
    
else
    GeomOnly=1;         %only geometry passed in
    num_iters = 1;
    Geometry=INPUT;
    if numel(Geometry) > 1 
        AnimateOn=1;
    else
        AnimateOn=0;
    end
   
end

%Parse additional inputs 
p=inputParser;

addRequired(p,'Geometry')

addParameter(p,'Config',0)
addParameter(p,'n_sides',12)

% Options to plot multiple geometries
addParameter(p,'TileOn',0)
addParameter(p,'TileConfig',[1 3]);
addParameter(p,'TileTitles',[]);

addParameter(p,'Generate',0);
addParameter(p,'FitnessVals',[]);
addParameter(p,'VarNames',[]);   %Has to be string array

parse(p,INPUT,varargin{:})

TileOn=p.Results.TileOn;
Tile=p.Results.TileConfig;
TileTitles=p.Results.TileTitles;

Config = p.Results.Config;
VarNames=p.Results.VarNames;
n_sides = p.Results.n_sides;

%% Organise data for generating visualisation
iter_data = struct;    %Structure to contain lattice data for generating visualisation
if p.Results.Generate == 1  % Generate Geometry from lattice_params - can change size etc..
    for a = 1:num_iters
        %Generate new lattice
        geom_a = GenerateLattice(Geometry(a).lattice_params) ;
        
        %Get data from new geometry
        % diameter of each element
        iter_data(a).diameters = geom_a.Diameters;
        % vector descibing orientation
        iter_data(a).el_orient_vec = geom_a.V(geom_a.E(:,end-1),:) ...
                                    - geom_a.V(geom_a.E(:,1),:);  
        % length of each element
        iter_data(a).el_lengths = sqrt( ...
                                (iter_data(a).el_orient_vec(:,1)).^2 + ...
                                (iter_data(a).el_orient_vec(:,2)).^2 + ...
                                (iter_data(a).el_orient_vec(:,3)).^2  );  
                            
        % Vector describing position          
        iter_data(a).el_pos_vec = geom_a.V(geom_a.E(:,1),:); 
        
    end 
    if a==1
        Geometry = geom_a;
    end
else    %% use exact input lattice in Geometry
    for a = 1:num_iters  
        %Get data from new geometry
        % diameter of each element
        iter_data(a).diameters = Geometry(a).Diameters;

        % vector descibing orientation
        iter_data(a).el_orient_vec = Geometry(a).V(Geometry(a).E(:,end-1),:) ...
                                   - Geometry(a).V(Geometry(a).E(:,1),:);  
        % length of each element
        iter_data(a).el_lengths = sqrt( ...
                                (iter_data(a).el_orient_vec(:,1)).^2 + ...
                                (iter_data(a).el_orient_vec(:,2)).^2 + ...
                                (iter_data(a).el_orient_vec(:,3)).^2  );  
        % Vector describing position                    
        iter_data(a).el_pos_vec = Geometry(a).V(Geometry(a).E(:,1),:); 

    end 
end
%% Loft Settings
%cPar.numSteps=17; 
cPar.closeLoopOpt=1; 
cPar.patchType='quad';
cPar.numSteps=1;

ns=n_sides;
t=linspace(0,2*pi,ns);
t=t(1:end-1);

F_all=cell(num_iters,1);
V_all=cell(num_iters,1);
     
     
%% Create patch data of a cylinder for each element
for k=1:num_iters
    clear Vt; clear Ft; clear Vs; clear Fs; clear Vq; clear Fq
    % Loop through every element
    for i=1:size(iter_data(k).diameters,1)
        ElHeight=iter_data(k).el_lengths(i);  %%
        orient_vec=iter_data(k).el_orient_vec(i,:);  %%
        if size(iter_data(k).diameters,2) == 1

            r=iter_data(k).diameters(i)/2;  

            % Sketching profile 1
            x=r*cos(t); 
            y=r*sin(t); 
            z=zeros(size(x));
            V_bottomBR=[x(:) y(:) z(:)];
    
            % Sketching profile 2
            z2=ones(size(x))*ElHeight;

            V_topBR=[x(:) y(:) z2(:)];
        else
            r=iter_data(k).diameters(i,:)/2;  
            
            % Sketching profile 1
            x=r(1)*cos(t); 
            y=r(1)*sin(t); 
            z=zeros(size(x));
            V_bottomBR=[x(:) y(:) z(:)];
    
            % Sketching profile 2
            x2=r(2)*cos(t); 
            y2=r(2)*sin(t); 
            z2=ones(size(x))*ElHeight;
            V_topBR=[x2(:) y2(:) z2(:)];

        end

        % Combine array
        VcNew=[V_bottomBR;  V_topBR] ;

        %Find euler angles from direction vector of element
        %heading_angle=    atan2(sqrt(dir_vec(2)^2 + dir_vec(3)^2) ,dir_vec(1));
%         heading_angle=    asin(dir_vec(1)/sqrt(dir_vec(2)^2 + dir_vec(3)^2));
%         pitch_angle =  asin(dir_vec(3));
% 
%         EuAng = [pitch_angle heading_angle    0];

        %Find euler angles from direction vector of element
        u = [0 0 1];
        v=[ 0 orient_vec(2:3)];

        EuAng(1) = acos(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));
        if orient_vec(2)> 0
            EuAng(1) = -EuAng(1);
        end
        if EuAng(1) ~=0
            u=v; 
        end
        EuAng(2) = acos(max(min(dot(u,orient_vec)/(norm(u)*norm(orient_vec)),1),-1));
        if orient_vec(1) < 0
            EuAng(2)=-EuAng(2);
        end

        EuAng(3)=0;

        % Create rotation matrix from euler angles
        [R,~]=euler2DCM(EuAng); %The rotation matrix

        % Rotate the ends 
        VcNew=(R*VcNew')'; 

        % Apply Displacement to each end
        VcNew = VcNew+iter_data(k).el_pos_vec(i,:);

        % Take out top and bottom circle and get patch data
        V_bottoms=VcNew(1:size(VcNew,1)/2,:);
        V_tops=VcNew((size(VcNew,1)/2)+1:end,:);

        %Get patch data for element end
        F_tops=1:size(V_tops,1);

        %Create patch data for element cylinder surface
        [F,V]=polyLoftLinear(V_bottoms,V_tops,cPar);

        if i==1  %for first loop
            Fs = F; Vs =V;
            Vt=V_tops; Ft=F_tops;
        else
            %Join new patch data for cylinder surface to previous array
            [Fs,Vs]=joinElementSets({F,Fs},{V,Vs});

            %Join new patch data for ends to previous array
            [Ft,Vt]=joinElementSets({F_tops,Ft},{V_tops,Vt});
        end
        
    end
    
    %Convert element end to tri patch data
    [Ft,Vt]=patch2tri(Ft,Vt);
    
    if iscell(Ft)
        Ft=Ft{2};
    end

    %Convert to quads in order to join with cylinder surface
    [Fq,Vq]=tri2quad(Ft,Vt);
    
    %Join surface and end patch data to one array
    [F_all{k},V_k]=joinElementSets({Fs,Fq},{Vs,Vq});
    
    % Move bottom surface to zero
%     V_k(:,1)=V_k(:,1) - min(V_k(:,1));
%     V_k(:,2)=V_k(:,2) - min(V_k(:,2));
%     V_k(:,3)=V_k(:,3) - min(V_k(:,3));
     V_all{k}=V_k;
    
end

%% Plotting

if TileOn==1   %Create sub plots of every lattice
    AnimateOn = 0;
    numFigs=ceil(num_iters/prod(Tile));
    c=1;
    for i=1:numFigs
        cFigure; %Open figure 
        if num_iters-(c-1) >= prod(Tile)
            numPlots=prod(Tile);
        else
            numPlots=num_iters-(c-1);
        end

        for j=1:numPlots
            subplot(Tile(1),Tile(2),j)
            if isempty(TileTitles)
                title(['Iteration ' num2str(c)])
            else
                title(TileTitles{c})
            end
            gpatch(F_all{c},V_all{c},'g','none',1); %Add graphics object to animate
            axisGeom(gca,15);   
            axis(axisLim(vertcat(V_all{c}))); %Set axis limits statically    
            view(140,30);
            camlight headlight;  
            c=c+1;
        end
    end


else  %Create plot of inital geometry/first entry
    
    if GeomOnly==1 % Only 1 Geometry passed in
        figure
        lat=gpatch(F_all{1},V_all{1},[0.8 0.8 0.8],'none',1); 
        hold on
        % light blue [0.302 0.7451 0.933]
        % Orange - [0.9294 0.6941 0.1255]  Green - [0.4667 0.6745 0.1882]
        % Blue -    [0 0.4471 0.7412]   Red - [0.8510 0.3255 0.0980]
        % Greys - [0.8 0.8 0.8]

        if isempty(Geometry.lattice_params.MaterialTransition) == 0
            % Plot plane seperating the two material zones
            % Min and max dimension
            max_size =  max(V_all{1}(:)) - min(V_all{1}(:));

            plane_0 = [-0.5*max_size -0.5*max_size 0;
                     -0.5*max_size 0.5*max_size 0;
                     0.5*max_size 0.5*max_size 0;
                     0.5*max_size -0.5*max_size 0];

           for i = 1:size(Geometry.lattice_params.MaterialTransition,1)
                d = Geometry.lattice_params.MaterialTransition(i,4:6);
                n = Geometry.lattice_params.MaterialTransition(i,1:3);
                
                %Find euler angles from direction vector of element
                u = [0 0 1];
                v=[ 0 n(2:3)];
        
                EuAng(1) = acos(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1));
                if n(2)> 0
                    EuAng(1) = -EuAng(1);
                end
                if EuAng(1) ~=0
                    u=v; 
                end
                EuAng(2) = acos(max(min(dot(u,n)/(norm(u)*norm(n)),1),-1));
                if n(1) < 0
                    EuAng(2)=-EuAng(2);
                end
        
                EuAng(3)=0;
        
                % Create rotation matrix from euler angles
                [R,~]=euler2DCM(EuAng); %The rotation matrix
        
                % Rotate the ends 
                plane_new=(R*plane_0')'; 
        
                % Apply Displacement to each end
                plane_new = plane_new+d;
    
                gpatch([1 2 3 4],plane_new,i,'k',0.5)

          
%             x = [min(V_all{1}(:,1)); max(V_all{1}(:,1));  max(V_all{1}(:,1)); min(V_all{1}(:,1))];
%             y = [min(V_all{1}(:,2)); min(V_all{1}(:,2));  max(V_all{1}(:,2)); max(V_all{1}(:,2))];
%             z = [min(V_all{1}(:,3)); min(V_all{1}(:,3));  max(V_all{1}(:,3)); max(V_all{1}(:,3))];
%             
           
%             for i = 1:size(Geometry.lattice_params.MaterialTransition,1)
%                 d = Geometry.lattice_params.MaterialTransition(i,4:6);
%                 n = Geometry.lattice_params.MaterialTransition(i,1:3);
%     
%     
%                 for j = 1:4
%                     z_plane(j,i) = ( n(1)*(x(i) - d(1)) + n(2)*(y(i) - d(2)) - n(3)*d(3) ) / -1*n(3);
%                     
%                     y_plane(j,i) = ( n(1)*(x(i) - d(1)) + n(2)*(z(i) - d(3)) - n(2)*d(2) ) / -1*n(2);
%                     x_plane(j,i) = ( n(1)*(y(i) - d(2)) + n(2)*(z(i) - d(3)) - n(1)*d(1) ) / -1*n(1);
%                 end
%    
%                 gpatch([1 2 3 4],[x_plane(:,i) y_plane(:,i) z_plane(:,i)],i,'k',0.5)


            end

        end

        % Plot options
        axisGeom(gca,30);   
        view(120,15);
        set(gca,'linewidth',5)
        box off
        camlight headlight;

    elseif GeomOnly==0   %If there is stress strain data - eg. optimisation results
        Error = vertcat(EvalData(:).OUT);  %RMSE For each evaluation
        
        hf=cFigure; %Open figure 
        
        %% Lattice Visualisation
        if Config == 1
            subplot(2,3,[1,4]) 
        elseif Config ==2
            subplot(1,2,1) 
        end
        % Plots base mesh / unit cells of lattice
        lat=gpatch(F_all{1},V_all{1},[0.8 0.8 0.8],'none',1);    % light blue [0.302 0.7451 0.933]
        % Orange - [0.9294 0.6941 0.1255]  Green - [0.4667 0.6745 0.1882]
        % Blue -    [0 0.4471 0.7412]   Red - [0.8510 0.3255 0.0980]
        % Greys - [0.8 0.8 0.8]

        % Plot options
        axisGeom(gca,20); 
        %axis([0 5 0 5 0 5])
        %axis(repmat([min(V_all{1},[],'all') max(V_all{1}),[],'all'],3))
        axis(axisLim(vertcat(V_all{1}))); %Set axis limits statically    
        view(120,15);
        %title('Lattice Geometry')
        %axis off
        box off
        camlight headlight;
        
        %% Target and Current Stress/strain 
        %restructure stress strain data
        allStrains=cell(numel(StressStrains),1);
        allStresses=cell(numel(StressStrains),1);
        %normStresses=cell(numel(StressStrains),1);

        for i=1:numel(StressStrains)
            curve=StressStrains{i};
            allStrains{i}=curve(:,1);
            allStresses{i}=curve(:,2);
            %normStresses{i}=curve(:,2).*(1/max(curve(:,2)));
        end
        
        
        
        if Config == 1
            subplot(2,3,[2,5])  
        elseif Config ==2
            subplot(1,2,2) 
        end
        target_static1=plot(target_curve(:,1),target_curve(:,2),'r','LineWidth',2);
        hold on
        stress_curve=plot(allStrains{1},allStresses{1},'k','LineWidth',3);
        xlabel('Strain')
        ylabel('Stress (MPa)')
        xlim([0 max(target_curve(:,1))]); %max(target_curve(:,1))
        %ylim([0 max(vertcat(allStresses{:}))*1.1]);
        ylim([0 max(target_curve(:,2))*2]);
        title('Compressive Stress-Strain Curve')
        legend('Target Curve','Current','Location','south')

        hold off

        %% Variables
        if Config == 1
            subplot(2,3,3) 
            current_vars=plot(ones(1,size(variable_vals,2)),variable_vals(1,:),'*r');
            hold on
            var_legend=cell(size(variable_vals,2),1);
            for i=1:size(variable_vals,2)
                % Plot each variable per iteration
                plot(variable_vals(:,i));

                % Names - Need to store these
                if isempty(VarNames) %Variables names undefined
                    var_legend{i+1}=['Variable' num2str(i)];
                else
                    var_legend{i+1}=VarNames(i);
                end
            end
            var_legend{1} = 'Current Lattice Variables';
            legend(var_legend(:)','Location','eastoutside')
            hold off
            title('Variable Values')
            xlabel('Iteration')
            ylabel('Value')
        elseif Config ==2
             
        end
        % Current Variable


        %% Plot the vale of objective function for each iteration
        if Config == 1
            subplot(2,3,6)
             % Plot Graph
            plot(Error)
            hold on
            % Current Iteration
            current_value=plot(Error(1),'*r');
            hold off
            title('RMSE Error (MPa)')
            xlabel('Evaluation Number')
            ylabel('RMSE Error (MPa)')
        elseif Config ==2
             
        end


 
     
    end
end
    
% Plotting the simulated results using |anim8| to visualize and animate
if AnimateOn==1
    % Set up animation features
    animStruct.Time=1:num_iters; %The time vector    
     for qt=1:num_iters %Loop over time increments 
        
        if Config == 1
            %Set entries in animation structure
            animStruct.Handles{qt}=[lat lat ...
                                    stress_curve stress_curve  ...
                                    current_vars current_vars ...
                                    current_value current_value]; %Handles of objects to animate
            animStruct.Props{qt}={'Faces','Vertices',...
                                  'XData','YData',...
                                  'XData','YData',...
                                  'XData','YData'}; %Properties of objects to animate
            animStruct.Set{qt}={F_all{qt},V_all{qt},...              % Lattice 
                                allStrains{qt},allStresses{qt},...   % Stress/Strain
                repmat(qt,1,size(variable_vals,2)),variable_vals(qt,:),...                % Variables
                qt,Error(qt)};                         % Objective Function

        elseif Config ==2
               %Set entries in animation structure
            animStruct.Handles{qt}=[lat lat ...
                                    stress_curve stress_curve];
            animStruct.Props{qt}={'Faces','Vertices',...
                                  'XData','YData'}; %Properties of objects to animate
            animStruct.Set{qt}={F_all{qt},V_all{qt},...              % Lattice 
                                allStrains{qt},allStresses{qt}};                         % Objective Function

        end
        
     end
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
end


end