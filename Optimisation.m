% This script 

clear all
% Change to correct directory
filepath = fileparts(mfilename('fullpath'));
cd filepath
% Declaration of global variables to store data
global Gen GenData Eval EvalData Pop PopulationData

% File Names and option to load in initial polulation to restart
load_initial_data = 0;          % Inital population loaded if == 1.
inital_data_file = 'TEST111.mat';
OptName='Optimisation_test1';     % Name of the optimum .obd file and the final data

% Optimisation Options for GA
Algorithm=4;            % 0-Lsqnonlin, 1-Fmincon, 2-Patternsearch, 3-Swarm, 4-GA 5 - GA MO
seed = 1;               % RNG seed
max_gen=10;             % Max number of generations for GA, Particle swarm an multi objective
pop_size=5;            % Population size

% Make folder to contain all data on optimisation case and save a copy of
% optimisation and optrun used.
mkdir(['mat_out/' OptName]);
copyfile('Optimisation.m',['mat_out/' OptName])
copyfile('OptRun.m',['mat_out/' OptName])

% Load in past results to re-start from previous optimisation case
if load_initial_data == 1 
    load(inital_data_file)    
    Eval = size(EvalData,2) ;
    Gen = size(GenData,2) ;
    pop_initial = GenData(end).x ;

    %Initialise Data Structures
    PopulationData = struct('variable_vals',[], 'variable_names',[],'StressStrain'...
        ,[],'target_curve',[],'error_points',[],'Strength',[],'Energy',[],'Time',...
        [],'Energies',[],'Geometry',[],'OUT',[],'attr_out',[]);
else
    %Initialize counters
    Gen=0; 
    Eval=0;
    Pop = 0;

    %Initialise Data Structures
    PopulationData = struct('variable_vals',[], 'variable_names',[],'StressStrain'...
         ,[],'target_curve',[],'error_points',[],'Strength',[],'Energy',[],'Time',...
         [],'Energies',[],'Geometry',[],'OUT',[],'attr_out',[]);
    EvalData = struct('variable_vals',[], 'variable_names',[],'StressStrain'...
        ,[],'target_curve',[],'error_points',[],'Strength',[],'Energy',[],'Time',...
        [],'Energies',[],'Geometry',[],'OUT',[],'attr_out',[]);
end


% Define problem structure
% Objective function
Objective=@OptRun;

% Initial parameters
% Number of Z cells - Taper - Nominal/Bottom Diam - Gradient Factors
x0 = [3 100 100 100 100 100];

nvars = numel(x0);

% Define upper and lower bounds - scaled to integers
lb = [ 3  70  50  50  50  50];
ub = [6 100 150 150 150 150];

% Inequality Constraints (Constraint applied gradient increases towards top
% of lattice)    A*x <= b

A = [0 0 0 0 -1 1;
     0 0 0 -1 1 0;
     0 0 -1 1 0 0];

b = [0;
     0;
     0];

% Equality Constraints   Aeq*x = beq
Aeq = [];
beq = [];

% Status file
msgfile = fopen('OptimisationStatus.txt','a+');
start = '%%%%%%%% Optimisation Beginning %%%%%%%%';
fprintf(msgfile,'\n%s\n\n',start);

 switch Algorithm
    case 0  %Lsqnonlin
        problem.objective = Objective;
        problem.x0=x0;
        problem.lb = lb;
        problem.ub = ub; 
        problem.solver='lsqnonlin';
        problem.options=optimoptions(@lsqnonlin,'TolFun',1e-6,'TolX',1e-6,...
                                    'StepTolerance',1e-6,'MaxGen',1000,...
                                    'DiffMinChange',5,...
                                    'Algorithm','trust-region-reflective',...
                                    'Display','Gen','Outputfcn',@outfun);
        [x, Resnorm,~,~,output]=lsqnonlin(problem);    

    case 1  %Fmincon
        problem.objective = Objective;
        problem.x0=x0;
        problem.Aineq = A;
        problem.bineq = b;
        problem.Aeq = Aeq;
        problem.beq = beq;
        problem.lb = lb;
        problem.ub = ub; 
        problem.nonlcon=[];
        problem.solver='fmincon'; 
        problem.options=optimoptions(@fmincon,...
                                    'Algorithm','sqp',...
                                    'DiffMinChange',5,...
                                    'Display','Gen','Outputfcn',@outfun);
        [x, fval,~,output]=fmincon(problem);    

    case 2  %Patternsearch
        problem.objective = Objective;
        problem.x0=x0;
        problem.Aineq = A;
        problem.bineq = b;
        problem.Aeq = Aeq;
        problem.beq = beq;
        problem.lb = lb;
        problem.ub = ub; 
        problem.nonlcon=[];           
        problem.solver='patternsearch';
        problem.rngstate=rng(seed);
        problem.options=optimoptions(@patternsearch,...
                                    'Display','Gen','Outputfcn',@outfun_pat...
                                    ,'MaxFunEvals',200,...
                                    'Cache','on','CacheTol',0.5);
        [x, fval,~,output]=patternsearch(problem);  

    case 3  % Particle Swarm
        problem.solver='particleswarm';
        problem.objective = Objective;
        problem.nvars=nvars;
        problem.lb = lb;
        problem.ub = ub;   
        problem.rngstate=[];
        problem.options=optimoptions(@particleswarm,...
                                    'Display','Gen','Outputfcn',@outfun_swarm,...
                                    'SwarmSize',10,'MaxGenations',50);
        [x, fval,~,output]=particleswarm(problem);   

    case 4  % Genetic Algorithm
        problem.fitnessfcn= Objective;
        problem.nvars=nvars;
        problem.Aineq = A;
        problem.bineq = b;
        problem.Aeq = Aeq;
        problem.beq = beq;
        problem.lb = lb;
        problem.ub = ub;    
        problem.nonlcon=[]; 
        problem.rngstate=rng(seed);
        problem.solver='ga';
        problem.intcon=1:nvars;
        problem.options=optimoptions(@ga,...
                                    'Display','Gen','Outputfcn',@outfun_ga,...
                                    'PopulationSize',pop_size,...                       
                                    'MaxGenerations',max_gen,...%   'TimeLimit', 60*60*20, ...
                                    'MaxStallGenerations',50 ,...  % Default = 50
                                    'FunctionTolerance',1e-6, ...  % Default = 1 e-6
                                    'PlotFcn', {@gaplotrange, @plotfun_ga});
        if load_initial_data == 1 
            problem.options = optimoptions(problem.options,...
                'InitialPopulationMatrix',pop_initial);
        end
        [x, fval,~,output]=ga(problem);   

     case 5  % Multi Objective Optimisation
        problem.fitnessfcn = Objective;
        problem.nvars=nvars;
        problem.Aineq = A;
        problem.bineq = b;
        problem.Aeq = Aeq;
        problem.beq = beq;
        problem.lb = lb;
        problem.ub = ub;    
        problem.nonlcon=[]; 
        problem.rngstate=rng(seed);
        problem.solver='gamultiobj';
        problem.intcon=1:nvars;
        problem.options=optimoptions(@gamultiobj,...
                                    'Display','Gen','Outputfcn',@outfun_gamultiobj,...
                                    'PopulationSize',pop_size,...                                  
                                    'MaxGenerations',max_gen,...%   'TimeLimit', 60*60*20, ...
                                    'MaxStallGenerations',50 ,...  % Default = 50
                                    'FunctionTolerance',1e-6, ...  % Default = 1 e-6
                                    'PlotFcn', @gaplotpareto);
        if load_initial_data == 1 
            problem.options = optimoptions(problem.options,...
                'InitialPopulationMatrix',pop_initial);
        end
        [x, fval,~,output]=gamultiobj(problem);   

 end

%Run with optimum parameters
[~]=OptRun(x,1);
movefile('abq_out/OptRun.odb',['abq_out/' OptName '.odb']);  %rename optimimum solution odb file

OptParam=x/100;    %Scale back optimim parameters    
 
finish = '%%%%%%%% Optimisation Complete %%%%%%%%%';
timenow = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(msgfile,'%23s\n',timenow);  
fprintf(msgfile,'\n%s\n\n',finish); 
fprintf(msgfile,'OPTIMUM: %s',num2str(OptParam,5));
fclose(msgfile);

% Save results to .mat file
save(['mat_out/' OptName '.mat'],'EvalData','GenData','x','fval')
 
%% Visualise Results - 1) Prepare common data
% Data common to all results Visualisations

% Variable names
Vars = ["Number of Cells"; "Taper"; "Diam 1"; "Diam 2"; "Diam 3";"Diam 4"];  
% Set Eval = number of objective function evaluations
Eval=size(EvalData,2);   

% According to Function Eval
fvals = vertcat(EvalData(:).OUT); % All raw fitness values
[best,min_i]=min(fvals(:,1),[],1);         % Minimum and it's index
[fv_sort, is]=sort(fvals(:,1));            % Sorted best to worst, index is in fvals/EvalData
[fv_uni, iu] = unique(fv_sort);       % Unique fitness values, sorted, index iu in fvals/EvalData

% Per generation
gen_fvals=horzcat(GenData(1:end-1).fval);   % Includes Elites
[mins, min_IDs] = min(gen_fvals);

Gen_ids = vertcat(GenData(1:end-1).Genation);  
best_fvals = [min(GenData(1).fval) GenData(end).bestfval];   

for i = 1:numel(Gen_ids)
    gen_best_data(i) = GenData(i).Data(min_IDs(i));
end
 %% Create graphs showing optimisation process - Per Function Evaluation - Normalised Result
%Optimal Stress-Strain Curve
figure
hold on
Colours=repmat(linspace(0.95,0.8,Eval-1),[3 1]);
for i=1:Eval-2
    plot(EvalData(i).StressStrain(:,1),...
         EvalData(i).StressStrain(:,2),...
         'LineWidth',0.01,'Color',Colours(:,i))  
end
Gen = plot(EvalData(i).StressStrain(:,1),...
            EvalData(i).StressStrain(:,2),...
            'LineWidth',0.01,'Color',Colours(:,i)) ;
opt=plot(EvalData(min_i).StressStrain(:,1),...
         EvalData(min_i).StressStrain(:,2),...
         '-k','LineWidth',4);
tar=plot(EvalData(1).target_curve(:,1),...
         EvalData(1).target_curve(:,2),...
         '-r','LineWidth',2);
hold off
xlabel('Strain')
ylabel('Stress (MPa)')
%title('Progress of Stress Strain Response')
legend([opt tar ],{'Optimal Solution','Target'},'Location','south')
legend('boxoff')
ylim([0 max(EvalData(1).target_curve(:,2))*2])
xlim([0 0.75])


%% Variables vs function evaluations
% Seperate into 3 graphs - cells, taper, diameters
figure                  
plot(1:Eval-1,vertcat(EvalData(1:end-1).variable_vals))   
xlabel('Function Evaluations')
ylabel('Parameters')
hold on                 %highlights the optimum parameters
plot(Eval,EvalData(Eval).variable_vals,'rx')
for i=1:numel(Vars)
    text(Eval,EvalData(Eval).variable_vals(i),...
        ['  ' num2str(EvalData(Eval).variable_vals(i))])
end
hold off
%legend([LayersLegend(:)'],'Location','eastoutside')
legend(Vars,'Location','eastoutside')

return

%% Objective Function Values vs Generation
figure
plot(best_fvals)
hold on
% plot(gen_fvals','.')
hold off
xlabel('Generation')
ylabel('Fitness Values (MPa)')
legend({'Best Point per generation'})


%% Stress Strain Curves for Best design in each generation
gen_colours=repmat(linspace(0.8,0.2,numel(Gen_ids)),[3 1]);
figure
target = plot(GenData(1).Data(1).target_curve(:,1),...
              GenData(1).Data(1).target_curve(:,2),...
              '-r','LineWidth', 3);
hold on
% Method using MATLAB data
% for i= 1:numel(Gen_ids)
%     plot(GenData(i).Data(min_IDs(i)).StressStrain(:,1),...
%          GenData(i).Data(min_IDs(i)).StressStrain(:,2),...
%          '-k','LineWidth', 1,'Color',gen_colours(:,i));
% end

% Method using data in GenData
for i = 1:numel(Gen_ids)
    plot(GenData(i).Data(min_IDs(i)).StressStrain(:,1),...
         GenData(i).Data(min_IDs(i)).StressStrain(:,2),...
        'LineWidth',1,'Color',gen_colours(:,i)); 
end
hold off
xlabel('Strain')
ylabel('Stress (MPa)')


%% Plot Best n solutions
n = 10;
figure
hold on
Colours=repmat(linspace(0.7,0.4,10),[3 1]);
for i = 1:n
    plot(EvalData(is(iu(i))).StressStrain(:,1),...
        EvalData(is(iu(i))).StressStrain(:,2),...
        'LineWidth',0.01,'Color',Colours(:,i))  
end
opt=plot(EvalData(is(iu(1))).StressStrain(:,1),...
         EvalData(is(iu(1))).StressStrain(:,2),...
         '-k','LineWidth',4);
tar=plot(EvalData(1).target_curve(:,1),...
         EvalData(1).target_curve(:,2),...
         '--r','LineWidth',2);
hold off

xlabel('Strain')
ylabel('Stress (MPa)')
%title('Progress of Stress Strain Response')
legend([opt tar],{'Optimal Solution','Target'},'Location','south')
legend('boxoff')
title("Best %s Designs")
%ylim([0 max(EvalData(1).target_curve(:,2))*1.5])
%xlim([0 0.5])


%% Plot residuals to target
i = min_i;

figure
plot(EvalData(1).target_curve(:,1),EvalData(1).target_curve(:,2),'r','LineWidth', 2)
hold on 
plot(EvalData(i).StressStrain(:,1),EvalData(i).StressStrain(:,2),'k','LineWidth', 2)
% plot([V V]',[zeros(size(V,1),1) ones(size(V,1),1)]','r:','LineWidth',0.2)  %plot bars
plot(EvalData(i).error_points{1},EvalData(i).error_points{2},'-b.','LineWidth',0.2,'MarkerSize',1)  %distances
hold off
legend(["Target" "Result" "Residuals"],'Location','south')
% xlim([-0.01 strain_vals(end)+0.01])
%xlim([0 0.75])
xlabel('Strain')
ylabel('Stress (MPa)')


%% Animate Geometry & Stress Strain Results - per function evaluation
% Figure window with animation of the optimisation progression
% Can take a few minutes to finish with many Genations
% Config = 1 or 2 . 2 -> Curve and Geom only
ViewLattice(EvalData,"Config",1,"VarNames",Vars) ;

%% Animate Geometry & Stress Strain Results - Best for each Generation
ViewLattice(gen_best_data,"Config",1,"VarNames",Vars);


