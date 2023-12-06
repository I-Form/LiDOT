function [state, options, optchanged] = outfun_ga(options, state, flag)
% stop = outfun(x,optimValues,state)
global Gen GenData EvalData PopulationData Pop

if strcmp('iter',flag)
    Gen=Gen+1;
    msgfile = fopen('OptimisationStatus.txt','a+');
    fprintf(msgfile,'Start of Generation %d\n',Gen); 
    fclose(msgfile);

    GenData(Gen).state=flag;
    GenData(Gen).iteration=state.Generation;
    GenData(Gen).funccount=state.FunEval;
    GenData(Gen).x=state.Population;
    GenData(Gen).fval=state.Score;
    GenData(Gen).bestfval=state.Best;
    GenData(Gen).Data = PopulationData;
    PopulationData(:) = [];
    Pop = 0;
    % Save temp_data in case matlab crashes
    save('mat_out/temp_data.mat','EvalData','GenData')
end
optchanged = false; 

% below options force every point to be evaluates
% Important this is true to ensure all data is collected in EvalData
state.EvalElites = true;
state.HaveDuplicates =false;

end

%%
% Lattice Inverse Design & Optimisation Tool 
% 06/12/2023 - B.McDonnell - University of Galway
% GNU AFFERO GENERAL PUBLIC LICENSE - See LICENSE file details