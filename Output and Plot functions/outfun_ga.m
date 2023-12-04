function [state, options, optchanged] = outfun_ga(options, state, flag)
% stop = outfun(x,optimValues,state)
global Gen GenData EvalData PopulationData Pop

Gen=Gen+1;
GenData(Gen).state=flag;
GenData(Gen).Genation=state.Generation;
GenData(Gen).funccount=state.FunEval;
GenData(Gen).x=state.Population;
GenData(Gen).fval=state.Score;
GenData(Gen).bestfval=state.Best;
GenData(Gen).Data = PopulationData;
state.EvalElites = true;
PopulationData(:) = [];
Pop = 0;
save('mat_out/temp_data.mat','EvalData','GenData')

optchanged = false; 
end

