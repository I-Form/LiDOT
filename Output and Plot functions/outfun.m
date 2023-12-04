function stop = outfun(x,optimValues,state)
global Gen GenData EvalData

Gen=Gen+1;
GenData(Gen).x=x;
GenData(Gen).state=state;
GenData(Gen).Genation=optimValues.Genation;
GenData(Gen).funccount=optimValues.funccount;

if exist('optimValues.gradient','var')
    GenData(Gen).gradient=optimValues.gradient;
end
if exist('optimValues.fval','var')
    GenData(Gen).fval=optimValues.fval;
end

if exist('optimValues.trustregionradius','var')
    GenData(Gen).radius=optimValues.trustregionradius;
end    
if exist('optimValues.residual','var')
    GenData(Gen).residual=optimValues.residual;
    GenData(Gen).fval=optimValues.resnorm;
end

save('mat_out/temp_data.mat','EvalData','GenData')

stop=false;
end