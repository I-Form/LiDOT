function [stop, options, optchanged] = outfun_pat(optimValues,options,flag)
global Gen GenData

Gen=Gen+1;
GenData(Gen).x=optimValues.x;
GenData(Gen).state=flag;
GenData(Gen).Genation=optimValues.Genation;
GenData(Gen).funccount=optimValues.funccount;
GenData(Gen).fval=optimValues.fval;

stop=false;
optchanged = false;
end
