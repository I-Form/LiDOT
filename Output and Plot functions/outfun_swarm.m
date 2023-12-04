function stop = outfun_swarm(optimValues,state)
% stop = outfun(x,optimValues,state)
global Gen GenData

Gen=Gen+1;
GenData(Gen).state=state;
GenData(Gen).Genation=optimValues.Genation;
GenData(Gen).funccount=optimValues.funccount;
GenData(Gen).x=optimValues.swarm;
GenData(Gen).fval=optimValues.swarmfvals;

GenData(Gen).bestx=optimValues.bestx;
GenData(Gen).bestfval=optimValues.bestfval;

stop=false;
end