function [state, options, optchanged] = plotfun_ga(options, state, flag)
% stop = outfun(x,optimValues,state)
global Gen GenData EvalData PopulationData Pop
n = 10;
%% Plot Best n solutions
% According to Function Eval
fvals = vertcat(EvalData(:).OUT); % All raw fitness values
[best,min_i]=min(fvals(:,1),[],1);         % Minimum and it's index
[fv_sort, is]=sort(fvals(:,1));            % Sorted best to worst, index is in fvals/EvalData
[fv_uni, iu] = unique(fv_sort);       % Unique fitness values, sorted, index iu in fvals/EvalData



Colours=repmat(linspace(0.7,0.4,n),[3 1]);
for i = 1:n
    plot(EvalData(is(iu(i))).StressStrain(:,1),...
        EvalData(is(iu(i))).StressStrain(:,2),...
        'LineWidth',0.01,'Color',Colours(:,i))  
end
hold on
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
title("Best 10 Designs")
subtitle(["RMSE of Best Design = " num2str(min(fv_uni))])

optchanged = false; 
end

