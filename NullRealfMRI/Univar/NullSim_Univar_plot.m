clear


PWmethods = {'AR-YW_AR1','AR-YW_AR2','AR-YW_AR5','AR-YW_AR10',...
    'ACF_AR5','ACF_AR10','ACF_AR15','ACF_AR30','ACF_AR60',...
    'AR-W_AR1','AR-W_AR2','AR-W_AR5','AR-W_AR10','AR-W_AR15'};

BetaBias = @(Truth_MC,Estimated_TH) ((Estimated_TH-Truth_MC)./Truth_MC).*100;

Biases=[]; ZCPS=[];
for pwm = PWmethods
    disp(pwm{1})
    BB = load(['R/NullSim_Univar_Biases_' pwm{1} '.mat'],...
        'Beta_PW_SE_S','Beta_PW_SE_T','Bhat_Naive_S','Beta_Naive_SE_T',...
        'CPZ_Naive','CPZ_PW'); 
    Beta_PW_SE_S_STD = std(BB.Beta_PW_SE_S);
    Beta_PW_SE_T_MEAN = mean(BB.Beta_PW_SE_T);
    Biases = [Biases BetaBias(Beta_PW_SE_S_STD,Beta_PW_SE_T_MEAN)];
    
    ZCPS = [ZCPS mean(BB.CPZ_PW)]; 
    
end
    
Biases = [BetaBias(std(BB.Bhat_Naive_S),mean(BB.Beta_Naive_SE_T)) Biases];
ZCPS   = [mean(BB.CPZ_Naive) ZCPS];

PWmethods = ['Naive' PWmethods];

fh = figure; 
fh_bias = subplot(2,1,1);
hold on; grid on; box on; 
title('[THEORETICAL - TRUE (SIMULATIONS)] / SIMULATION (SIMULATIONS) x 100')
bh0 = bar(Biases);
fh_bias.XTick=1:numel(Biases);
fh_bias.XTickLabel = PWmethods;
fh_bias.XTickLabelRotation = 45;
fh_bias.FontSize = 12; 
set(gca,'TickLabelInterpreter','none')
ylabel('Bias in SE of $\hat\beta$','Interpreter','latex')

fh_CSP=subplot(2,1,2);
hold on; grid on; box on; 
title('Cumulative Periodgram Z-scores')
bh1 = bar(ZCPS);
fh_CSP.XTick=1:numel(Biases);
fh_CSP.XTickLabel = PWmethods;
fh_CSP.XTickLabelRotation = 45;
fh_CSP.FontSize = 12; 
set(gca,'TickLabelInterpreter','none')
ylabel('CPS (Z-scores)','Interpreter','latex')

set(fh,'color','w')






