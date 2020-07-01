function [Sen,Spc,Acc,TP,FP,FN,TN] = SenSpec(e,g)

    e = e(:); 
    g = g(:);

    CG  = numel(g);  

    POS = sum(g); %number of all TRUE edges in GROUND TRUTH
    NEG = CG-POS; %number of all Negative edges in the Ground Truth

    TP  = sum( e .* g ); %ground is 1, effect is 1
    FN  = sum(~e .* g ); %ground is 1, noneffect 1

    FP  = sum( e .* ~g);
    TN  = NEG-FP;

    Sen = TP./(TP+FN);
    Spc = TN./(TN+FP);
    Acc = (TP+TN)./(TP+FP+FN+TN);
end