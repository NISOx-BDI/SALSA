function MCp = GenMotionParam(MCp,MParamNum)
% MCp: T x 6 motion parameter
% MParamNum: scalar
%
% SA, Ox, 2020

if MParamNum     == 6
    MCp = MCp;
elseif MParamNum == 12
    MCp = [MCp,MCp.^2]; % 12 parameter motion 
elseif MParamNum == 24
    o6MCp   = MCp;            % 6 orig param
    so6MCp  = o6MCp.^2;       % 6 square orig
    do6MCp  = diff(o6MCp);    % 6 diff
    do6MCp  = [ do6MCp; zeros(1,size(do6MCp,2)) ];
    sd6MCp  = do6MCp.^2;      % 6 square of diff
    MCp     = [o6MCp,so6MCp,do6MCp,sd6MCp];
else
    MCp = [];
end
disp(['Number of motion parameter: ' num2str(MParamNum)])


    % z-score the motion parameters
    MCp = MCp./std(MCp); 
    MCp = MCp-mean(MCp); 

end