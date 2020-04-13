function [Trends,TrndOrd] = GenerateTemporalTrends(T,TR,ttmethod,numTremdOrCTFreq)
% Generates trends using splines. 
% bits and pieces from fmrilm.m of Keith Worsely. 
% 
% SA, Ox, 2020
    
    if nargin==2; ttmethod='spline'; end; 

    if ~exist('numTremdOrCTFreq','var') || isempty(numTremdOrCTFreq); numTremdOrCTFreq=[]; end
    switch lower(ttmethod)
        case{'spline'}
            [Trends,TrndOrd] = SplineTrendGen(T,TR,numTremdOrCTFreq);
        case{'dct'}    
            [Trends,TrndOrd] = DCTTrendGen(T,TR,numTremdOrCTFreq);
        case{'poly','polynomial'}
            [Trends,TrndOrd] = PolyTrendGen(T,TR,numTremdOrCTFreq);
        otherwise
            error('The methods of trends is not legit.')
    end
    disp(['> GenerateTemporalTrends:: ' ttmethod ', ' num2str(TrndOrd)])
end


function [Trends,numDCT] = DCTTrendGen(T,TR,CutOffFreq)
    %if nargin<3; CutOffFreq = 128; end
    if isempty(CutOffFreq); CutOffFreq = 128; end; 
    
    % built-in SPM
    %K = struct('RT',TR,'row',1:T,'HParam',SPMdefaultHPCutOff);
    %Ydct = spm_filter(K,Y'); %DCT is happenening withing the function 
    % Similar:        
    
    numDCT  = fix(2*(T*TR)/CutOffFreq+1);

    t       = (0:(T-1))';
    Trends  = zeros(size(t,1),numDCT);
    Trends(:,1)  = ones(size(t,1),1)/sqrt(T);
    for k=2:numDCT
        Trends(:,k) = sqrt(2/T)*cos(pi*(2*t+1)*(k-1)/(2*T));
    end
    
    % If I wanted to continue and detrend
    %dctY = Y' - DCT_hpf*(DCT_hpf'*Y');      
end

function [Trends,TrndOrd] = PolyTrendGen(T,TR,numPoly)

    if isempty(numPoly); numPoly = 1+floor(TR.*T/150); end;  
    Y         = (1:T); 
    trend     = ((2*Y-(max(Y)+min(Y)))./(max(Y)-min(Y)))';
    Trends    = (trend*ones(1,numPoly+1)).^(ones(T,1)*(0:numPoly));
    TrndOrd   = numPoly;
end


function [Trends,numSpline] = SplineTrendGen(T,TR,numTremd)

    %if nargin<3; numTremd = 3; end
    if isempty(numTremd); numTremd=3; end; 
    
    Y = 1:T;
    numSpline    = round(numTremd*TR*T/360);   
    trend        = ((2*Y-(max(Y)+min(Y)))./(max(Y)-min(Y)))';
    if numSpline<=3
      Trends    = (trend*ones(1,numSpline+1)).^(ones(T,1)*(0:numSpline));
    else
      Trends    = (trend*ones(1,4)).^(ones(T,1)*(0:3));
      knot      = (1:(numSpline-3))/(numSpline-2)*(max(Y)-min(Y))+min(Y);
      for k=1:length(knot)
         cut    = Y'-knot(k);
         Trends = [Trends (cut>0).*(cut./max(cut)).^3];
      end
    end
    % temporal_trend=temporal_trend*inv(chol(temporal_trend'*temporal_trend));
end