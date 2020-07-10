    
T    = 50; % length of time series
TR   = 1; % dummy TR
p    = 3; % order of model

%% a dummy design
EV1     = randn(T,1);
EV1     = EV1 - mean(EV1); 
X       = [ones(T,1),EV1]; 

%% Residual forming matrix
pinvX   = pinv(X); 
R       = eye(T)-X*pinvX;

%% HP Filter 
K = hp_fsl(T,100,TR); % get the filter 


%% FMRISTAT implementation
% http://www.math.mcgill.ca/keith/fmristat/toolbox/fmrilm.m
% search for 'invM='
M_biasred   = zeros(p+1);
for i=1:(p+1)
    Di                  = (diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1)); %symmetric selector
    for j=1:(p+1)
       Dj               = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1)); % symmetric selector
       M0(i,j)          = trace(R*Di*R*Dj)/(1+(i>1));
    end
end

%% Appendix A:
% tr( R D_{l} ) this will be illiminated using (1+(j==1))
% tr(RD_{l} R ( D_{j} + D_{j}')) 
for l=1:(p+1)
    Dl = diag(ones(1,T-l+1),l-1); % upper triangle
   for j=1:(p+1)
        Dj = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
        M1(l,j) = trace(R*Dl*R'*Dj);
   end
end

% The M0 and M1 is identical 

%% Filter K was added into the M

RK  = R*K;
for i=1:(p+1)
    Di = (diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1)); %symmetric selector
   for j=1:(p+1)
        Dj = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
        MK0(i,j) = trace(RK'*Di*RK*Dj)/(1+(i>1));
   end
end


RK  = R*K;
for l=1:(p+1)
    Dl = diag(ones(1,T-l+1),l-1); % upper triangle
   for j=1:(p+1)
        Dj = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
        MK1(l,j) = trace(RK'*Dl*RK*Dj);
   end
end