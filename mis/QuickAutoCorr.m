% genYWARparam=[0.3412 0.0384 0.0249 0.1734 0.2236 0.0386 -0.0262 -0.0177 -0.0233 0.0787 0.1605 -0.1093 0.0340 -0.0050 0.0277 -0.0249 0.0237 0.0212 -0.0720 0.0342];
% %Generate a time series -------------------------
% % ac_tmp       = dYacorr(:,vi);    
% % R_tmp        = toeplitz(ac_tmp(1:SimMord));
% % r_tmp        = ac_tmp(2:SimMord+1);
% % genYWARparam = R_tmp\r_tmp;
% % ------------------------------------------------        
% ACMat                  = full(spm_Q(genYWARparam,T));
% [ACMatDecomp,spdflag]  = chol(ACMat);
% 
% if spdflag
%     disp([num2str(vi) ', used nearestSPD.'])
%     spdCOVmat    = nearestSPD(ACMat);
%     ACMatDecomp  = chol(spdCOVmat);
% end    
% 

function [sqrtmV,p1] = QuickAutoCorr(dRESac_adj,T)

[Ainvt,posdef]  = chol(toeplitz(dRESac_adj)); 
p1              = size(Ainvt,1); 

if p1<2; disp('QuickAutoCorr:: AR1!'); end;

A               = Ainvt'; 
sqrtmV          = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 

end

