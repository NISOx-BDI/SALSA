function where2stop = FindBreakPoint(acf,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018
    if ~sum(ismember(size(acf),T)); error('There is something wrong, mate!'); end
    if size(acf,2) ~= T; acf = acf'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T);
    %idx        = find(abs(acs)>bnd);
    isit       = abs(acf)>bnd;
    %where2stop = find(isit==0); %finds the break point -- intercept 

    where2stop = zeros(1,size(acf,1));
    for i = 1:size(acf,1)
        where2stop(i) = find(~isit(i,:),1)-1;
    end
%     if where2stop(1)==1
%         where2stop = 0; 
%     else
%         where2stop = where2stop(1)-1; 
%     end;
end