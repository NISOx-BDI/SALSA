function V = ARMACovMat(par,T,p,q)
% a: the MA parameters
% b: the AR parameters

    b=par(1:p);
    a=[1 par(p+1:end)];
    
    B=diag(ones(p-1,1),-1);
    B(1,:)=b(1:p); 
    
    % Returns with LL=10^10 when parameters lead to non-stationarity
    eigB=max(abs(eig(B)));
    if eigB>=1 | max(abs(a))>10
        LL=10^10;
        disp('ARMACovMat:: parameters lead to non-stationary.');
        %return
    end
 
 
    m=max(p,q);
    if q>0
        
        % Compute the stationary distribution for dimension max(p,q)
        A=zeros(m,q);
        A(1,:)=a(2:end);
 
        B=diag(ones(m-1,1),-1);
        B(1,1:p)=b;
 
        CC=[B A;zeros(q,m) diag(ones(q-1,1),-1)];
 
        Ve=zeros(m+q,m+q);
        Ve(1,1)=1;
        Ve(1,m+1)=1;
        Ve(m+1,1)=1;
        Ve(m+1,m+1)=1;
 
        C=(eye((m+q)^2)-kron(CC,CC))\Ve(:);
        C=reshape(C,m+q,m+q);
        Vstationary=C(1:m,1:m);
    else
        B=diag(ones(p-1,1),-1);
        B(1,:)=b;
 
        Ve=zeros(p,p);
        Ve(1,1)=1;
 
        C=(eye(p^2)-kron(B,B))\Ve(:);
        Vstationary=reshape(C,p,p);
    end


    % MA component
    AA=eye(q+1);
    for i=1:q
        AA=AA+a(i+1)*diag(ones(q-i+1,1),-i);
    end
    AA=AA'*AA;
 
    aa=AA(1,:);
    AA=zeros(T,T);
    for i=1:q
        AA=AA+aa(i+1)*diag(ones(T-i,1),-i);
    end
    AA=AA+AA'+aa(1)*eye(T);
 
    % Infill the top max(p,q) submatrix of AA assuming the process is
    % initially at stationarity
    BB=eye(m);
    for i=1:p
        BB=BB-b(i)*diag(ones(m-i,1),-i);
    end
 
    A1=BB*Vstationary*BB';
    AA(1:m,1:m)=A1;
 
    % AR component
    W=eye(T);
    for i=1:p
        W=W-b(i)*diag(ones(T-i,1),-i);
    end
    invW=W\eye(T);
 
    V=invW*AA*invW';