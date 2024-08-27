function [S, SUB_M_SIZE, WN, phi, BIN_TABLE_NEW] = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,PHI_MAX,phi,BIN_TABLE)

S=[];BIN_TABLE_NEW=[];
N=size(M,1);

if isempty(SUB_M_SIZE_FIX)
    SUB_M_SIZE=round(N/CHRSPLIT);
else
    SUB_M_SIZE = SUB_M_SIZE_FIX;
end

WN=1+floor((N-SUB_M_SIZE)./phi);
while WN>PHI_MAX
    phi=phi+1;
    WN=1+floor((N-SUB_M_SIZE)./phi);
end

R1=1:phi:size(M,1);R1=R1(1:WN);
R2=R1+SUB_M_SIZE-1;R2=R2(1:WN);
R=[R1',R2'];

M=log(M);
for rr=1:WN
    m=M(R(rr,1):R(rr,2),R(rr,1):R(rr,2));
    if all(isnan(m(:)))
        ENT=nan;
    else
        %E=[E;sum(isnan(m(:)))/numel(m)];
        m(isnan(m))=nanmin(m(:));

        P = corrcoef(m,'rows','complete');
        P(isnan(P))=0;
        % to match julia: cor() --> diagonal elements in P are always one. 
        % for EV computation, elements where entire row and column == nan/zero --> P set to 0
        P(logical(eye(size(P))))=1; 
        P(all(m==0,1) & all(m==0,2))=0;

        rho=P./size(P,1);

        lam = eig(full(rho));
        lam = lam(lam>0);

        ENT = -sum(real(lam.*log(lam)));
    end
    S=[S;ENT];

    BIN_TABLE_NEW = [BIN_TABLE_NEW;...
        [BIN_TABLE.binNr(R(rr,1)),BIN_TABLE.binNr(R(rr,2)),...
        BIN_TABLE.START(R(rr,1)),BIN_TABLE.END(R(rr,2))]];

end








