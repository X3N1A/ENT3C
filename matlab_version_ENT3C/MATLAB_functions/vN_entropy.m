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

R1=1:phi:N;R1=R1(1:WN);
R2=R1+SUB_M_SIZE-1;R2=R2(1:WN);
R=[R1',R2'];

M(M==0)=nan;
M=log(M);

fid = fopen('outputM.txt', 'w');
for rr=1:WN

    m=M(R(rr,1):R(rr,2),R(rr,1):R(rr,2));

    if any(sum(isnan(m)|m==0)>=(SUB_M_SIZE-1))% if any column is all nan/0, ignore this submatrix (0 comes from log(1)!)
        ENT=nan;
    else
        %E=[E;sum(isnan(m(:)))/numel(m)];
        m(isnan(m))=nanmin(m(:));

        P = corrcoef(m);%,'rows','complete');

        % to match julia: cor() --> diagonal elements in P are always one.
        % should be the case anyways? ... yes but this can still happen if
        % a column has constant values!!!
        P(logical(eye(size(P))))=1; 
        % for EV computation, elements where entire row and column == nan/zero --> P set to 0
        % P(all(m==0,1) & all(m==0,2))=0; --> cannot happen if diag is one...

        % set P to 0 for constant columns
        SDs = std(m)<eps; % m=[1,2,3;11,2,33;111,2,333]; ... std in 2nd column is zero. But not for float: m(2,:)=0.4
        P(SDs,:)=0;
        P(:,SDs)=0;
        P(isnan(P))=0;% can still happen if two columns identical...(also float int difference) m=[1,2,3;11,2,33;111,2,333]


        rho=P./size(P,1);

        lam = eig(full(rho));
        lam = lam(lam>eps);

        ENT = -sum(real(lam.*log(lam)));
        fprintf(fid,'Step %d: S = %4.7f\n',rr-1,ENT);
        clear rho
    end
    S=[S;ENT];
    

    BIN_TABLE_NEW = [BIN_TABLE_NEW;...
        [BIN_TABLE.binNr(R(rr,1)),BIN_TABLE.binNr(R(rr,2)),...
        BIN_TABLE.START(R(rr,1)),BIN_TABLE.END(R(rr,2))]];

end







