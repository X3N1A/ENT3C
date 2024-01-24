include("load_cooler.jl")
function [SIMILARITY,VN_ENT] = ENT3C(IN1,IN2,SHORT_NAMES,OUT,Resolution,ChrNrs,CHRSPLIT)

	" INPUT VARS
	
	       IN1 & IN2 ... filenames of two contact matrices to be compared
                 IN1="/HDD3/DocumentsHDD3/07_2022-/Mapping/ENCSR079VIJ_G401/PreProcessed/BioRep1/ENCSR079VIJ.BioRep1.30000000.40kb.cool"
                 IN2="/HDD3/DocumentsHDD3/07_2022-/Mapping/ENCSR079VIJ_G401/PreProcessed/BioRep2/ENCSR079VIJ.BioRep2.30000000.40kb.cool"
	       SHORT_NAMES=[fdg,df] ... names of contact matrices to save in output table
	       OUT='temp.txt' ... filename of output table
	       Resolution=40e3; ... resolution of cool/to extract from mcool file
	       ChrNrs=22 ... for which chromasome to
	       CHRSPLIT=10; ... determines window/sub-matrix size on which entropy values S(window) are calculated on
	       MAX_WINDOWS=500; ... maximum number of entropy values to compute (window shift is reduced until desired window number is reached)
	
	 OUTPUT VARS
	       SIMILARITY ... similarity score of two input matrices
	       VN_ENT ... output table with entropy values and other information
        "	
	
	VN_ENT=[];
	for ChrNr=ChrNrs
	
	INTERSECT_EMPTY_BINS=[];
	for k=1:2
        	if k==1
                        FN=IN1; 
                else
                        FN=IN2;
                end
	
        	C_MAP, BIN_TABLE=load_cooler(FN,ChrNr,Resolution);
	
        	eval(Meta.parse("BIN_TABLE_$k=BIN_TABLE;"))
                eval(Meta.parse("C_MAP_$k=C_MAP;"))
	
                INTERSECT_EMPTY_BINS=vcat(INTERSECT_EMPTY_BINS, BIN_TABLE.binNrCHRS[isnan.(BIN_TABLE.weights) .|| isnan.(BIN_TABLE.CONTACT)]); 
                print(size(INTERSECT_EMPTY_BINS))
	end
	
	for k=1:2
	
        eval(Meta.parse("BIN_TABLE=BIN_TABLE_$k;"))
        eval(Meta.parse("C_MAP=C_MAP_$k;"))
        BIN_TABLE_RED=BIN_TABLE[setdiff(1:size(BIN_TABLE, 1), INTERSECT_EMPTY_BINS), :];
        "f=intersect(BIN_TABLE.binNrCHRS, INTERSECT_EMPTY_BINS)"
        C_MAP=C_MAP[BIN_TABLE_RED.binNrCHRS,BIN_TABLE_RED.binNrCHRS];
	SUB_M_SIZE=round(size(C_MAP,1)/CHRSPLIT);
	
	WS=1;
	WN=1+floor((size(C_MAP,1)-SUB_M_SIZE)./WS);
	while WN>MAX_WINDOWS
        	WS=WS+1;
	        WN=1+floor((size(C_MAP,1)-SUB_M_SIZE)./WS);
	end
	
	WN=1+floor((size(C_MAP,1)-SUB_M_SIZE)./WS);
        WN=convert(Int, WN)
	R1=collect(1:WS:size(C_MAP,1));R1=R1[1:WN];
        R2=R1.+SUB_M_SIZE.-1;R2=R2[1:WN];R2 = map(Int, R2);
        R=vcat(R1,R2);
	
	for rr=1:WN
        	m=log(C_MAP(R(rr,1):R(rr,2),R(rr,1):R(rr,2)));
	        f=find(isnan(m));
        	m(f)=nanmin(m(:));
	        P = corrcoef(m,'rows','complete');
        	rho=P./size(P,1);
	
        	lam = eig(full(rho));
	        lam = lam(lam>0); % handle zero entries better: we want 0*log(0) = 0, not NaN
	
        	S = [log(size(rho,1)) -sum(real(lam.*log(lam)))];
	
        	VN_ENT=[VN_ENT;table(SHORT_NAMES(k),ChrNr,Resolution,SUB_M_SIZE,WN,WS,BIN_TABLE_RED.START(rr),BIN_TABLE_RED.END(rr),S(1),S(2),...
	        'VariableNames',{'Name','ChrNr','Resolution','Sub_M_Size','WinNrs','WS','start','end','Smax','S'})];
	
        	end
	end
	
	SIMILARITY=corrcoef(VN_ENT.S(strcmp(VN_ENT.Name,SHORT_NAMES(1))),VN_ENT.S(strcmp(VN_ENT.Name,SHORT_NAMES(2))));
	SIMILARITY=SIMILARITY(1,2);
	
	sprintf('Chr%d Resolution %d WN=%d windows of size %dx%d (WS=%d)\n\n%s\n%s\n\nSimScore=%4.4f',...
	ChrNr,Resolution,WN,SUB_M_SIZE,SUB_M_SIZE,WS,IN1,IN2,SIMILARITY)
	
	writetable(VN_ENT,OUT,'Delimiter','tab')
	
	return SIMILARITY, VN_ENT
end
		
