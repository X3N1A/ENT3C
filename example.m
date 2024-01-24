clear all
%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
MAX_WINDOWS=1000;
CHRSPLIT=nan;
SUB_M_SIZE_FIX=300;
ChrNr=15;
Resolution=40e3;
WS=1;
NormM=0;
%%% example files are 
% ENCSR079VIJ_G401 BR1 and 2 DS to 60000000 interactions 
% 4DNESWST3UBH_HFFc6 BR1 and 2 DS to 800000000 interactions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_PATH='DATA_30e6';

CELL_TYPEs={'G401','HFFc6'};

NUMMEL_BRs=struct();N=0;
for CELL_TYPE=CELL_TYPEs
    CELL_TYPE=CELL_TYPE{1};
    eval(sprintf('NUMMEL_BRs.%s=1:%d;',CELL_TYPE,2))
end
    
INTERSECT_EMPTY_BINS=[];k=0;
for CELL_TYPE=CELL_TYPEs
    CELL_TYPE=CELL_TYPE{1};
    for BR=1:2
        IN=sprintf('%s/%s.BR%d.%dkb.cool',DATA_PATH,CELL_TYPE,BR,Resolution/1e3);

        [M,BIN_TABLE]=load_cooler(IN,ChrNr,Resolution,NormM);
        eval(sprintf('BIN_TABLE_%s%d=BIN_TABLE;',CELL_TYPE,BR))
        eval(sprintf('CMatrix_%s%d=M;',CELL_TYPE,BR))
        INTERSECT_EMPTY_BINS=[INTERSECT_EMPTY_BINS;BIN_TABLE.binNrCHRS(isnan(BIN_TABLE.weights)|isnan(BIN_TABLE.CONTACT))];
    end
end

VN_ENT=[];
for CELL_TYPE=CELL_TYPEs
    CELL_TYPE=CELL_TYPE{1};
    for BR=1:2
        eval(sprintf('BIN_TABLE=BIN_TABLE_%s%d;',CELL_TYPE,BR))
        eval(sprintf('CMatrix=CMatrix_%s%d;',CELL_TYPE,BR))

        BIN_TABLE_RED=BIN_TABLE;
        f = find(ismember(BIN_TABLE_RED.binNrCHRS,INTERSECT_EMPTY_BINS, 'rows'));
        BIN_TABLE_RED(f,:)=[];
        m=CMatrix(BIN_TABLE_RED.binNrCHRS,BIN_TABLE_RED.binNrCHRS);      
        %figure(1);clf;imagesc(log(m));
        tic
        [VN_ENT_TEMP] = ENT3C(m,CELL_TYPE,BR,ChrNr,BIN_TABLE_RED,Resolution,MAX_WINDOWS,CHRSPLIT,SUB_M_SIZE_FIX,WS);
        VN_ENT=[VN_ENT;VN_ENT_TEMP]
        toc
        
    end
end

[combinations] = get_pairwise_combs(CELL_TYPEs,NUMMEL_BRs);
figure(1);clf
for ct=1:size(combinations,1)

    COMBS=combinations(ct,:);
    CELL_TYPE1 = COMBS(1);BioRep1 = COMBS{2};
    CELL_TYPE2 = COMBS(3);BioRep2 = COMBS{4};

    S=[VN_ENT.S(strcmp(VN_ENT.Name,CELL_TYPE1)&VN_ENT.BR==BioRep1),...
        VN_ENT.S(strcmp(VN_ENT.Name,CELL_TYPE2)&VN_ENT.BR==BioRep2)];

    Q=corrcoef(S(:,1),S(:,2));
    Q=Q(1,2);

    figure(1);nexttile;plot(S(:,1));hold on;plot(S(:,2));
    title(sprintf('Chr%d %sBR%d vs %sBR%d: Q=%4.4f',...
        ChrNr,CELL_TYPE1{1},BioRep1,CELL_TYPE2{1},BioRep2,Q),'interpreter','latex')


    sprintf('Chr%d Resolution %d WN=%d windows of size %dx%d (WS=%d)\n\n%s BR%d\n%s BR%d\n\nQ=%4.4f',...
        ChrNr,Resolution,size(S,1),SUB_M_SIZE_FIX,WS,CELL_TYPE1{1},BioRep1,CELL_TYPE2{1},BioRep2,Q)
end




