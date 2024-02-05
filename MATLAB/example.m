clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WN_MAX=1000;
CHRSPLIT=7;
SUB_M_SIZE_FIX=nan;
ChrNr=14;
Resolution=40e3;
WS=1;
NormM=0;
DATA_PATH='DATA_30e6';
FILES = {sprintf("%s/ENCSR079VIJ.BioRep1.mcool",DATA_PATH),"G401_BR1",
         sprintf("%s/ENCSR079VIJ.BioRep2.mcool",DATA_PATH),"G401_BR2",
         sprintf("%s/ENCSR444WCZ.BioRep1.mcool",DATA_PATH),"A549_BR1",
         sprintf("%s/ENCSR444WCZ.BioRep2.mcool",DATA_PATH),"A549_BR2"};
OUT_DIR='OUTPUT/MATLAB';
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INFO = struct('FN', '', 'META', '');
FNs = [];
for f = 1:size(FILES,1)
    FNs = [FNs, INFO];
    FNs(end).FN = FILES{f, 1};
    FNs(end).META = FILES{f, 2};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract common empty bins of input matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXCLUDE=[];
for f = 1:numel(FNs)
        FN=FNs(f).FN;
        [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM);
        EXCLUDE=[EXCLUDE;BIN_TABLE.binNr(isnan(BIN_TABLE.weights)|isnan(BIN_TABLE.CONTACT))];
end
EXCLUDE=unique(EXCLUDE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% produce ENT3C data frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENT3C_OUT=[];
for f = 1:numel(FNs)
    FN=FNs(f).FN;

    [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM);
    INCLUDE = 1:size(M,1);
    INCLUDE = setdiff(INCLUDE,EXCLUDE);

    M = M(INCLUDE,INCLUDE);
    BIN_TABLE = BIN_TABLE(INCLUDE,:);

    [S, SUB_M_SIZE1, WN1, WS1, BIN_TABLE_NEW] = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,BIN_TABLE);

    N = length(S);

    OUT1 = table(repmat(FNs(f).META,N,1),...
        repmat(ChrNr,N,1),repmat(Resolution,N,1),...
        repmat(WN1,N,1),repmat(WS1,N,1),...
        BIN_TABLE_NEW(:,1),BIN_TABLE_NEW(:,2),...
        BIN_TABLE_NEW(:,3),BIN_TABLE_NEW(:,4),S,...
        'VariableNames',{'Name','ChrNr','Resolution','WN','WS',...
        'binNrStart','binNrEND','START','END','S'});

    ENT3C_OUT=[ENT3C_OUT;OUT1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% similarity table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAMPLES=unique(ENT3C_OUT.Name);
comparisons = get_pairwise_combs(SAMPLES);
Similarity = [];plotted = "tempstring";
for f=1:size(comparisons,1)

    S1 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,1}),:);
    S2 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,2}),:);
    Q = corrcoef(S1.S,S2.S);Q=Q(1,2);

    Similarity=[Similarity;...
        table(comparisons{f,1},comparisons{f,2},Q,...
        'VariableNames',{'Sample1','Sample2','Q'})];

    if ~ismember(comparisons{f,1}, plotted)
        figure(1);plot(S1.S);hold on;
        plotted = [plotted;comparisons{f,1}];
    end
    if ~ismember(comparisons{f,2}, plotted)
        figure(1);plot(S2.S);hold on;
        plotted = [plotted;comparisons{f,2}];
    end
end

L = arrayfun(@(s) setfield(s, 'META', strrep(s.META, '_', ' ')), FNs);

legend(L.META)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);saveas(gcf,sprintf('%s/ENT3C_OUT.png',OUT_DIR))
writetable(ENT3C_OUT,sprintf('%s/ENT3C_OUT.csv',OUT_DIR),'Delimiter','tab')
writetable(ENT3C_OUT,sprintf('%s/ENT3C_similarity.csv',OUT_DIR),'Delimiter','tab')



