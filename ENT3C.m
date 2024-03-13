clear all
addpath('MATLAB_functions/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load json file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = fileread('config/config.json');
config = jsondecode(config);

WN_MAX=config.WN_MAX;
CHRSPLIT=config.CHRSPLIT;
SUB_M_SIZE_FIX=config.SUB_M_SIZE_FIX;
ChrNr=config.ChrNr;
Resolution=config.Resolution;
WS=config.WS;
NormM=config.NormM;
DATA_PATH=config.DATA_PATH;
FILES = config.FILES;
FILES = reshape(FILES,2,size(FILES,1)/2);
FILES = [cellfun(@(file) fullfile(DATA_PATH, file), FILES(1,:), 'UniformOutput', false); FILES(2,:)];
OUT_DIR=config.OUT_DIR;
OUT_DIR=[OUT_DIR,'MATLAB'];
OUT_PREFIX=config.OUT_PREFIX;
if isempty(OUT_PREFIX)
    OUT_PREFIX=sprintf('Chr%d_%dkb',ChrNr,Resolution/1e3);
end
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INFO = struct('FN', '', 'META', '');
FNs = [];
for f = 1:size(FILES,2)
    FNs = [FNs, INFO];
    FNs(end).FN = FILES{1,f};
    FNs(end).META = FILES{2,f};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract common empty bins of input matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXCLUDE=[];
for f = 1:numel(FNs)
        FN=FNs(f).FN;
        [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM);
        EXCLUDE=[EXCLUDE;BIN_TABLE.binNr(isnan(BIN_TABLE.CONTACT))];
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

    OUT1 = table(repmat({FNs(f).META},N,1),...
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
        figure(1);plot(S1.S);hold on;axis tight
        plotted = [plotted;comparisons{f,1}];
    end
    if ~ismember(comparisons{f,2}, plotted)
        figure(1);plot(S2.S);hold on;axis tight
        plotted = [plotted;comparisons{f,2}];
    end


end
L = arrayfun(@(s) setfield(s, 'META', strrep(s.META, '_', ' ')), FNs);
legend(L.META)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);saveas(gcf,sprintf('%s/%s_ENT3C_signals.png',OUT_DIR,OUT_PREFIX))
writetable(ENT3C_OUT,sprintf('%s/%s_ENT3C_OUT.csv',OUT_DIR,OUT_PREFIX),'Delimiter','tab')
writetable(Similarity,sprintf('%s/%s_ENT3C_similarity.csv',OUT_DIR,OUT_PREFIX),'Delimiter','tab')





