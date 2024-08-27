function [microC,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,Norm,weights_name)
% microC ... balanced matrix at desired Resolution
% BINIDs ... start and end binNr of non-zero upper diagonal elements  
% BINS ... BinNr, start and end of all elements
% EXCLUDE ... index of empty bins
% INCLUDE ... index of non-empty bins
% weights_name ... name of weights vector in cooler 

% h5disp(FN);

if contains(FN,'mcool')
    preStr=sprintf('/resolutions/%d',Resolution);
else
    preStr=[];
end

% h5disp(FN,sprintf('%s',preStr))
% h5disp(FN,'/resolutions')

% non-zero uppder diagonal elements
BINIDs = [h5read(FN,sprintf('%s/pixels/bin1_id',preStr)) h5read(FN,sprintf('%s/pixels/bin2_id',preStr))];
counts = h5read(FN,sprintf('%s/pixels/count',preStr));

BINS = [h5read(FN,sprintf('%s/bins/start',preStr)) h5read(FN,sprintf('%s/bins/end',preStr))];
BINS = [(1:size(BINS,1))',BINS];
if Norm==1
    weights = h5read(FN,sprintf('%s/bins/%s',preStr,weights_name));
else
    weights = nan(size(BINS,1),1);
end
chrs = string([h5read(FN,sprintf('%s/bins/chrom',preStr))]);

BIN_TABLE = table(chrs,BINS(:,1),BINS(:,2),BINS(:,3),weights,...
    'VariableNames',{'chrs','BINS_ALL','START','END','weights'});
binNr=[];

BIN_TABLE=BIN_TABLE((strcmp(BIN_TABLE.chrs,ChrNr{1}) | strcmp(BIN_TABLE.chrs,['chr',ChrNr{1}])),:);
binNr=[binNr;transpose(1:sum((strcmp(BIN_TABLE.chrs,ChrNr{1}) | strcmp(BIN_TABLE.chrs,['chr',ChrNr{1}]))))];

BIN_TABLE = addvars(BIN_TABLE,binNr,'After','BINS_ALL');

%get individual chromosome
f=find(BINIDs(:,1)>=min(BIN_TABLE.BINS_ALL(:))&BINIDs(:,2)>=min(BIN_TABLE.BINS_ALL(:))&...
    BINIDs(:,1)<=max(BIN_TABLE.BINS_ALL(:))&BINIDs(:,2)<=max(BIN_TABLE.BINS_ALL(:)));

BINIDs=BINIDs(f,:);
counts=counts(f,:);
BINIDs=BINIDs-BINIDs(1,1)+1;
%%%
microC=nan(size(BIN_TABLE,1),size(BIN_TABLE,1));
index = sub2ind(size(microC),BINIDs(:,1),BINIDs(:,2));
microC(index)=counts;
if Norm==1
    microC=microC.*BIN_TABLE.weights;
end

microC=triu(microC)+triu(microC,1)';

INCLUDE = unique([find(sum(isnan(microC),2)<size(microC,2))]);

EMPTY_BIN=nan(size(BIN_TABLE,1),1);
EMPTY_BIN(INCLUDE)=1;

BIN_TABLE = addvars(BIN_TABLE,EMPTY_BIN,'After','BINS_ALL','NewVariableNames',{'CONTACT'});
