function [microC,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution)
% microC ... balanced matrix at desired Resolution
% BINIDs ... start and end binNr of non-zero upper diagonal elements  
% BINS ... BinNr, start and end of all elements
% EXCLUDE ... index of empty bins
% INCLUDE ... index of non-empty bins

% h5disp(FN);
% info = h5info(FN);

if contains(FN,'mcool')
    preStr=sprintf('/resolutions/%d',Resolution);
else
    preStr=[];
end

ChrSizes=ChrInfo(FN,preStr);
% h5disp(FN,sprintf('%s',preStr))
% h5disp(FN,'/resolutions')

% non-zero uppder diagonal elements
BINIDs = [h5read(FN,sprintf('%s/pixels/bin1_id',preStr)) h5read(FN,sprintf('%s/pixels/bin2_id',preStr))];
counts = h5read(FN,sprintf('%s/pixels/count',preStr));

BINS = [h5read(FN,sprintf('%s/bins/start',preStr)) h5read(FN,sprintf('%s/bins/end',preStr))];
BINS = [(1:size(BINS,1))',BINS];
weights = h5read(FN,sprintf('%s/bins/weight',preStr));

chrs = string([h5read(FN,sprintf('%s/bins/chrom',preStr))]);
BIN_TABLE = table(chrs,BINS(:,1),BINS(:,2),BINS(:,3),weights,...
    'VariableNames',{'chrs','binNrALL','START','END','weights'});
BT=table();binNrCHRS=[];
for c=1:length(ChrNr)
    BT = [BT;BIN_TABLE(BIN_TABLE.chrs==string(['chr',num2str(ChrNr(c))]),:)];
    binNrCHRS=[binNrCHRS;transpose(1:sum(BIN_TABLE.chrs==string(['chr',num2str(ChrNr(c))])))];
end
BIN_TABLE=BT;clear BT

BIN_TABLE = addvars(BIN_TABLE,binNrCHRS,'After','binNrALL');

f=find(BINIDs(:,1)>=min(BIN_TABLE.binNrALL(:))&BINIDs(:,2)>=min(BIN_TABLE.binNrALL(:))&...
    BINIDs(:,1)<=max(BIN_TABLE.binNrALL(:))&BINIDs(:,2)<=max(BIN_TABLE.binNrALL(:)));

BINIDs=BINIDs(f,:);
counts=counts(f,:);
BINIDs=BINIDs-BINIDs(1,1)+1;
%%%
microC=nan(size(BIN_TABLE,1),size(BIN_TABLE,1));
index = sub2ind(size(microC),BINIDs(:,1),BINIDs(:,2));
microC(index)=counts;
microC=microC.*BIN_TABLE.weights;

microC=triu(microC)+triu(microC,1)';

INCLUDE = unique([find(sum(isnan(microC),2)<size(microC,2))]);

EMPTY_BIN=nan(size(BIN_TABLE,1),1);
EMPTY_BIN(INCLUDE)=1;
EMPTY_BIN(isnan(BIN_TABLE.weights))=nan;

BIN_TABLE = addvars(BIN_TABLE,EMPTY_BIN,'After','binNrALL','NewVariableNames',{'CONTACT'});
