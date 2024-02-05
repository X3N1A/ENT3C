function [ChrSizes]=ChrInfo(FN,preStr)

if strcmp(FN(end-1:end),'h5')==1
    chroms = [deblank(h5read(FN,sprintf('%s',preStr)))];
    chromsN = [deblank(strrep(h5read(FN,sprintf('%s',preStr)),'chr',''))];
    RFN='/HDD1/DocumentsHDD1/REF_FILES/genome/Human/hg19/hg19.chrom.sizes';
    L=nan(size(chromsN,1),1);
    for k=1:size(chromsN,1)
        [~,l]=unix(sprintf("cat %s | grep -Po 'chr%s\t[0-9]+' | awk '{print $2}' ",RFN,chromsN(k)));
        L(k)=str2num(l);
    end
    ChrSizes = table(chroms,L);
else
    chroms = [deblank(strrep(h5read(FN,sprintf('%s/chroms/name',preStr)),'chr',''))];
    L = h5read(FN,sprintf('%s/chroms/length',preStr));
    ChrSizes = table(chroms,L);
end

% ChrSize=nan(length(ChrNrs),2);
% 
% for k=1:length(ChrNrs)
%     [~,CS]=unix(sprintf("cat /HDD1/DocumentsHDD1/REF_FILES/genome/Human/38p13/hg38.fa.sizes | grep -v 'chrX'| grep -v 'chrY' | grep -E 'chr%d\t' | awk '{print $2}'",ChrNrs(k)));
%     ChrSize(k,:)=[k str2double(CS)];
% end


% if strcmp(ChrNr,'ALL')==1
%     [~,ChrSize]=unix(sprintf("cat /HDD1/DocumentsHDD1/REF_FILES/genome/Human/38p13/hg38.fa.sizes | grep -v 'chrX'| grep -v 'chrY' | awk '{print $2}'"));
%     ChrSize=str2num(ChrSize);
%     ChrSize=sum(ChrSize);
% else
%     [~,ChrSize]=unix(sprintf("cat /HDD1/DocumentsHDD1/REF_FILES/genome/Human/38p13/hg38.fa.sizes | grep -v 'chrX'| grep -v 'chrY' | grep -E 'chr%d\t' | awk '{print $2}'",ChrNr));
%     ChrSize=str2num(ChrSize);
% end
