function save_configurations(PHI_MAX,CHRSPLIT,SUB_M_SIZE_FIX,Resolutions,ChrNrs,...
    phi,NormM,DATA_PATH,FILES,OUT_DIR,OUT_PREFIX)

fid=fopen(sprintf('%s/CONFIG.txt',OUT_DIR),'w');
fprintf('Resolutions=%s, Chrs=%s\nPHI_MAX=%d phi=%d c=%d n=%d normalize=%d\n',Resolutions,ChrNrs,PHI_MAX,phi,CHRSPLIT,SUB_M_SIZE_FIX,NormM)
fprintf('DATA_PATH=%s\n',DATA_PATH)
fprintf('FILES=%s \n',FILES{:})
fprintf('OUT=%s/%s\n',OUT_DIR,OUT_PREFIX)
fclose(fid)


end