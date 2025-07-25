function [] = ENT3C(config_file)
addpath('MATLAB_functions/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load json file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%config_file='config/config.json';
config = fileread(config_file);
config = jsondecode(config);
PHI_MAX=config.PHI_MAX;
CHRSPLIT=config.CHRSPLIT;
SUB_M_SIZE_FIX=config.SUB_M_SIZE_FIX;
Resolutions=str2num(config.Resolution);
phi=config.phi;
NormM=config.NormM;
weights_name=config.WEIGHTS_NAME;
DATA_PATH=config.DATA_PATH;
FILES = config.FILES;
FILES = reshape(FILES,2,size(FILES,1)/2);
FILES = [cellfun(@(file) fullfile(DATA_PATH, file), FILES(1,:), 'UniformOutput', false); FILES(2,:)];
if any(contains(FILES(2,:),'_'))% only important in case of computing mean(Q_BRs) mean(Q_nonBRs)
    Biologicap_replicates=true;
else
    Biologicap_replicates=false;
end
OUT_DIR=config.OUT_DIR;
OUT_DIR=[OUT_DIR,'/MATLAB'];
OUT_PREFIX=config.OUT_PREFIX;
if isempty(OUT_PREFIX)
    OUT_PREFIX=sprintf('%dkb',Resolution/1e3);
end
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

ChrNrs=config.ChrNr;
ChrNrs=strsplit(ChrNrs, ',');

ENT3C_OUT_FN=sprintf('%s/%s_ENT3C_OUT.csv',OUT_DIR,OUT_PREFIX);

if ~exist(ENT3C_OUT_FN,'file')
    sprintf("Generating new file: %s",ENT3C_OUT_FN)
    ENT3C_OUT=table('Size',[0 11],...
    'VariableTypes', {'cell','cell','int64','int64','int64','int64','int64','int64','int64','int64','double'},...
    'VariableNames', {'Name','ChrNr','Resolution','n','PHI','phi','binNrStart','binNrEND','START','END','S'});
    for Resolution=Resolutions
        for ChrNr=ChrNrs

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
            CM=jet(size(FNs,2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract common empty bins of input matrices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            EXCLUDE=[];
            for f = 1:numel(FNs)
                FN=FNs(f).FN;
                [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM,weights_name);
                if NormM==0
                    EXCLUDE=[EXCLUDE;BIN_TABLE.binNr(isnan(BIN_TABLE.CONTACT))];
                elseif NormM==1
                    EXCLUDE=[EXCLUDE;BIN_TABLE.binNr(isnan(BIN_TABLE.CONTACT));BIN_TABLE.binNr(isnan(BIN_TABLE.weights))];
                end
            end
            EXCLUDE=unique(EXCLUDE);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fill ENT3C table
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure(1);tiledlayout(2,2,'tilespacing','compact');
            for f = 1:numel(FNs)
                FN=FNs(f).FN;
                [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM,weights_name);
                INCLUDE = 1:size(M,1);
                INCLUDE = setdiff(INCLUDE,EXCLUDE);

                M = M(INCLUDE,INCLUDE);

                %{
            PERC_0 = 100*sum(isnan(M(:)))/numel(M(:))
            nexttile
            imagesc(M);imagesc(log(M));axis square;set(gca,'YDir','normal');title(sprintf('Chr%s %s NormM%d %d%%',...
                    ChrNr{1},FNs(f).META,NormM,round(100*sum(isnan(M(:)))/numel(M(:)))))
                %}
                BIN_TABLE = BIN_TABLE(INCLUDE,:);

                % writetable(BIN_TABLE,sprintf('%s_chr%d_BINMATRIXMATLAB.csv',FN,ChrNr),'Delimiter','tab')

                [S, SUB_M_SIZE1, PHI_1, phi_1, BIN_TABLE_NEW] = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,PHI_MAX,phi,BIN_TABLE);
                N = length(S);

                OUT1 = table(repmat({FNs(f).META},N,1),...
                    repmat({sprintf('%s',ChrNr{1})},N,1),repmat(Resolution,N,1),...
                    repmat(SUB_M_SIZE1,N,1),repmat(PHI_1,N,1),repmat(phi_1,N,1),...
                    BIN_TABLE_NEW(:,1),BIN_TABLE_NEW(:,2),...
                    BIN_TABLE_NEW(:,3),BIN_TABLE_NEW(:,4),S,...
                    'VariableNames',{'Name','ChrNr','Resolution','n','PHI','phi',...
                    'binNrStart','binNrEND','START','END','S'});

                ENT3C_OUT=[ENT3C_OUT;OUT1];
              
            end
        end
    end

    writetable(ENT3C_OUT,ENT3C_OUT_FN,'Delimiter','tab')

else

    opts = detectImportOptions(ENT3C_OUT_FN, 'TextType', 'string');
    opts = setvartype(opts, 'ChrNr', 'string'); % or 'char'
    ENT3C_OUT=readtable(sprintf('%s',ENT3C_OUT_FN),opts);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill similarity table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAMPLES=unique(ENT3C_OUT.Name);

Similarity = [];
if numel(SAMPLES)>1
    comparisons = get_pairwise_combs(SAMPLES);
    for Resolution=Resolutions
        figure(Resolution);hold on;tiledlayout('flow');hold on%'flow');hold on
        for ChrNr=ChrNrs

            figure(Resolution);nexttile

            plotted = "tempstring";
            c=1;
            for f=1:size(comparisons,1)

                S1 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,1})&strcmp(ENT3C_OUT.ChrNr,ChrNr)&ENT3C_OUT.Resolution==Resolution,:);
                S2 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,2})&strcmp(ENT3C_OUT.ChrNr,ChrNr)&ENT3C_OUT.Resolution==Resolution,:);
                non_nan_idx = ~isnan(S1.S)&~isnan(S2.S);
                Q = corrcoef(S1.S(non_nan_idx),S2.S(non_nan_idx));Q=Q(1,2);

                Similarity=[Similarity;...
                    table(Resolution,ChrNr,comparisons{f,1},comparisons{f,2},Q,...
                    'VariableNames',{'Resolution','ChrNr','Sample1','Sample2','Q'})];

                if ~ismember(comparisons{f,1}, plotted)
                    figure(Resolution);plot(S1.S,'color',CM(c,:),'linewidth',1.3);hold on;axis tight
                    plotted = [plotted;comparisons{f,1}];c=c+1;
                end
                if ~ismember(comparisons{f,2}, plotted)
                    figure(Resolution);plot(S2.S,'color',CM(c,:),'linewidth',1.3);hold on;axis tight
                    plotted = [plotted;comparisons{f,2}];c=c+1;
                end
            end
            if Biologicap_replicates
                title(sprintf('Chr%s %dkb\n$\\overline{Q}_{BR}=%4.2f$ $\\overline{Q}_{nonBR}=%4.2f$',...
                    ChrNr{1},Resolution/1e3,...
                    mean(Similarity.Q(strcmp(Similarity.ChrNr,ChrNr{1})&Similarity.Resolution==Resolution&...
                    strcmp(extractBefore(Similarity.Sample1,'_'),extractBefore(Similarity.Sample2,'_')))),...
                    mean(Similarity.Q(strcmp(Similarity.ChrNr,ChrNr{1})&Similarity.Resolution==Resolution&...
                    ~strcmp(extractBefore(Similarity.Sample1,'_'),extractBefore(Similarity.Sample2,'_'))))),...
                    'interpreter','latex','fontsize',25)
            else
                title(sprintf('Chr%s %dkb',ChrNr{1},Resolution/1e3),'interpreter','latex','fontsize',25)
            end

            last=ChrNrs(end);
            if strcmp(ChrNr{1},last{1})
                plotted(1)=[];plotted=strrep(plotted,'_',' ');
                legend(plotted,'location','northeastoutside')
            end
            set(gca,'FontSize',15)
        end

        figure(Resolution);
        set(gcf,'Position',[52 215 1717 813])
        saveas(gcf,sprintf('%s/%s_%d_ENT3C_signals.svg',OUT_DIR,OUT_PREFIX,Resolution))

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    writetable(Similarity,sprintf('%s/%s_ENT3C_similarity.csv',OUT_DIR,OUT_PREFIX),'Delimiter','tab')
end

