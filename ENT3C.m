clear all;close all
addpath('MATLAB_functions/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load json file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = fileread('config/config.json');
config = jsondecode(config);
PHI_MAX=config.PHI_MAX;
CHRSPLIT=config.CHRSPLIT;
SUB_M_SIZE_FIX=config.SUB_M_SIZE_FIX;
Resolutions=str2num(config.Resolution);
phi=config.phi;
NormM=config.NormM;
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
OUT_DIR=[OUT_DIR,'MATLAB'];
OUT_PREFIX=config.OUT_PREFIX;
if isempty(OUT_PREFIX)
    OUT_PREFIX=sprintf('%dkb',Resolution/1e3);
end
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end
ChrNrs=str2num(config.ChrNr);

ENT3C_OUT=[];
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
            [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM);
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
        for f = 1:numel(FNs)
            FN=FNs(f).FN;

            [M,BIN_TABLE]=load_cooler(FN,ChrNr,Resolution,NormM);
            INCLUDE = 1:size(M,1);
            INCLUDE = setdiff(INCLUDE,EXCLUDE);

            M = M(INCLUDE,INCLUDE);
            BIN_TABLE = BIN_TABLE(INCLUDE,:);

            writetable(BIN_TABLE,sprintf('%s_chr%d_BINMATRIXMATLAB.csv',FN,ChrNr),'Delimiter','tab')

            [S, SUB_M_SIZE1, WN1, WS1, BIN_TABLE_NEW] = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,PHI_MAX,phi,BIN_TABLE);
            N = length(S);

            OUT1 = table(repmat({FNs(f).META},N,1),...
                repmat(ChrNr,N,1),repmat(Resolution,N,1),...
                repmat(SUB_M_SIZE1,N,1),repmat(WN1,N,1),repmat(WS1,N,1),...
                BIN_TABLE_NEW(:,1),BIN_TABLE_NEW(:,2),...
                BIN_TABLE_NEW(:,3),BIN_TABLE_NEW(:,4),S,...
                'VariableNames',{'Name','ChrNr','Resolution','sub_m_dim','WN','phi',...
                'binNrStart','binNrEND','START','END','S'});

            ENT3C_OUT=[ENT3C_OUT;OUT1];
            ENT3C_OUT
        end
    end
end

writetable(ENT3C_OUT,sprintf('%s/%s_ENT3C_OUT.csv',OUT_DIR,OUT_PREFIX),'Delimiter','tab')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill similarity table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAMPLES=unique(ENT3C_OUT.Name);

Similarity = [];
if numel(SAMPLES)>1
    for Resolution=Resolutions
        figure(Resolution);hold on;tiledlayout('flow');hold on
        for ChrNr=ChrNrs
            figure(Resolution);nexttile

            comparisons = get_pairwise_combs(SAMPLES);

            plotted = "tempstring";
            c=1;
            for f=1:size(comparisons,1)

                S1 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,1})&ENT3C_OUT.ChrNr==ChrNr&ENT3C_OUT.Resolution==Resolution,:);
                S2 = ENT3C_OUT(strcmp(ENT3C_OUT.Name,comparisons{f,2})&ENT3C_OUT.ChrNr==ChrNr&ENT3C_OUT.Resolution==Resolution,:);
                Q = corrcoef(S1.S,S2.S);Q=Q(1,2);

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
                title(sprintf('Chr%d %dkb\n$\\overline{Q}_{BR}=%4.2f$ $\\overline{Q}_{nonBR}=%4.2f$',...
                    ChrNr,Resolution/1e3,...
                    mean(Similarity.Q(Similarity.ChrNr==ChrNr&Similarity.Resolution==Resolution&...
                    strcmp(extractBefore(Similarity.Sample1,'_'),extractBefore(Similarity.Sample2,'_')))),...
                    mean(Similarity.Q(Similarity.ChrNr==ChrNr&Similarity.Resolution==Resolution&...
                    ~strcmp(extractBefore(Similarity.Sample1,'_'),extractBefore(Similarity.Sample2,'_'))))),...
                    'interpreter','latex','fontsize',10)
            else
                title(sprintf('Chr%d %dkb',ChrNr,Resolution/1e3),'interpreter','latex','fontsize',9)
            end

            if ChrNr==ChrNrs(end)
                plotted(1)=[];plotted=strrep(plotted,'_',' ');
                legend(plotted,'location','northeastoutside')
            end
            set(gca,'FontSize',9)
            plotted
        end

        figure(Resolution);
        set(gcf,'Position',[52 215 1717 813])
        saveas(gcf,sprintf('%s/%s_ENT3C_signals.svg',OUT_DIR,OUT_PREFIX))

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    writetable(Similarity,sprintf('%s/%s_ENT3C_similarity.csv',OUT_DIR,OUT_PREFIX),'Delimiter','tab')
end
