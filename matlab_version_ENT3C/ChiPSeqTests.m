%TEST=sprintf('NaN_ResolutionTest.all_cells.noNorm.WLMAX1000.06.12.2023.mat');
%TEST=sprintf('40kb.fixedChrSplit.7.ALLINTs.NoNorm.WLMAX1000.30.11.2023.mat')
%FN_RES=sprintf('../../MAT/VN_ENT_%s',TEST);
%ENT3C = readtable(FN_RES);

TEST=sprintf('OUTPUT/MATLAB/HFFc6_normalized_pooled_ENT3C_OUT.csv')
FN_RES=TEST

%cooltools insulation -o /HDD1/DocumentsHDD1/DISSERTATION/07_2022-/scripts/3C_SIMILARITY/ENT3C/local_git_31_01_24/ENT3C_local_git/insulation_score/insulation.score.HFF.cooltools.2.txt 4DNESWST3UBH_HFFcombined.5000bp.cool 50000



%%
close all
ENT3C = readtable(FN_RES);

% ENT3C = ENT3C(ENT3C.Resolution==40e3&ENT3C.sub_m_dim==500,:)
% c3=repmat(sprintf('chr%d',3),sum(ENT3C.ChrNr==3),1);
% c6=repmat(sprintf('chr%d',6),sum(ENT3C.ChrNr==6),1);
% ENT3C.ChrNr=[c3;c6];
% ENT3C = ENT3C(:,[2,9,10])
%writetable(ENT3C,'/HDD3/DocumentsHDD3/07_2022-/Mapping/4DNESWST3UBH_HFFc6/PreProcessed/ENT3CBED_regions_insulation.csv','Delimiter','tab')


%ENT3C = ENT3C(isnan(ENT3C.NR_INT),:);
DIR='/HDD3/DocumentsHDD3/07_2022-';
Resolutions=unique(ENT3C.Resolution)';%Resolutions=Resolutions(3:5);

ChrNrs=unique(ENT3C.ChrNr)';%ChrNrs=ChrNrs(2:2:end);


for C={'HFF'}%unique(ENT3C.Name)'
    n=n+1;
    C=C{1};

    IS = readtable('insulation_score/insulation.score.HFF.cooltools.200000.txt');

    counts = readtable('RNASeqData/PC_CODING_GENE_COUNTS_HFF.txt');

    chip=sprintf('H3K4me3/%s.txt',C);
    chip=readtable(chip);

    % H=chip(strcmp(chip.Var1,'chr3'),:);
    % figure;plot(H.Var4);
    % H.Var4= rescale(H.Var4,0,1);

    for R=5e3%Resolutions
        for ChrNr=ChrNrs
            for sub_n=[500 1e3]%unique(ENT3C.sub_m_dim)'
                figure(1);nexttile;

                %S=ENT3C(ENT3C.ChrNr==ChrNr & strcmp(C,ENT3C.Name) & ENT3C.Resolution==R & ENT3C.BioRep==1,:);
                S=ENT3C(ENT3C.ChrNr==ChrNr & strcmp([C,'c6'],ENT3C.Name{1}) & ENT3C.Resolution==R & ENT3C.sub_m_dim==sub_n,:);

                % H=counts;
                S_H=nan(length(S.S),5);
                % for k=1:length(S.S)
                %     if length(H.HFF(H.Transcription_start_site__TSS_>=S.START(k,:)&H.Gene_end__bp_<=S.END(k,:)))~=0
                %         S_H(k,1:3)=[S.S(k),...
                %             nanmean(H.HFF(H.Transcription_start_site__TSS_>=S.START(k,:)&H.Gene_end__bp_<=S.END(k,:))),...
                %             nanvar(H.HFF(H.Transcription_start_site__TSS_>=S.START(k,:)&H.Gene_end__bp_<=S.END(k,:)))];
                %     else
                %         S_H(k,2:3)=0;
                %     end
                % end
                % H=chip;
                % for k=1:length(S.S)
                %     if length(H.Var4(H.Var2>=S.START(k,:)&H.Var3<=S.END(k,:)))~=0
                %         S_H(k,4:5)=[nanmean(H.Var4(H.Var2>=S.START(k,:)&H.Var3<=S.END(k,:))),...
                %             nanmedian(H.Var4(H.Var2>=S.START(k,:)&H.Var3<=S.END(k,:)))];
                %     else
                %         S_H(k,4:END)=0;
                %     end
                % end
                H=IS;S_IS=[];S_BS=[];

                for k=1:length(S.S)
                    REG=(H.start>=S.START(k,:)&H.xEnd<=S.END(k,:));
                    if sum(REG==0)~=length(REG)
                        S_H(k,:)=[S.S(k),nanmean(H.boundary_strength_200000(REG)),...
                            nanmedian(H.boundary_strength_200000(REG)),...
                            nanmean(H.log2_insulation_score_200000(REG)),...
                            nanmedian(H.log2_insulation_score_200000(REG))];
                        S_IS=[S_IS;S.START(k,:),S.END(k,:),nanmean(H.boundary_strength_200000(REG))];
                        S_BS=[S_BS;S.START(k,:),S.END(k,:),nanmean(H.log2_insulation_score_200000(REG))];
                    else
                        S_H(k,2:end)=0;
                        S_IS=[S_IS;S.START(k,:),S.END(k,:),0];
                        S_BS=[S_BS;S.START(k,:),S.END(k,:),0];
                    end
                end

writematrix(S_IS,'insulation_score/insulation_score.csv',"Delimiter","\t")
writematrix(S_BS,'insulation_score/boundary_strength.csv',"Delimiter","\t")

                for k=1:5
                    S_H(:,k)=movmean(S_H(:,k),3);
                    S_H(:,k)=rescale(S_H(:,k),0,1);
                end
                %subplot(3,1,1);

                % plot(S_H(:,2),'b','LineWidth',1);axis tight;hold on;grid on;
                % plot(S_H(:,3),'c','LineWidth',1);axis tight;hold on;grid on;
                % plot(S_H(:,4),'r','LineWidth',1);axis tight;hold on;grid on;
                % plot(S_H(:,5),'y','LineWidth',1);axis tight;hold on;grid on;

                plot(S_H(:,2),'b','LineWidth',1);axis tight;hold on;grid on;
                plot(S_H(:,4),'r','LineWidth',1);axis tight;hold on;grid on
                plot(S_H(:,1),'k--','LineWidth',2);axis tight;hold on;grid on;

                if R==40e3 && ChrNr==6
                    legend('mean BR','mean IS','ENT3C')
                end

                S_H(isnan(S_H(:,1)),:)=[];
                rr_BS=corrcoef(S_H(:,1),S_H(:,2));
                rr_IS=corrcoef(S_H(:,1),S_H(:,4));
                r=[rr_BS(1,2) rr_IS(1,2)];
                set(gca,'Xtick',[],'Ytick',[])
                title(sprintf('Chr%d %dkb $n=%d$ $r_{BS}=%4.2f$ $r_{IS}=%4.2f$',ChrNr,R/1e3,sub_n,r),...
                    'interpreter','latex','fontsize',9)

            end
        end
    end
end

figure(1);
set(gcf,'Position',[1921 69 1920 1013])
saveas(gcf,sprintf('insulation_score/correlations.means.200000.svg'))

