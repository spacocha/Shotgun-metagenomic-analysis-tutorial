% Need to have run process_tpm and loaded haotaxat.txt
% and run command [HAOtxnum HAOtx]=sorttaxa(haotaxa,50,indenv+1)
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command


%mk_amos
%can't subtract log values, so use these
nlamoAt=metabolicgenesnoatpHGnorm(konoatp==10944.1,:)';
nlamoBt=metabolicgenesnoatpHGnorm(konoatp==10945.1,:)';
amoABt=log(nlamoAt + nlamoBt);

Cruise=(metadata.Yearday<125)*1+(metadata.Yearday>125&metadata.Yearday<175)*2+(metadata.Yearday>175&metadata.Yearday<225)*3+(metadata.Yearday>225)*4;
Surf=metadata.Depth<2;
Bot=metadata.Depth>2;

%snp_hydro
clf
subplot(121)
scatter(log(nitri1_freference(Surf)),amoABt(Surf),80,Cruise(Surf),'filled');hold on
scatter(log(nitri1_freference(Bot)),amoABt(Bot),80,Cruise(Bot),'filled','d');
xlabel('log ammonia oxidation (mmol m^{-3} day^{-1})', 'FontWeight','bold')
ylabel('log normalized gene abundance', 'FontWeight','bold')
title("(a) {\itamoAB} ({\itr}="+ num2str(corr_obs(24), '%5.2f') + " {\itp}=" + num2str(pval(24), '%5.2f') + ")", 'fontweight','bold','fontsize',12)
%hc=colorbar('Location',['eastoutside']);
grid on


colormap("default")
cmap=colormap;
colormap(cmap([1 86 171 255],:))
hc1=colorbar();
caxis([0.5 4.5])
hc1.Position =[0.04 0.35 0.02 0.3];
xlabel(hc1,'Cruise Month');
labelsc={'Apr','Jun','Jul','Aug'}
hc1.Label.Rotation=0;
hc1.Label.Position=[0.05 0.35 0];
hc1.Direction="reverse";
hc1.Ticks=(1:4);
hc1.TickLabels=labelsc;