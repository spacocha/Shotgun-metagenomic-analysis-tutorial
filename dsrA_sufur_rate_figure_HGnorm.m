% Need to have run process_tpm and loaded haotaxat.txt
% run process_TPM_chesapeake_SPP, ches_freference_SPP, relative_rates_figure command

Surf=metadata.Depth<2;
Bot=metadata.Depth>2;

%mk_adsrox
%can't subtract log values, so use these
nldsrAoxt=metabolicgenesnoatpHGnorm(konoatp==11180.1,:)';
nldsrBoxt=metabolicgenesnoatpHGnorm(konoatp==11181.1,:)';
dsrABoxt=log(nldsrAoxt + nldsrBoxt);

%mk_dsrred
%can't subtract log values, so use these
nldsrAredt=metabolicgenesnoatpHGnorm(konoatp==11180.2,:)';
nldsrBredt=metabolicgenesnoatpHGnorm(konoatp==11181.2,:)';
dsrABredt=log(nldsrAredt + nldsrBredt);

%snp_hydro
clf
colormap hot

subplot(121)
scatter(log((soxo_freference(Surf)+soxno3_freference(Surf)+soxno2_freference(Surf))),dsrABoxt(Surf),40,DO_freference(Surf),'filled');hold on;
scatter(log((soxo_freference(Bot)+soxno3_freference(Bot)+soxno2_freference(Bot))),dsrABoxt(Bot),40,DO_freference(Bot),'filled','d');
xlabel('log sulfur oxidation (mmol m^-^3 d^-^1)')
ylabel('log normalized gene abundance')
title("(a) Ox. {\itdsrAB} ({\itr}="+ num2str(corr_obs(31), '%5.2f') + ", {\itp}=" + num2str(pval(31), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on


subplot(122)
scatter(log((srra_freference(Surf))'),dsrABredt(Surf),40,DO_freference(Surf),'filled');hold on;
scatter(log((srra_freference(Bot))'),dsrABredt(Bot),40,DO_freference(Bot),'filled','d');
xlabel('log sulfate reduction. (mmol m^-^3 d^-^1)')
ylabel('log normalized gene abundance')
title("(b) Red. {\itdsrAB} ({\itr}="+ num2str(corr_obs(32), '%5.2f') + ", {\itp}=" + num2str(pval(32), '%5.2f') + ")")
%hc=colorbar('Location',['eastoutside']);
grid on



hc1=colorbar;
hc1.Position =[0.05 0.35 0.02 0.3];
ylabel(hc1,'oxygen (\mum)');