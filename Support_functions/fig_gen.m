load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_adj_R12_nobg.mat')
ADJ = imFinal;
load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_CS_R12_nobg.mat')
CS = imFinal;
load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_SENSE_R12_nobg.mat')
SENSE = imFinal;
%% compare three
% imshow(imcomplement(cat(2,imcrop(mean(abs(CS(:,:,:,:)),4),[80 67 200 220]),imcrop(mean(abs(SENSE(:,:,:,:)),4),[80 67 200 220]))),[])
imshow(imcomplement(cat(2,imcrop(mean(abs(CS(:,:,:,:)),4),[80 67 200 220]),0.8*imcrop(mean(abs(ADJ(:,:,:,:)),4),[80 67 200 220])-0.000002,imcrop(mean(abs(SENSE(:,:,:,:)),4),[80 67 200 220])-0.000001)),[0.99999 1])
F = getframe(gca);
imwrite(F.cdata, ['temp_mean_R12_whole.png']);

imshow(imcomplement(cat(2,imcrop(mean(abs(CS(:,:,:,:)),4),[210 180 50 50]),0.8*imcrop(mean(abs(ADJ(:,:,:,:)),4),[210 180 50 50])-0.000002,imcrop(mean(abs(SENSE(:,:,:,:)),4),[210 180 50 50])-0.000001)),[0.99999 1])
F = getframe(gca);
imwrite(F.cdata, ['temp_mean_R12_detail.png']);

%% time courses
imshow(imcomplement(0.8*cat(2,imcrop(mean(abs(ADJ(:,:,:,2)),4),[80 67 200 220]),imcrop(mean(abs(ADJ(:,:,:,5)),4),[80 67 200 220]),imcrop(mean(abs(ADJ(:,:,:,8)),4),[80 67 200 220])-0.000002)),[0.99998 1])
drawnow
F = getframe(gca);
imwrite(F.cdata, ['R12_ADJ_time.png']);
imshow(imcomplement(cat(2,imcrop(mean(abs(SENSE(:,:,:,2)),4),[80 67 200 220]),imcrop(mean(abs(SENSE(:,:,:,5)),4),[80 67 200 220]),imcrop(mean(abs(SENSE(:,:,:,8)),4),[80 67 200 220])-0.000001)),[0.99998 1])
drawnow
F = getframe(gca);
imwrite(F.cdata, ['R12_SENSE_time.png']);
imshow(imcomplement(cat(2,imcrop(mean(abs(CS(:,:,:,2)),4),[80 67 200 220]),imcrop(mean(abs(CS(:,:,:,5)),4),[80 67 200 220]),imcrop(mean(abs(CS(:,:,:,8)),4),[80 67 200 220]))),[0.99998 1])
F = getframe(gca);
imwrite(F.cdata, ['R12_CS_time.png']);

%% acceleration factors
load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_CS_R23_nobg.mat')
CS_23 = imFinal;
load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_CS_R46_nobg.mat')
CS_46 = imFinal;
load('jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_offlineRecon_CS_R1_nobg.mat')
CS_12 = CS;
CS_1 = imFinal;

imshow(imcomplement(cat(2,imcrop(mean(abs(CS_1(:,:,:,:)),4),[80 67 200 220]),imcrop(mean(abs(CS_12(:,:,:,:)),4),[80 67 200 220]),imcrop(mean(abs(CS_23(:,:,:,:)),4),[80 67 200 220]),imcrop(mean(abs(CS_46(:,:,:,:)),4),[80 67 200 220]))),[0.99999 1])
F = getframe(gca);
imwrite(F.cdata, ['CS_accelerated.png']);

