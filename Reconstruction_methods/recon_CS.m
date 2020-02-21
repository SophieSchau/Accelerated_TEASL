clear

%% 0. Add paths and setup FesslerLibrary adds

% irtdir = '/Users/jwoods/Documents/ASL/Recon/FesslerLibrary/irt/';
% run([irtdir 'setup.m'])

%% 1. Setup
% File ID

folder = '/Users/schauman/Documents/ONBI/DoctoralProject/RawData/Angio_Joe';

%%%%%% Seq %%%%%

% MeasFileID = [folder '/meas_MID276_jw_CAPIASL_CV_nce_angio_SEQ_LD360_9Phase_VFA_HighRes_FID328.dat'];
% phases = 1:9; % Number of time points
% numSpokes = 1104;
% averages = 1;
% TE = false;
% encMat = [-1; 1];
% 
% R = 12;
% acc_factor = R/2;
% spokes(:,1) = ceil(1:R:numSpokes);
% spokes(:,2) = ceil(R/3:R:numSpokes);
% spokes(:,3) = ceil(2*R/3:R:numSpokes);
% spokes(:,4) = ceil(1:R:numSpokes);
% spokes(:,5) = ceil(R/3:R:numSpokes);
% spokes(:,6) = ceil(2*R/3:R:numSpokes);
% spokes(:,7) = ceil(1:R:numSpokes);
% spokes(:,8) = ceil(R/3:R:numSpokes);
% spokes(:,9) = ceil(2*R/3:R:numSpokes);


% %%%%%%%% TE %%%%%%%%

MeasFileID = [folder '/meas_MID277_jw_CAPIASL_CV_nce_angio_TE3_LD360_3Phase_VFA_highRes_FID329.dat'];
phases = 1:3; % Number of time points
numSpokes = 552;
averages = 1;
TE = true;
% encMat = normHadamard(4,averages);
encMat = [-1 -1 -1;...
1  -1  1;...
-1 1  1;...
1  1 -1];

R = 12;
acc_factor = R;
spokes(:,1) = ceil(1:R:numSpokes);
spokes(:,2) = ceil(R/3:R:numSpokes);
spokes(:,3) = ceil(2*R/3:R:numSpokes);

%% List of acquired angles for non-golden ratio
angles = zeros( numSpokes/R, length(phases) );
for ph = phases
    lineNum = spokes(:,ph)-1;
    angles(:,ph) = pi * lineNum / (numSpokes);
    
    if mod(numSpokes,2) == 1
        angles(:,ph) = angles(:,ph) * 2;
    end
    % figure; hold on
    polarplot(angles(:,ph),ones(size(angles(:,ph))),'x'); hold on
    polarplot(angles(:,ph),-ones(size(angles(:,ph))),'x'); hold off
    pause(0.5)
end

%% 2. Read in data + set up trajectory



for ph = phases
    [kdataframe, kinfo] = get_K_data_standardRadial( MeasFileID, 'Spokes', spokes(:,ph), 'Phases', ph);
    kdata(:,ph,:,:) = reshape(permute(kdataframe,[1,3,7,6,2,4,5]),[],1, size(encMat,1),32);
end

% Construct the trajectory
for ph = phases
    kspaceframe = gen_radial_traj(angles(:,ph),kinfo.hdr.Config.RawCol, []);
    kspaceframe = [-kspaceframe(:,1) kspaceframe(:,2)];
    kspace(:,ph,:,:) = repmat(reshape(kspaceframe, [],1,1,2), 1, 1, size(encMat,1));
end

% Reconstruction matrix size
xdims = kinfo.hdr.Config.NImageCols;
ydims = kinfo.hdr.Config.NImageLins;
zdims = kinfo.hdr.Config.NPar;

%% 3. coil sensitivity estimates
% Estimate sensitivity maps
if exist(['coils/' kinfo.hdr.Dicom.tProtocolName '_coils_R' num2str(acc_factor) '.mat'], 'file')
    load(['coils/' kinfo.hdr.Dicom.tProtocolName '_coils_R' num2str(acc_factor) '.mat'], 'sens')
else

    E0 = xfm_NUFFT_VEASL([xdims ydims zdims phases(end) 1, 1],ones([xdims ydims zdims]),[],kspace(:,:,1,:),1);

    ims = zeros([xdims, ydims, zdims, size(kdata,4)]);

    for c = 1:size(kdata,4)
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(size(kdata,4))])
        ims(:,:,:,c) =   mean( E0'*( mean(kdata(:,:,:,c),3) .* E0.w .* E0.norm ) , 4 );
    end

    disp('Adaptive combine estimation of coil sensitivities')
    sens = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]), 'kernel', 10, 'thresh', 0.05);

    clear('ims','E0')
    save(['coils/' kinfo.hdr.Dicom.tProtocolName '_coils_R' num2str(acc_factor) '.mat'], 'sens')
end
%% 5. Generate operator

disp('Generating operator')
E = xfm_NUFFT_LinearEncoding([xdims ydims zdims phases(end), 4], sens, [], kspace, encMat);

%% 3. phase correction + pre-weighting

disp('Phase correction')

% % Find the phase relative to the first repeat/average
% for ii = 1:size(kdata,2) % Loop through phases
%     for jj = 1:size(kdata,3) % Loop through Repeats/Averages (measurements)
%         % Complex number times it's complex conjugate = sqaure of magnitude of complex number
%         p = angle( mean( reshape( kdata(:,ii,jj,:) .* conj(kdata(:,ii,1,:)) , [] , 1 ) ) ); % Mean phase difference of each kspace to that of the first measurement (mean over kspace dims and channels)
%         kdata(:,ii,jj,:) = kdata(:,ii,jj,:).*exp(-1j*p); % Subtract this phase difference from the data
%     end
% end
% weight = repmat(E.w, 1, 1, 1, E.Nc);
% kdata = weight .* kdata;

Siz = size(kdata);
kdata_pc = reshape(kdata, [xdims+ydims, numSpokes/R, Siz(2:4)]);

% Find the phase relative to the first repeat/average
for ii = 1:size(kdata_pc,2) % Loop through radial spokes
    for jj = 1:size(kdata_pc,3) % Loop through phases
        for kk = 1:size(kdata_pc,4) % Loop through Repeats/Averages (measurements)
            % Complex number times it's complex conjugate = sqaure of magnitude of complex number
            p = angle( mean( reshape( kdata_pc(:,ii,jj,kk,:) .* conj(kdata_pc(:,ii,jj,1,:)) , [] , 1 ) ) ); % Mean phase difference of each kspace to that of the first measurement (mean over kspace dims and channels)
            kdata_pc(:,ii,jj,kk,:) = kdata_pc(:,ii,jj,kk,:).*exp(-1j*p); % Subtract this phase difference from the data
        end
    end
end

weight = repmat(E.w, 1, 1, 1, E.Nc);
kdata_pc = weight .* reshape(kdata_pc,Siz);
 
%% 5. Reconstruct

disp('Reconstruct data')

% setup l2 regularisers
lambdas = [];
l2_operators = {};
function_inputs = {};

lambdas = cat(2,lambdas, 0.1);
l2_operators = cat(2,l2_operators, {@temporal_finite_difference_tenc});
function_inputs = cat(2,function_inputs, {{}});


imFinal = fista_general(kdata_pc, E, 1, 0.000001, 100, 0.01, ...
    'showProgress', 1, ...
    'lambdas', lambdas, ...
    'l2_operators', l2_operators, ...
    'function_inputs', function_inputs);


%% 6. Reorder TE images

if TE
    imFinal = reshape(imFinal(:,:,:,:,end:-1:1), size(imFinal,1), size(imFinal,2), size(imFinal,3), []);
end

%% Save

% Get resolution information for NIFTI file
baseRes = kinfo.hdr.Meas.ReadFoV / kinfo.hdr.Meas.BaseResolution;
sliceThick = kinfo.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
TR = kinfo.hdr.MeasYaps.alTR{1} / 1e6; % in seconds

% Save in TWIX folder 
name = ['Data_tempsmoothall/' kinfo.hdr.Dicom.tProtocolName '_offlineRecon_CS_R' num2str(acc_factor) '_nobg.mat'];
save(name, 'imFinal')
% % save_avw(rot90(abs(imFinal),-1),name,'d',[baseRes,baseRes,sliceThick,TR])
% save_avw(rot90(abs(imFinal),-1),name,'d',[baseRes,baseRes,sliceThick,TR])
% 
% % Copy geometry information from NIFTI file
% dataToCopyGeom = dir([folder 'NIFTI/*' kinfo.hdr.Dicom.tProtocolName '.nii.gz']);
% system(['fslcpgeom ' dataToCopyGeom.folder '/' dataToCopyGeom.name ' ' name ' -d'])



