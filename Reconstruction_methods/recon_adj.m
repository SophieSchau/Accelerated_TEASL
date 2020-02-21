function [] = recon_adj(datafilename, R, encMat, datafolder, savefolder)
%RECON_ADJ Recoinstructs radially acquired TE-ASL data by an adjoint operation
%   explanation
% 
%
%
%
%

%% 0. setup + sanity checks

% Check mandatory input and set default values to others
if nargin < 1
    error('no data file entered')
elseif nargin < 2
    R = 1;
    warning('No R provided. R is set to 1')
elseif nargin < 3
    encMat = 1;
    warning('No encoding matrix provided. encMat set to 1')
end

if ~exist('datafolder', 'var')
    datafolder = './';
end

if ~exist('savefolder', 'var')
    savefolder = './';
end

% check formatting of provided strings
if datafolder(end)~='/'
    datafolder = [datafolder '/'];
end

if savefolder(end)~='/'
    savefolder = [savefolder '/'];
end


%% 1. Read in data + set up trajectory

twix_obj = mapVBVD([datafolder, datafilename],'ignoreSeg');
kdata_all =  twix_obj.image;

numSpokes = twix_obj.hdr.Config.NLin;
phases = 1:twix_obj.hdr.Config.NPhs;
numChannels = twix_obj.hdr.Config.NChaMeas;
numEncodings = twix_obj.hdr.Config.NAve;
xdims = twix_obj.hdr.Config.NImageCols;
ydims = twix_obj.hdr.Config.NImageLins;
zdims = twix_obj.hdr.Config.NPar;

for ph = phases
    spokes = ceil((ph-1)/length(phases)*R+1:R:numSpokes);
    kdataframe = kdata_all(:,:,spokes,:,:,:,ph,:,:,:,:,:,:,:,:,:);
    kdata(:,ph,:,:) = reshape(permute(kdataframe,[1,3,7,6,2,4,5]),[],1, numEncodings,numChannels);

    angles = pi * (spokes-1) / (numSpokes);
    
    if mod(numSpokes,2) == 1
        angles = angles * 2;
    end
    % figure; hold on
    polarplot(angles(:),ones(size(angles(:))),'x'); hold on
    polarplot(angles(:),-ones(size(angles(:))),'x'); hold off
    pause(0.5)
    
    % Construct the trajectory
    kspaceframe = gen_radial_traj(angles(:),twix_obj.hdr.Config.RawCol, []);
    kspaceframe = [-kspaceframe(:,1) kspaceframe(:,2)];
    kspace(:,ph,:,:) = repmat(reshape(kspaceframe, [],1,1,2), 1, 1, size(encMat,1));
    
end


%% 2. coil sensitivity estimates
% Estimate sensitivity maps

E0 = xfm_NUFFT_VEASL([xdims ydims zdims phases(end) 1, 1],ones([xdims ydims zdims]),[],kspace(:,:,1,:),1);

ims = zeros([xdims, ydims, zdims, size(kdata,4)]);

for c = 1:size(kdata,4)
    disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(size(kdata,4))])
    ims(:,:,:,c) =   mean( E0'*( mean(kdata(:,:,:,c),3) .* E0.w .* E0.norm ) , 4 );
end

disp('Adaptive combine estimation of coil sensitivities')
sens = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]), 'kernel', 10, 'thresh', 0.05);

clear('ims','E0')
%% 5. Generate operator

disp('Generating operator')
E = xfm_NUFFT_VEASL([xdims ydims zdims phases(end), 4], sens, [], kspace, encMat);

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

imFinal = E'*kdata_pc;

%% 6. Reorder TE images

if TE
    Siz = size(imFinal);
    imFinal = reshape(imFinal, xdims, ydims, zdims, prod(Siz(4:end)));
    
    imFinal = imFinal(:,:,:,end:-1:1);
    
    for ii = 1 : size(encMat,2)
        indNow  = (1 :  1 : phases(end)) + (ii-1)*phases(end);
        indThen = (phases(end) : -1 : 1) + (ii-1)*phases(end);
        imFinal(:,:,:,indNow) = imFinal(:,:,:,indThen);
    end
end

%% Save

% Get resolution information for NIFTI file
baseRes = kinfo.hdr.Meas.ReadFoV / kinfo.hdr.Meas.BaseResolution;
sliceThick = kinfo.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
TR = kinfo.hdr.MeasYaps.alTR{1} / 1e6; % in seconds

% Save in TWIX folder 
name = ['Data/' kinfo.hdr.Dicom.tProtocolName '_offlineRecon_adj_R' num2str(R) '_nobg.mat'];
save(name, 'imFinal')
% % save_avw(rot90(abs(imFinal),-1),name,'d',[baseRes,baseRes,sliceThick,TR])
% save_avw(rot90(abs(imFinal),-1),name,'d',[baseRes,baseRes,sliceThick,TR])
% 
% % Copy geometry information from NIFTI file
% dataToCopyGeom = dir([folder 'NIFTI/*' kinfo.hdr.Dicom.tProtocolName '.nii.gz']);
% system(['fslcpgeom ' dataToCopyGeom.folder '/' dataToCopyGeom.name ' ' name ' -d'])







end

