function [] = recon_adj(datafilename, encMat, R, sensfile, datafolder, savefolder)
%RECON_ADJ Recoinstructs radially acquired TE-ASL data by an adjoint
%operation
%   This function reads in radially acquired time encoded ASL data from a
%   Siemens twix file (.dat). It then keeps only every Rth spoke, uniformly
%   distributed among the Look Locker readouts. 
%
%   As pre-processing, first the coil sensitivities are estimated from the
%   retrospectively undersampled data by creating an average image friom
%   all spokes, the phase was corrected on a line by line basis to correct
%   for scanner drift during the scan. 
%
%   The reconstruction was performed as a simple adjoint operation of the
%   forward linear operator E that contains the linear time encoding
%   matrix, the coil sensitivity profile and the non-uniform Fourier
%   transform.
% 
%   The reconstruction is saved with the same name as the original data
%   file + the reconstruction method and acceleration factor (R).
%
%   INPUTS:
%           datafilename        string with the file name of the .dat file
%           encMat              encoding matrix, encodings x timepoints
%           R                   "acceleration factor" from the total number
%                               of acquired spokes (default 1)
%           sensfile            .mat file containing sensitivity estimate.
%                               If not provided, the sensitivity profile
%                               will be estimated from the data, and saved
%                               in a file with the provided name.
%           datafolder          If the raw data is not in the current
%                               directory the path to the file can be 
%                               provided. (optional)
%           savefolder          Where to save the output. Defaults to the
%                               current directory. (optional)
%
%   DEPENDENCIES:
%           The Michigan Image Reconstruction Toolbox for the NUFFT
%           (https://web.eecs.umich.edu/~fessler/code/)
%
%           MapVBVD to read twix data
%           (https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD)
%
%
%   Sophie Schauman 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. setup + sanity checks
% Check mandatory input and set default values to others
if nargin < 1
    error('no data file entered')
elseif nargin < 2
    error('No encoding matrix provided')
elseif nargin < 3
    R = 1;
    warning('No R provided. R is set to 1')
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
    spokes = floor((ph-1)/length(phases)*R)+1:R:numSpokes;
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
    kspace(:,ph,:,:) = repmat(reshape(kspaceframe, [],1,1,2), 1, 1, numEncodings);
    
end


%% 2. coil sensitivity estimates
% Estimate sensitivity maps
if exist(sensfile,'file')
    load(sensfile, 'sens')
else
    E0 = xfm_NUFFT_LinearEncoding([xdims ydims zdims phases(end) 1, 1],ones([xdims ydims zdims]),[],kspace(:,:,1,:),1);

    ims = zeros([xdims, ydims, zdims, size(kdata,4)]);

    for c = 1:size(kdata,4)
        disp(['Calculating coil sensitivity map ' num2str(c) ' of ' num2str(size(kdata,4))])
        ims(:,:,:,c) =   mean( E0'*( mean(kdata(:,:,:,c),3) .* E0.w .* E0.norm ) , 4 );
    end

    disp('Adaptive combine estimation of coil sensitivities')
    sens = adaptive_estimate_sens('data', permute(ims,[4,1,2,3]), 'kernel', 10, 'thresh', 0.05);

    clear('ims','E0')
    save(sensfile, 'sens')
end
%% 3. Generate operator
disp('Generating operator')
E = xfm_NUFFT_LinearEncoding([xdims ydims zdims phases(end), 4], sens, [], kspace, encMat);

%% 3. phase correction + pre-weighting
disp('Phase correction')
Siz = size(kdata);
kdata_pc = reshape(kdata, [xdims+ydims, numSpokes/R, Siz(2:4)]);

% Find the phase relative to the first repeat/average
for ii = 1:size(kdata_pc,2) % Loop through radial spokes
    for jj = 1:size(kdata_pc,3) % Loop through phases
        for kk = 1:size(kdata_pc,4) % Loop through encodings
            % Complex number times it's complex conjugate = sqaure of magnitude of complex number
            p = angle( mean( reshape( kdata_pc(:,ii,jj,kk,:) .* conj(kdata_pc(:,ii,jj,1,:)) , [] , 1 ) ) ); % Mean phase difference of each kspace to that of the first measurement (mean over kspace dims and channels)
            kdata_pc(:,ii,jj,kk,:) = kdata_pc(:,ii,jj,kk,:).*exp(-1j*p); % Subtract this phase difference from the data
        end
    end
end

weight = repmat(E.w, 1, 1, 1, E.Nc);
kdata_pc = weight .* reshape(kdata_pc,Siz);
 
%% 4. Reconstruct
disp('Reconstruct data')
imFinal = E'*kdata_pc;

%% 5. Reorder TE images
imFinal = reshape(imFinal(:,:,:,:,end:-1:1), size(imFinal,1), size(imFinal,2), size(imFinal,3), []);

%% 6. Save
name = [savefolder datafilename(1:end-4) '_offlineRecon_adj_R' num2str(R) '_nobg.mat'];
save(name, 'imFinal')
disp('Done')
end

