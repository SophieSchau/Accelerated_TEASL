function [ kdata, kinfo, R ] = get_K_data( MeasFileID, varargin)
%GET_K_DATA Get suitably undersampled kdata, from 2D radial golden angle 
%vessel encoded ASL images from Siemens raw .dat file 
%
%   Sophie Schauman July 2018 - sophie.schauman@dtc.ox.ac.uk
%   Based partially on code by Thomas W Okell.
%
%   INPUTS: 
%   MeasFileID - .dat file identifier. Can be filename with path or just
%                   the MID from the .dat filename if it is on the defined path.
%   (R            - acceleration factor)
%   (PrepNos      - list of specified prep numbers to use)
%   (Frames       - indecies of frames to read in)
%
%   OUTPUTS:
%   kdata         - 4D matrix (nKPoints x nTPoints x  nEncodings x nCoils)
%   kinfo         - twix object with parsed information from .dat file
%   R             - R used (might have changed from input if input doesn't
%                   work with data.
%
%   DEPENDENCIES:
%   mapVBVD by Philipp Ehses (philipp.ehses@tuebingen.mpg.de) for reading 
%   .dat files. Available on:
%   https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
    % Read in extra inputs for reading data 
    p   =   inputParser;
    p.addParameter('Spokes',       []);
    p.addParameter('Phases', []);
    p.parse(varargin{:});
    p   =   p.Results;
    
    Spokes = p.Spokes;
    Phases = p.Phases;
    
    % Read the raw data headers
    twix_obj = mapVBVD(MeasFileID,'ignoreSeg');
    kinfo = twix_obj;
kdata =  twix_obj.image(:,:,Spokes,:,:,:,Phases,:,:,:,:,:,:,:,:,:);

        % Order of raw data:
        %  1) Columns
        %  2) Channels/Coils
        %  3) Lines
        %  4) Partitions
        %  5) Slices
        %  6) Averages
        %  7) (Cardiac-) Phases
        %  8) Contrasts/Echoes
        %  9) Measurements
        % 10) Sets
        % 11) Segments
        % 12) Ida
        % 13) Idb
        % 14) Idc
        % 15) Idd
        % 16) Ide
        
end

