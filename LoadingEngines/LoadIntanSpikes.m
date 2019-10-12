function [varargout] = LoadIntanSpikes(varargin)
% MClust loading engine for Intan spike data.
% 
% INPUT
%   fname - filename of spikes file
%   records_to_get - (optional) which spikes to retreive
%   record_units - (optional) code for the format of records_to_get
%       1 - records_to_get is a timestamp list
%       2 - records_to_get is a record number list
%       3 - records_to_get is range of timestamps (a vector with 2 
%               elements: a start and an end timestamp)
%       4 - records_to_get is a range of records (a vector with 2
%               elements: a start and an end record number)
%       5 - return only the count of spikes (records_to_get should be [])
%       default - return all spikes
%
% OUTPUT
%   T - vector of spike times
%   WV - spike waveforms (#spikes x #channels x #samplesPerSpike matrix)
%
% Jan 2016
% Updated Apr 2017 w/ new spike file format (doubles instead of int16s)
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

%% ADDED code to fill keys
% ADR 2017-10-02
if nargin == 1
    fname = varargin{1};
elseif nargin == 2
    if strcmp(varargin{1}, 'get')
        switch (varargin{2})
            case 'ExpectedExtension'
                varargout{1} = '.spikes'; return;
            case 'AverageWaveform_ylim'
                varargout{1} = [-1000 1000]; return                
            otherwise
                error('Unknown get condition.');
        end
    else
        error('2 argins requires "get" as the first argument.');
    end
elseif nargin == 3
    fname = varargin{1};
    varargin = varargin(2:end);
end



%%


% Open file and read file header
fid = fopen(fname, 'r');
num_spikes = fread(fid, 1, 'uint32');   %first 4 bytes are # spikes
num_channels = fread(fid, 1, 'uint16'); %next 2 bytes are #channels
wind_len = fread(fid, 1, 'uint16');     %next 2 bytes are #samples/spike

% Figure out format waveform values (int16 or double?)
sz = dir(fname);
if sz.bytes == 8+8*num_spikes+2*num_spikes*num_channels*wind_len
    FORMAT = 0; %old format
    disp('LoadIntanSpikes: Using old int16 format')
else %otherwise assume new format which stores waveforms w/ doubles
    FORMAT = 1;
end

% Do we just want spike count?
if length(varargin)==2 && varargin{2}==5 %(5=code for just spike count)
    varargout{1} = num_spikes;  %then just return spike count
    varargout{2} = [];          %and don't worry about waveforms
    return
end

% Read the spiketimes
T = fread(fid, num_spikes, 'double');

% Read the spike waveforms
S = zeros(num_spikes, num_channels, wind_len); %to store spike waveforms
if FORMAT==0 %old format w/ int16
    Sna = fread(fid, num_spikes*num_channels*wind_len, 'int16');%waveforms
else
    Sna = fread(fid, num_spikes*num_channels*wind_len, 'double');%waveforms
end
fclose(fid); %close file

% Reorder data in memory
for s = 1:num_spikes
    for c = 1:num_channels
        ind = (s-1)*num_channels*wind_len+(c-1)*wind_len+1;
        S(s, c, :) = double(Sna(ind:(ind+wind_len-1)));
    end
end
S = -S; %flip so spikes are positive

% Set return values according to what was requested
if length(varargin)==2
    rtg = varargin{1};
    record_units = varargin{2};
    switch record_units
        case 1 %timestamp list
            rinds = arrayfun(@(x)find(T==x,1),rtg); %requested spikes
            varargout{1} = T(rinds);
            varargout{2} = S(rinds,:,:);
        case 2 %record number list
            varargout{1} = T(rtg);
            varargout{2} = S(rtg,:,:);
        case 3 %range of timestamps
            rinds = T>rtg(1) & T<rtg(2); %requested spikes
            varargout{1} = T(rinds);
            varargout{2} = S(rinds,:,:);
        case 4 %range of record numbers
            varargout{1} = T(rtg(1):rtg(2));
            varargout{2} = S(rtg(1):rtg(2),:,:);
        otherwise %return all spikes by default
            varargout{1} = T;
            varargout{2} = S;
    end
else % Default: return everything
    varargout{1} = T;
    varargout{2} = S;
end

end