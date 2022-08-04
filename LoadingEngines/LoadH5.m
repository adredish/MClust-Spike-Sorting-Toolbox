function [t, wv] = LoadH5(varargin)
% MClust loading engine for simple h5 files with spikes data.
% 
% Aug 2022
% Fernando Chaure

if nargin == 2
    % New "get" construction"
    if strcmp(varargin{1}, 'get')
        switch (varargin{2})
            case 'ExpectedExtension'
                t = '.h5'; return;
            otherwise
                error('Unknown get condition.');
        end
    else
        error('2 argins requires "get" as the first argument.');
    end
end

info = h5info(varargin{1});

if ~any(strcmp({info.Datasets.Name},'t')) || ~any(strcmp({info.Datasets.Name},'wv'))
    error('Invalid h5 file. The file should contain at least t and wv.')
end

if nargin == 3 && varargin{3}==5 %with the metadata is enough
    t = info.Datasets(strcmp({info.Datasets.Name},'t')).Dataspace.Size;
    return
end

wv = double(h5read(varargin{1},'/wv'));

t = double(h5read(varargin{1},'/t'));
t = t(:); %to column vector

if nargin == 1
    return
else % other nargin == 3
    list = varargin{2};
    switch (varargin{3})
        case 1 % timestamp list [they must be in the data!]
            sel = zeros(size(list));
            for i=1:numel(list)
                sel(i)= find(t==list(i),1,'first');
            end            
        case 2 % record number list
            sel = list;
        case 3 % range of timestamps
            assert(numel(list)==2, 'second argument must have 2 elements') 
            sel = t>=list(1) && t<=list(2);
        case 4 % range of records
            assert(numel(list)==2, 'second argument must have 2 elements') 
            sel = list(1):list(2);
        otherwise
            error('Unknown unit expectation.');
    end
end
    
t = t(sel);
wv = wv(sel,:,:);
