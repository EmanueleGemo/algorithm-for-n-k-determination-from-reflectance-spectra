function structure = load_structure_data(structure,range,kparam)
% loads the structure data
%
% Code author: 
% Emanuele Gemo - University of Exeter - 04/06/2020
% 
% cite as:
% E. Gemo, S.V. Kesava, C. Ruíz De Galarreta, L. Trimby, S. García-Cuevas
% Carrillo, M. Riede, A. Baldycheva, A.M. Alexeev, and C.D. Wright, " A
% simple technique for determination of the refractive index of
% phase-change materials using near-infrared reflectometry", Optical
% Material Express, under peer review (2020)

    % adjust param for low k
    if nargin<3
        kparam= 0.001;
    end
    
    structure.investigated = find(cellfun(@isempty,structure.material))';
    if isempty(structure.investigated)
        error('The material stack does not contain empty slots assigned to the investigated material. Check the structure.material cell variable');
    end
    
    % here we assume standard wavelength/n/k columns assignment for
    % each dataset.
    if ~isfield(structure,'l_column')
        structure.l_column = [1]; % wl Position
    end
    if ~isfield(structure,'n_column')
        structure.n_column = [2]; % n column Positions
    end
    if ~isfield(structure,'k_column')
        structure.k_column = [3]; % k column Positions
    end

    % find size of structure
    ll = [length(structure.material);...
        length(structure.l_column);...
        length(structure.n_column);...
        length(structure.k_column);...
        length(structure.thickness)];
    structure.N = ll(1);
    if ll(2) == 1
        structure.l_column = ones(1,structure.N)*structure.l_column;
        ll(2) = structure.N;
    end
    if ll(3) == 1
        structure.n_column = ones(1,structure.N)*structure.n_column;
        ll(3) = structure.N;
    end
    if ll(4) == 1
        structure.k_column = ones(1,structure.N)*structure.k_column;
        ll(4) = structure.N;
    end
    if numel(unique(ll))>1
        error('STRUCTURE variables not consistent - please check');
    end
    
    % load variables from .mat files
    structure.ndata = zeros(numel(range),structure.N);
    structure.kdata = zeros(numel(range),structure.N);
    for ii = 1:structure.N
        if ii~=structure.investigated
            file = load([structure.material{ii},'.mat']);
            structfields = fieldnames(file);
            if file.(structfields{1})(structure.l_column(ii),1) > 1 % check if it's in meters or nanometers
                rr = range*1e9; % nm
            else
                rr = range; % m
            end
            structure.ndata(:,ii) = pchip(file.(structfields{1})(:,structure.l_column(ii)),...
                file.(structfields{1})(:,structure.n_column(ii)),rr);
            file.(structfields{1})(file.(structfields{1})(:,structure.k_column(ii))<kparam,structure.k_column(ii)) = kparam;
            structure.kdata(:,ii) = pchip(file.(structfields{1})(:,structure.l_column(ii)),...
                file.(structfields{1})(:,structure.k_column(ii)),rr);
            
        else
            % here lies investigated material - stays 0 for now.
        end
    end
    structure.edata = nk2e(structure.ndata,structure.kdata);
    structure.range = range;
end