function reference = load_reference_data(reference,range,kparam)
% loads the reference data
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

    % adjust parameter for k if too low
    if nargin<3
        kparam = 0.001;
    end
    
    % here we assume standard wavelength/n/k columns assignment for
    % each dataset.
    if ~isfield(reference,'l_column')
        reference.l_column = [1]; % wl Position
    end
    if ~isfield(reference,'n_column')
        reference.n_column = [2]; % n column Positions
    end
    if ~isfield(reference,'k_column')
        reference.k_column = [3]; % k column Positions
    end
    
    % find size of reference
    ll = [length(reference.material);...
        length(reference.l_column);...
        length(reference.n_column);...
        length(reference.k_column)];
    reference.N = ll(1);
    if ll(2) == 1
        reference.l_column = ones(1,reference.N)*reference.l_column;
        ll(2) = reference.N;
    end
    if ll(3) == 1
        reference.n_column = ones(1,reference.N)*reference.n_column;
        ll(3) = reference.N;
    end
    if ll(4) == 1
        reference.k_column = ones(1,reference.N)*reference.k_column;
        ll(4) = reference.N;
    end
    if numel(unique(ll))>1
        error('REFERENCE variables not consistent - please check');
    end

    % check for adjust parameter
    if isfield(reference,'adjust_parameter')
        adjust_parameter = reference.adjust_parameter;
    else
        adjust_parameter = 1;
    end
    if numel(adjust_parameter)<reference.N
        adjust_parameter = adjust_parameter.*ones(1,reference.N);
    end
    
    % load data from .mat files
    reference.ndata = zeros(numel(range),reference.N);
    reference.kdata = reference.ndata;
    reference.edata = reference.ndata;
    for ii = 1: reference.N
        file = load([reference.material{ii},'.mat']);
        structfields = fieldnames(file);
        if file.(structfields{1})(reference.l_column(ii),1) > 1  % check if it's in meters or nanometers
            rr = range*1e9; % nm
        else
            rr = range; % m
        end
        reference.ndata(:,ii) = pchip(...
            file.(structfields{1})(:,reference.l_column(ii)),...
            adjust_parameter(ii)*file.(structfields{1})(:,reference.n_column(ii)),...
            rr);
        file.(structfields{1})(file.(structfields{1})(:,reference.k_column(ii))<kparam,reference.k_column(ii)) = kparam;
        reference.kdata(:,ii) = pchip(...
            file.(structfields{1})(:,reference.l_column(ii)),...
            adjust_parameter(ii)*file.(structfields{1})(:,reference.k_column(ii)),...
            rr);        
        reference.edata(:,ii) = nk2e(reference.ndata(:,ii),reference.kdata(:,ii));
    end
    reference.range = range;
end