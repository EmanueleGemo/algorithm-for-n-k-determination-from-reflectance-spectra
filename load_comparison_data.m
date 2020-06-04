function comparison = load_comparison_data(comparison,range,kparam)
% loads the data of the optical properties that the output solution 
%would be compared with
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
    % adjust parameter for low k
    if nargin<3
        kparam = 0.001;
    end

    % check for adjust parameter
    if isfield(comparison,'adjust_parameter')
        adjust_parameter = comparison.adjust_parameter;
    else
        adjust_parameter = 1;
    end
    
    % here we assume standard wavelength/n/k columns assignment for
    % each dataset.
    if ~isfield(comparison,'l_column')
        comparison.l_column = [1]; % wl Position
    end
    if ~isfield(comparison,'n_column')
        comparison.n_column = [2]; % n column Positions
    end
    if ~isfield(comparison,'k_column')
        comparison.k_column = [3]; % k column Positions
    end
    
    % load data from .mat files
    if isfield(comparison,'material') && ~isempty(comparison.material{1})
        comparison.ndata = zeros(numel(range),1);
        comparison.kdata = comparison.ndata;
        comparison.edata = comparison.ndata;
        % n,k data
        file = load([comparison.material{1},'.mat']);
        structfields = fieldnames(file);
        if file.(structfields{1})(comparison.l_column,1) > 1  % check if it's in meters or nanometers
            rr = range*1e9; % nm
        else
            rr = range; % m
        end
        comparison.ndata = pchip(...
            file.(structfields{1})(:,comparison.l_column),...
            adjust_parameter*file.(structfields{1})(:,comparison.n_column),...
            rr);
                
        file.(structfields{1})(file.(structfields{1})(:,comparison.k_column)<kparam,comparison.k_column) = kparam;
        comparison.kdata = pchip(...
            file.(structfields{1})(:,comparison.l_column),...
            adjust_parameter*file.(structfields{1})(:,comparison.k_column),...
            rr);
        comparison.kdata(comparison.kdata(:,1) < kparam,:) = kparam;
        comparison.edata = nk2e(comparison.ndata,comparison.kdata);
    end
    if isfield(comparison,'rdata') && ~isempty(comparsion.rdata)
        comparison.rdata = zeros(numel(range),comparison.N);
        % r data
        file = load([comparison.reflectance,'.mat']);
        structfields = fieldnames(file);
        if file.(structfields{1})(comparison.l_r_column(ii),1) > 1  % check if it's in meters or nanometers
            rr = range*1e9; % nm
        else
            rr = range; % m
        end
        comparison.rdata = pchip(...
            file.(structfields{1})(:,comparison.l_r_column),...
            adjust_parameter*file.(structfields{1})(:,comparison.r_column),...
            rr);
    end
    comparison.range = range;  
end