function [spectrum,range] = load_reflectance_data(spectrum,range)
% load spectrum file, finds lambda res (minimum of R), and integrate
% range by such datapoint.
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

    % constants 
    res = 1e6;
    
    if isrow(range)
        range = range';
    end

    % here we assume standard wavelength/n/k columns assignment for
    % each dataset.
    if ~isfield(spectrum,'l_column')
        spectrum.l_column = [1]; % wl Position
    end
    if ~isfield(spectrum,'r_column')
        spectrum.r_column = [2]; % n column Positions
    end
    
    
    % load file
	file = load([spectrum.file,'.mat']);
    
    structfields = fieldnames(file);  
    if file.(structfields{1})(1,spectrum.l_column) > 1
        file.(structfields{1})(:,spectrum.l_column) = ...
            file.(structfields{1})(:,spectrum.l_column)*1e-9;
    end
    
    l = file.(structfields{1})(:,spectrum.l_column);
    r = file.(structfields{1})(:,spectrum.r_column);
    
    % find r min (and lambda_res)
    li = (l(1):(l(end)-l(1))/(res-1):l(end))';
    ri = pchip(l,r,li);
    [~,idx] = min(ri); 
    idx = idx(1); %dumb fix of possible issues
    lambda_res = li(idx);
    
    % adjust range and data to fit lowest value of ri (and lambda_res)
    if all(idx~=range)
        pos = find(range>lambda_res,1,'first');
        range = [range(1:pos-1);lambda_res;range(pos:end)];
    end
    
    spectrum.range = range;
    spectrum.N = length(spectrum.range); % CAREFUL - spectrum.N is different from the .N field of structure and reference
    
    spectrum.rdata = pchip(file.(structfields{1})(:,spectrum.l_column), ...
        file.(structfields{1})(:,spectrum.r_column),range);
    if isrow(spectrum.rdata)
        spectrum.rdata = spectrum.rdata';
    end
    
end