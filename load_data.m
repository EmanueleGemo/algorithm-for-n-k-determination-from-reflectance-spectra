function [range,structure,spectrum,reference,comparison] = load_data(...
    range,structure,spectrum,reference,comparison)
% this function loads the data form .mat files
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

    if isrow(range)
        range = range';
    end
    
    [spectrum,range] = load_reflectance_data(spectrum,range);
    %   material data (structure)
    structure = load_structure_data(structure,range);
    %   reference data (investigated material)
    reference = load_reference_data(reference,range);
    %   comparison data (IF PRESENT, it will also be plotted with the solutions)
    comparison = load_comparison_data(comparison,range);
end