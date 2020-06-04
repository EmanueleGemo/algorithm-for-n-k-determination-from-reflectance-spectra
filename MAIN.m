% example script
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

close all;
clear;
clc;

%% variables
% investigation range 
range = (800e-9 : 5e-9 : 1600e-9)'; % wavelength [m];
% structure data
structure.material = {'SiO2',[],'Al_v','Si'}; % material filenames; investigated material - EMPTY
structure.thickness = [8e-9,177e-9,100e-9,525e-5]; % thickness [m] 
% experimental data
spectrum.file = 'ref_cGeTe_205_10min300deg'; % measurement filename (.mat)
% reference data
reference.material = {'cGeTe_Jafari','cGeTe_Shportko','cGeTe_Park','cGeTe_Peinado_3','cGeTe_Robertson'}; % refs filenames
% comparsion data
comparison.material = {'cGeTe_liam'}; % reference or ellipsometry data, to compare with calculation output

%% calculate
% load data
[range,structure,spectrum,reference,comparison] = load_data(range,structure,spectrum,reference,comparison);

% build artificial reference
[artificial_reference,reference] = build_artificial_reference(...
    range,structure,spectrum,reference,3);

% calculate solution
solution = cavity_patternsearch_v1(structure,spectrum,reference,comparison,1,true);

disp('done')