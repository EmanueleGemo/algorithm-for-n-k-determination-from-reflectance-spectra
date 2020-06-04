function delta = search_fun_reflectance(p,cavity,idx,investigated,reference,spectrum)
% objective function to minimize
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
    
    cavity.ER(:,investigated) = reduced_weighted_average(p,reference.edata(idx,:));
    Rx = TMM_fun_reduced(reference.range(idx),cavity.L,cavity.ER);
    % approach: minimise difference between reflectivities
    delta = sum((spectrum.rdata(idx)-Rx).^2);
end