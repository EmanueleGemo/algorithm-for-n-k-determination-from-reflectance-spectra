function str = optimal_thickness(range,reference) 
% yields the approximated optimal thickness for a given material
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

    m = 0;
    for ii = 1:reference.N
        m = m + mean(1/4*range'./sqrt((reference.ndata(:,ii)+ 1i*reference.kdata(:,ii)).*(reference.ndata(:,ii) - 1i*reference.kdata(:,ii))));
    end
    str = sprintf('optimal cavity thickness = %0.00f nm\n',m/reference.N*1e9);
end