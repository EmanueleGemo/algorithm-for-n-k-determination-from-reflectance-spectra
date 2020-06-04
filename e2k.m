function n = e2k(e)
% calculates n from epsilon_r
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

    rr = real(e);
    ii = imag(e);
    n = (1/2.*((rr.^2+ii.^2).^(1/2) - rr)).^(1/2);
end