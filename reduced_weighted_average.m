function a = reduced_weighted_average(p,val)
% weighted average function
% it works on rows, NOT on columns
% can perform vectorised operations with multiple rows
% constraint: works with n-1 parameters, with sum of weights = 1.

% weighted average:
% V -> Values       length_of_V = n
% w -> weights      length_of_w = n
% weighted_average = sum(w.*V) / sum(w)

% reduced weighted average:
% V -> Values       length_of_V = n
% p -> weights      length_of_p = n-1
% w = [p,1-sum(p)]  sum_of_w = 1
% weighted_average = sum([p,1-sum(p)].*V)
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

    a = sum(val.*[p,1-sum(p,2)],2);
end