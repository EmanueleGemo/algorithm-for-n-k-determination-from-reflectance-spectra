function Rx = TMM_fun_reduced(lambda,thickness,epsilon_relative)
% Reflectance calculation Transfer Matrix Method, E. Gemo, UoExeter
% Model adapted from R. C. Rumpf, DOI: pierb35/13.11083107
% Implementation adapted from S.Rao, mathworks - fileexchange - 47637
%
% The function calculates the reflection of a layer stack.
% highly optimized version - we here include the following approximations:
% non - magnetic materials only (mu_relative for each layer = 1+0i);
% non - polarized light hitting at normal incidence (theta, phi = 0) 
% environment material: air (epsilon_relative = 1+0i)
%
% inputs:
% lambda - vector(Lx1) containing the values of the incoming light
% wavelengths (# of values = L), in meters
% thickness - vector(1xT) containing the thickness of each layer composing
% the reflective device investigated (# of values = T), in meters
% epsilon_r - vector(LxT), containing the complex value of each layer
% dielectric function; each column T contains a set of values for the same
% layer pointed out in the thickness vector at the position T; each row
% refers to the same incoming beam wavelength, as specified in the lambda
% vector, at the respective position L.
%
% output:
% Rx - vector(Lx1) containing the reflectance of the layer stack at the
% wavelengths specified in lambda vector, for each respective position L.

    % propagation vector
    k0=(2*pi)./lambda;
    % eigenmode vector
    Vh = [0,-1i; 1i,0];
    % Rx init
    Rx = zeros(numel(k0),1);
    % Kz
    Kz = sqrt(epsilon_relative);
    
    % loop over ER
    for LL = 1:length(k0)  
        
        % initialize scattering matrices
        Sg11 = zeros(2);
        Sg12 = eye(2);
        Sg21 = eye(2);
        Sg22 = zeros(2);

        % sequentially updating the scattering matrices for each layer,
        for TT = 1:numel(thickness)
            
            % build Sd and So matrices (So = S21 = S12, Sd = S22 = S11)
            Q = [0, epsilon_relative(LL,TT);...
                -epsilon_relative(LL,TT), 0];
            Om = 1i*Kz(LL,TT)*eye(2);
            V = Q/Om;
            A = eye(2)+V\Vh;
            B = eye(2)-V\Vh;
            % X = expm(-Om*k0(jj)*thickness(layer)); below we use the 
            % direct expression as found in expm.m
%             X = diag(exp(full(diag(-Om*k0(LL)*thickness(TT)))));
            X = diag(exp(full(diag(-Om*k0(LL)*thickness(TT)))));
            Sd = (A-((X*B/A*X*B)))\((X*B/A*X*A)-B);
            So = ((A-((X*(B/A)*X*B)))\X)*(A-B/A*B);

            % updating scattering matrices
            Sg11=Sg11+(Sg12/(eye(2)-(Sd*Sg22))*Sd*Sg21);
            Sg12=Sg12/(eye(2)-(Sd*Sg22))*So;
            Sg21=So/(eye(2)-(Sg22*Sd))*Sg21;
            Sg22=Sd+(So/(eye(2)-(Sg22*Sd))*Sg22*So);
        end
        % retrieving Reflectance
        Rx(LL) = abs(Sg11(1)^(2));     
    end
end