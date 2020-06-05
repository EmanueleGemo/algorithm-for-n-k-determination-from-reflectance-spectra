function [artificial_reference,reference] = build_artificial_reference(...
    range,structure,spectrum,reference,fit_degree,force_e2_above_0,flag_plot)
% this function combines the provided references, structure and reflectance
% spectrum to infer a possible cavity material parameter dispersion which
% matches the reflectance.
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

    if nargin<7
        flag_plot = false;
    end
    if nargin<6
        force_e2_above_0 = true;
    end

    [e1,e2,dR] = calculate_dR_space(...
        structure,spectrum,reference,force_e2_above_0,flag_plot);
    [artificial_reference] = calculate_artificial_reference(...
    range,structure,spectrum,reference,e1,e2,dR,fit_degree,force_e2_above_0,flag_plot);

    reference.N = reference.N+1;
    reference.ndata = [reference.ndata,artificial_reference(:,2)];
    reference.kdata = [reference.kdata,artificial_reference(:,3)];
    reference.edata = [reference.edata,nk2e(artificial_reference(:,2),artificial_reference(:,3))];
    reference.material = [reference.material,{'artificial_reference'}];
end


function [e1,e2,dR] = calculate_dR_space(structure,spectrum,reference,force_e2_above_0,flag_plot)
% this function calculates dR and its e2/e1 coordinates, for the artificial
% reference building
% it performs an initial low res calculation to then anneal towards
% high-interest areas, and then interpolates to obtain a high-res dR space.
% to compensate for possible missing points, the lookup function anneals by
% areas, so it always looks for solutions on a set of quadrants in a raster

    % constants - do not change unless needed
    init_precision = 50;
    min_point_distance = 3;
    iter = 10;
    density = 50;
    quadrants = 5;
    depth_threshold = 7;
    fit_resolution = 1000;    
    expansion = 1.2;

    % calculate global search area by epslion min and max
    ar = floor(min(min(real(reference.edata(:,:))))); %left lim
    Ar = ceil(max(max(real(reference.edata(:,:))))); %right lim
    ai = floor(min(min(-imag(reference.edata(:,:))))); %down lim
    Ai = ceil(max(max(-imag(reference.edata(:,:))))); % up lim
    sar = floor(expansion*(Ar-ar)/2); % left-right corrected half-span
    car = round((Ar+ar)/2); % left-right span center
    sai = floor(expansion*(Ai-ai)/2); % u-d corrected half-span
    cai = round((Ai+ai)/2); % u-d span center
    car_sar = car-sar; % l-r corrected span left lim
    carsar = car+sar; % l-r corrected span right lim
    cai_sai = cai-sai; % u-d corrected span down lim
    caisai = cai+sai;% u-d corrected span up lim
    init_point_distance = 2*max([sar,sai])/init_precision; % initial span min resolution
    
    % here we force e2 to be bigger than 0
    if force_e2_above_0 && cai_sai<0
        cai_sai = 0;
    end

    % calculation coordinates
    [e1M,e2M] = meshgrid(...
        (car_sar:(carsar-car_sar)/(init_precision-1):carsar)',...
        (cai_sai:(caisai-cai_sai)/(init_precision-1):caisai)');
    e1M = e1M(:)+rand(size(e1M(:)));
    e2M = e2M(:)+rand(size(e1M));
    
    % interpolation coordinates    
    [e1,e2] = meshgrid(...
        (car_sar:(carsar-car_sar)/(fit_resolution-1):carsar)',...
        (cai_sai:(caisai-cai_sai)/(fit_resolution-1):caisai)');
    
    rdata = spectrum.rdata;
    
    WaitMessage = parfor_wait(spectrum.N, 'Waitbar', true);
    parfor idx = 1:spectrum.N
%     for idx = 1:spectrum.N
  
        max_point_distance = init_point_distance;
        
        % raster
        X = e1M;
        Y = e2M;
        
        % get cavity
        ER = structure.edata(idx*ones(init_precision^2,1),:);
        ER(:,structure.investigated) = X-1i*Y;
        
        R = TMM_fun_reduced(ones(init_precision^2,1)*structure.range(idx),structure.thickness,ER);

        Delta = ((rdata(idx)-R).^2);

%         if flag_plot
%             figure(1); 
%             hold off;
%             sz = 10;
%             scatter(X,Y,sz,log10(Delta),'s','filled'); 
%             xlim([car_sar carsar]);
%             ylim([cai_sai caisai]);
%             xlabel('Real({\epsilon})');
%             ylabel('Imag({\epsilon})');
%             title(sprintf('wl = %0.0f nm',structure.range(idx)*1e9));
%             colormap jet; 
%             hold on;
%             drawnow;        
%             ss = gobjects(quadrants^2,1);
%         end

        delta = cell(quadrants^2,1);
        log10Delta = delta;
        x = cell(quadrants^2,1);
        y = cell(quadrants^2,1);
        xx = min(X):(max(X)-min(X))/(quadrants):max(X)*1.01;
        yy = min(Y):(max(Y)-min(Y))/(quadrants):max(Y)*1.01;
        kk = 0;
        for ii = 1:quadrants  
            xmap = X>=xx(ii) & X<xx(ii+1);
            for jj = 1:quadrants
                ymap = Y>=yy(jj) & Y<yy(jj+1);
                map = xmap & ymap;
                kk = kk+1;
                delta{kk} = Delta(map);
                x{kk} = X(map);
                y{kk} = Y(map);
                
%                 if flag_plot
%                     ss(kk) = scatter(x{kk},y{kk},sz,ones(nnz(map),1),'ro');
%                     drawnow;
%                 end
                
            end
        end
        
%         if flag_plot
%         	delete(ss);
%         end
        
        jj = 0;
        while jj<iter
            jj = jj+1;
            
            % determine which areas to investigate
            
            max_point_distance = max_point_distance-(max_point_distance-min_point_distance)*jj/iter; 
            mm = 0; nn = 0;
            for kk = 1: quadrants^2
                log10Delta{kk} = -log10(delta{kk});
                log10Delta{kk}(log10Delta{kk}>depth_threshold) = depth_threshold;
                mm = mm+ sum(log10Delta{kk});
                nn = nn + numel(log10Delta{kk});
            end
            mm = mm/nn;
            for kk = 1: quadrants^2
                
                map = log10Delta{kk}>mm;%*(.5*jj/iter);
                newsites = nnz(map);
                if newsites < 1
                    continue
                end
                
%                 if flag_plot
%                     ss(kk) = scatter(x{kk}(map),y{kk}(map),sz,ones(newsites,1),'ro');
%                     drawnow;
%                 end

                % calculate new coords
                [x{kk},y{kk}] = randcoords_points([x{kk}(map),y{kk}(map)],density,max_point_distance);

                % calculate reflectivities
                n = length(x{kk});                
                ER = structure.edata(idx*ones(n,1),:);
                ER(:,structure.investigated) = x{kk}-1i*y{kk};
                r = TMM_fun_reduced(ones(n,1)*structure.range(idx),structure.thickness,ER);
                delta{kk} = ((rdata(idx)-r).^2);

                % integrate new variables
                R = [R;r];
                X = [X;x{kk}];
                Y = [Y;y{kk}];
                Delta = [Delta;delta{kk}]; 
            
%                 if flag_plot
%                     delete(ss)
%                     scatter(x{kk},y{kk},sz,log10(delta{kk}),'s','filled')
%                     drawnow;
%                 end
            end
        end

        
        % interpolates
        F = scatteredInterpolant(X,Y,Delta,'linear');
        zq = zeros(size(e1));
        zq(:) = F(e1(:),e2(:));
        dR(:,:,idx) = abs(zq);
        
%         if flag_plot
%             figure(1); 
%             hold off;
%             sz = 10;
%             
%             [xqp,yqp] = meshgrid(...
%                 (car_sar:(carsar-car_sar)/(999):carsar)',...
%                 (cai_sai:(caisai-cai_sai)/(999):caisai)');
%             zqp = zeros(size(xqp));
%             zqp(:) = F(xqp(:),yqp(:));
%             surf(xqp,yqp,zqp);    
%             view(0,90);
%             shading flat;
% 
%             xlim([car_sar carsar]);
%             ylim([cai_sai caisai]);
%             xlabel('Real({\epsilon})');
%             ylabel('Imag({\epsilon})');
%             title(sprintf('wl = %0.0f nm',structure.range(idx)*1e9));
%             colormap jet; 
%             hold on;
%             drawnow;
%         end
        WaitMessage.Send; %#ok<PFBNS>
    end
    WaitMessage.Destroy;

   
end
function [artificial_reference] = calculate_artificial_reference(...
    range,structure,spectrum,reference,e1,e2,delta_R,fit_degree,force_e2_above_0,flag_plot)
% end
% function [sol,nksol,nksol_centerofmass,nksol_fitted] = FUNCTION_FIT_cavity_prototype_4_postprocess(partial_solution_fitted,solution,spectrum,structure,reference,comparison,range)
    
    % constants - do not change unless needed
    circle_threshold = .5;
    ec_threshold = 0.85;
    
    %% find center of solution @ cavity resonant wavelength - rough method
    [~,lres_idx] = min(spectrum.rdata);
    
    normdR = normalize(-log10(delta_R(:,:,lres_idx))); % also - taking the logarithm to enhance differences
    
% %     normdR = normdR - min(min(normdR));
% %     normdR  = normdR./max(max(normdR));
    normdR(normdR<circle_threshold) = 0;
    normdR(normdR>=circle_threshold) = 1;

    if flag_plot
        surf(e1,e2,normdR);
        title('solution centroid at resonance')
        shading flat;
        view(0,90);
    end    
    
    stats = regionprops(normdR,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
    e_dR = [e1(1,round(stats.Centroid(1))),...
        e2(round(stats.Centroid(2)),1)];
    
    tmp = [round(stats.Centroid(1)-0.5*stats.MajorAxisLength),...
        round(stats.Centroid(1)+0.5*stats.MajorAxisLength)];
    delta_e1 = abs(diff(interp1(1:size(e1,2),e1(1,:),tmp,'linear','extrap')));
    tmp = [round(stats.Centroid(2)-0.5*stats.MinorAxisLength),...
        round(stats.Centroid(2)+0.5*stats.MinorAxisLength)];
    delta_e2 = abs(diff(interp1(1:size(e2,1),e2(:,1),tmp,'linear','extrap')));

    delta_e = (delta_e1+delta_e2)/2;

    %% calculating derivatives
    % calculation of references derivative
    de_dl_refs = cdiff(reference.edata,'E'); 
    % calculation of reference average
    refs_average = [mean(real(reference.edata),2),mean(-imag(reference.edata),2)];
    % calculation of references average derivative
    de_dl_average = cdiff(refs_average,'E'); 
    
    e_guess = zeros(spectrum.N,2);
    e_guess(lres_idx,:) = e_dR;
    
    
    % vector generation and application to e_guess
    t = zeros(spectrum.N,1); 
    
    % forward
    u = @(idx) de_dl_average(idx-1,:);
    v = @(idx) de_dl_average(idx,:);
    theta_0 = acos(dot(u(lres_idx),v(lres_idx))/(norm(u(lres_idx))*norm(v(lres_idx))));
    v0 = refs_average(lres_idx,:)-e_dR;
    SF = @(idx) sqrt(sum((refs_average(idx,:)+refs_average(idx-1)).^2))...
        /sqrt(sum((refs_average(lres_idx,:)+refs_average(lres_idx-1)).^2));
    theta = @(idx) acos(dot(u(idx),v(idx))/(norm(u(idx))*norm(v(idx)))) ...
        - theta_0;
    R = @(t) [cos(t), -sin(t);...
        sin(t), cos(t)];
    
    for idx = lres_idx+1:spectrum.N
        t(idx) = t(idx-1)+theta(idx);
        e_guess(idx,:) = refs_average(idx,:) - (R(t(idx))*v0'*SF(idx))';
    end 
    
    % backward
    u = @(idx) de_dl_average(idx+1,:);
    v = @(idx) de_dl_average(idx,:);
    theta_0 = acos(dot(u(lres_idx),v(lres_idx))/(norm(u(lres_idx))*norm(v(lres_idx))));
    SF = @(idx) sqrt(sum((refs_average(idx,:)+refs_average(idx+1)).^2))...
        /sqrt(sum((refs_average(lres_idx,:)+refs_average(lres_idx+1)).^2));
    theta = @(idx) acos(dot(u(idx),v(idx))/(norm(u(idx))*norm(v(idx)))) ...
        - theta_0;
    for idx = lres_idx-1:-1:1
        t(idx) = t(idx+1)+theta(idx);
        e_guess(idx,:) = refs_average(idx,:) - (R(t(idx))*v0'*SF(idx))';
    end
    
    % application of gaussian to probability plots, and determination of
    % probability space
    
    %% generation of G space and C space
    G = zeros(size(delta_R));
    C = G;
    e_c = ones(spectrum.N,1)+1i*ones(spectrum.N,1);
    
    for idx = 1:spectrum.N
        % build G and normalize
        G(:,:,idx)= normalize(...
            exp(-1/2*((e1-e_guess(idx,1))./delta_e).^2)...
            .*exp(-1/2*((e2-e_guess(idx,2))./delta_e).^2));
        % retrieve dR and normalize
        normdR = normalize(-log10(delta_R(:,:,idx)));
        % build C and normalize
        C(:,:,idx) = normalize(G(:,:,idx).*normdR);        
        % determining C center
        condition = true;
        tmp = ec_threshold;
        while condition
            A = C(:,:,idx)>tmp;        
            if nnz(A) == 0
                tmp = tmp-0.01;
            else
                condition = false;
            end
        end
        ec_idxs = find(A);
        e_c(idx) = (sum(e1(ec_idxs).*C(ec_idxs))...
            -1i.*sum(e2(ec_idxs).*C(ec_idxs)))./sum(C(ec_idxs));        
        
        if flag_plot
            surf(e1,e2,C(:,:,idx));
            title(sprintf('ec calculation @ %0.3g m',range(idx)))
            shading flat;
            view([0 90]);
            hold on;
            scatter3(real(e_c(idx)),-imag(e_c(idx)),30,'kx','LineWidth',1);
            hold off;
            drawnow
        end
    end
    
    % getting nk values

    n_c = e2n(e_c); %(1/2*((real(e_c).^2+(-imag(e_c)).^2).^(1/2)+real(e_c))).^(1/2);
    k_c = e2k(e_c); %(1/2*((real(e_c).^2+(-imag(e_c)).^2).^(1/2)-real(e_c))).^(1/2);
    
    % interpolating those by pinning values around lambda res
    n = n_c(lres_idx)+polyval(polyfit_0_intercept(1e9*(range-range(lres_idx)),...
        n_c-n_c(lres_idx),fit_degree),1e9*(range-range(lres_idx)));
    k = k_c(lres_idx)+polyval(polyfit_0_intercept(1e9*(range-range(lres_idx)),...
        k_c-k_c(lres_idx),fit_degree),1e9*(range-range(lres_idx)));
    
    if flag_plot
        cla;
        hold on;
        yyaxis left
        plot(range,n_c,'b+',range,n,'g-');
        ylabel('n');
        xlabel('wl [m]');
        yyaxis right
        plot(range,k_c,'r+',range,k,'g--');
        ylabel('k');        
    end
    
    
    artificial_reference = [range,n,k];

end


%% sub-sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = cdiff(a,cmd)
    if isempty(a)
        d = [];
    else
        n = numel(a);
        if n == 1
            d = a;
        else
            d = (a(3:end,:)-a(1:end-2,:))./2;
            if nargin>1
                switch cmd
                    case 'F' % forward diff for 1st
                        tmp = a;
                        tmp(1,:) = fdiff(a(1:2,:));
                        tmp(2:end-1,:) = d;
                        tmp = tmp(1:end-1,:);
                    case 'B' % backward diff for last
                        tmp = a;
                        tmp(end,:) = bdiff(a(end-1:end,:));
                        tmp(2:end-1,:) = d;
                        tmp = tmp(2:end,:);
                    case 'E' % both ends corrected
                        tmp = a;
                        tmp(1,:) = fdiff(a(1:2,:));
                        tmp(end,:) = bdiff(a(end-1:end,:));
                        tmp(2:end-1,:) = d;
                end
                d = tmp;
            end
        end
    end
end
function d = fdiff(a)
    if isempty(a)
        d = [];
    else
        n = numel(a);
        if n == 1
            d = a;
        else
            d = a(2:end,:)-a(1:end-1,:);
        end
    end
end
function d = bdiff(a)
    if isempty(a)
        d = [];
    else
        n = numel(a);
        if n == 1
            d = a;
        else
            d = a(2:end,:)-a(1:end-1,:);
        end
    end
end
function V = normalize(V)
% normalize the vector A (1D or 2D)
    m = min(min(V));
    M = max(max(V));
    V = (V-m)./(M-m);
end
function p = polyfit_0_intercept(x,y,n)
    % polyfit function, tweaked to solve for 0 intercept    
    a = x.^(n:-1:1); % calc x^n matrix excluding x^0
    p = a\y; % solve    
    p = [p;0]'; % set also intercept to zero
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = randcoords_points(coords,density,r)

    % divide into 2 batches to compensate for dense areas
    d1 = floor(density/2);
    d2 = density-d1;
    
    % 1st batch
    [X,idx] = sort(coords(:,1));
    Y = coords(idx,2);
    [x1,nearest] = fx(X,d1,r);
    y1 = Y(nearest)+r/2*(2*rand(d1,1)-1);
    
    % 2nd batch
    [Y,idx] = sort(coords(:,2));
    X = coords(idx,1);
    [y2,nearest] = fx(Y,d2,r);
    x2 = X(nearest)+r/2*(2*rand(d2,1)-1);
    
    % union of batches
    x = [x1;x2];
    y = [y1;y2];
    [x,idx] = sort(x);
    y = y(idx);
end
function [Y,nearest] = fx(X,density,r)
    
    % locate X distances greater than r
    D = diff(X); 
    S = (D>r);
    
    % calculate singular points
    
    % boundaries
    left = X([true;S])-r/2;
    right = X([S;true])+r/2;
    s_length = (right-left);
    cum_s_length = [0;cumsum(s_length)];
    total_length = sum(s_length);

    % generate random distribution (output x axis)
%     y = cumsum(rand(density,1));
%     y = y./y(end)*total_length;
    
    y = sort(rand(density,1))*total_length;
    
    % calculate Y coords
    Y = injective(y,left,cum_s_length);

    % locate nearest X point
    n = numel(X);
    if n == 1
        nearest = ones(n,1);
    else
        nearest = interp1(X,(1:n)',Y,'nearest','extrap');
    end
%     figure(1);hold off;
%     scatter(x,X,'k'); hold on;
%     scatter(y,Y,'sb','Linewidth',3);
end
function y = injective(x,left,cs)
    % find indices of segment to which each x belongs to
    
    idx_cs = (cs'>=x);
    idx_left = sum(~idx_cs,2);
    y = left(idx_left)+x-(cs(idx_left));
end