function solution = cavity_patternsearch_v1(structure,spectrum,reference,comparison,search_space,flag_plot)
% this function calculates the nk solution based on the structure, spectrum
% and reference provided, via pattern search method
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

    % set up patternsearch problem
    if nargin<4 || ~isnumeric(search_space) || search_space(1)<1
        search_space = 1;
    end
    if nargin<5 || ~islogical(flag_plot)
        flag_plot = false;
    end
    problem.Aineq = [ones(1,reference.N-1)];
    problem.bineq = [search_space];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = (1-search_space) * ones(1,length(reference.material)-1);
    problem.ub = search_space * ones(1,length(reference.material)-1);
    problem.nonlcon = [];
    problem.solver = 'patternsearch';
    problem.options = optimoptions(problem.solver,'UseParallel',true);

    % set up solution variable
    solution.N = reference.N-1;
    solution.p = ones(spectrum.N,solution.N);
    p = zeros(1,solution.N)-0.25;

    %% fit and solution calculation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % search - global
    idx = 1:spectrum.N;
    % initiate cavity
    cavity = create_cavity(structure,1:spectrum.N);
    % redefine solver problem objective function and x0
    problem.objective = @(p) search_fun_reflectance(p,cavity,idx,structure.investigated,reference,spectrum);
    problem.x0 = p;
    % calculate (can't use 'solve(problem)' as in documentation)
    solution.p = patternsearch(problem);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % calculate n and k solution from weight values
    solution.range = structure.range;
    solution.e = reduced_weighted_average(solution.p,reference.edata);
    solution.n = e2n(solution.e);
    solution.k = e2k(solution.e);
    
    % plot solution diagram
    if flag_plot
        f = figure('Name','solver','Units','normalized','Position',[.5 .5 .4 .4]);
        plot_weights(f,solution.p,reference);
        f = figure('Name','solution','Units','normalized','Position',[.1 .5 .4 .4]);
        plot_sols(f,solution,reference,comparison,structure,spectrum);
    end
    
end

function h = plot_weights(h,p,reference)
    if isempty(h)
        h = figure('Name','solver');
    else
        figure(h.Number);
    end
    delete(get(gca,'Children'));
    
    w = mean([p,1-sum(p,2)],1);
    
    % plot diagram
    theta = pi/2+(0:2*pi/reference.N:2*pi);
    polarplot(theta,ones(1,length(w)+1),'k:');
    hold on;
    polarplot(theta,[w,w(1)],'b','Linewidth',1);
      
    hold off;
    nticks = ceil(max(abs(w))*2)/2;
    if nticks<1
        nticks = 1;
    end
    rlim([0 nticks]);    
    rticks(0:.5:nticks); 

    [tks,idx,~] = unique(mod((90+180*(0:2/reference.N:2)),360));
    thetaticks(tks)
    
    labels = rm_(reference.material(idx));
    for ii = 1:numel(labels)
        labels{ii} = sprintf('%s, w = %0.3f',labels{ii},w(idx(ii)));
    end
    thetaticklabels(labels);
    drawnow
end
function h = plot_sols(h,solution,reference,comparison,structure,spectrum)

    if isempty(h)
        h = figure('Name','solution');
    else
        figure(h.Number);
    end
    delete(get(gca,'Children'));
    delete(gca);
    
    % constants
%     GREY = ([128,128,128])/255;
    LGREY = ([160 160 160])/255;
    BLACK = [0 0 0];
    RED = [255 8 8]/255;
    MAGENTA = [255,165,255]/255;
%     AZUL = [0 255 255]/255;
    BLUE = [70 70 255]/255;
%     WHITE = [1 1 1];
%     DGREEN = [50 126 41]/255;
%     LGREEN = [77 193 64]/255;
%     DORANGE = [222 118 3]/255 *.5;    
    Smooth_Constant = 40;
    LW = 3;
    dim_color = 1/2;


    a = axes('position',[.0613 .1833 .248 .6867]);
    delta = .33;
    b = axes('position',[.0613+delta .1833 .248 .6867]);
    delta = 2*delta;
    c = axes('position',[.0613+delta .1833 .248 .6867]);

    range = solution.range;
    p = [solution.p,1-sum(solution.p)];
    p = p-min(p); p = p./max(p).*dim_color+dim_color;
    
    % refractive index
    ng = reference.N+1+1;
    pn = gobjects(ng);
    hold(a,'on');
    for ii = 1:reference.N
        pn(ii) = plot(a,range,smooth(reference.ndata(:,ii),Smooth_Constant,'loess'),'linewidth',LW,'color',BLUE.*p(ii));
    end
    lg = reference.material;
    pn(ii+1) = plot(a,range,smooth(solution.n,Smooth_Constant,'loess'),'-','linewidth',LW,'color',RED);
    lg(reference.N+1) = {'Fitting'};
    if ~isempty(comparison.material) && ~isempty(comparison.material{1})     
        ii = ii+2;
        pn(ii) = plot(a,range,smooth(comparison.ndata,Smooth_Constant,'loess'),'--','LineWidth',LW,'color',MAGENTA);
        lg(ii) = {'Comparison'};
    end
    lg = rm_(lg);
    legend(a,lg,'Location','Best','fontweight','bold','fontname','Arial','fontsize',8','color',BLACK,'box','off')
    ylabel('n')
    xlabel('wavelength [m]');
    set(a,'fontweight','bold','fontname','Arial','fontsize',12,...
        'box','on','linewidth',1.5, 'XColor',LGREY,'YColor',LGREY);
    
    % extinction coefficient
    pn = gobjects(ng);
    hold(b,'on');
    for ii = 1:reference.N
        pn(ii) = plot(b,range,smooth(reference.kdata(:,ii),Smooth_Constant,'loess'),'linewidth',LW,'color',BLUE.*p(ii));
    end
    pn(ii+1) = plot(b,range,smooth(solution.k,Smooth_Constant,'loess'),'-','linewidth',LW,'color',RED);
    if ~isempty(comparison.material) && ~isempty(comparison.material{1})
        ii = ii+2;
        pn(ii) = plot(b,range,smooth(comparison.kdata,Smooth_Constant,'loess'),'--','LineWidth',LW,'color',MAGENTA);
    end
    legend(b,lg,'Location','Best','fontweight','bold','fontname','Arial','fontsize',8','color',BLACK,'box','off')
    ylabel('k')
    xlabel('wavelength [m]');
    set(b,'fontweight','bold','fontname','Arial','fontsize',12,...
        'box','on','linewidth',1.5, 'XColor',LGREY,'YColor',LGREY);
    
    % reflection
    pn = gobjects(ng);
    hold(c,'on');    
    for ii = 1:reference.N
        cavity = create_cavity(structure,reference.edata(:,ii),(1:numel(range))');
        pn(ii) = plot(c,range,TMM_fun_reduced(range,cavity.L,cavity.ER),'linewidth',LW,'color',BLUE.*p(ii));
    end
    cavity = mod_cavity(cavity,solution.e);
    pn(ii+1) = plot(c,range,TMM_fun_reduced(range,cavity.L,cavity.ER),'-','linewidth',LW,'color',RED);
    if ~isempty(comparison.material) && ~isempty(comparison.material{1})
        ii = ii+2;
        cavity = mod_cavity(cavity,comparison.edata);
        pn(ii) = plot(c,range,TMM_fun_reduced(range,cavity.L,cavity.ER),'--','LineWidth',LW,'color',MAGENTA);
    end
    pn(ii+1) = plot(c,range,spectrum.rdata,'-','LineWidth',LW,'color',BLACK);
    lg(ii+1) = {'Measured'};
    
    legend(c,lg,'Location','Best','fontweight','bold','fontname','Arial','fontsize',8','color',BLACK,'box','off')
    ylabel('k')
    xlabel('wavelength [m]');
    set(c,'fontweight','bold','fontname','Arial','fontsize',12,...
        'ylim',[0,1],'YTick',[0:.2:1],...
        'box','on','linewidth',1.5, 'XColor',LGREY,'YColor',LGREY);

    drawnow
end