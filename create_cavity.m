function cavity = create_cavity(structure,investigated_data,idx)
    if nargin<3
        idx = ( 1:size(structure.edata,1) )';
    end
    cavity.ER = structure.edata(idx,:);
    cavity.L = structure.thickness;
    cavity.idx = idx;
    cavity.investigated = structure.investigated;
    if nargin==3
        cavity.ER(:,cavity.investigated) = investigated_data(idx);
    end
end