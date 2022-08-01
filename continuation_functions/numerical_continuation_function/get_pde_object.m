function p=get_pde_object(dir,varargin)
% get_pde_object: get/load pde object from dir/fname.mat
% effective method to save pde information before/after branch swiching 
%
%  USAGE:
%  p=get_pde_object(dir)               - load file with largest label
%  p=get_pde_object(dir,fname)
%  p=get_pde_object(dir,fname,newname) - load and set problem name to newname
%
% Called by switch_branch.m

% Check the existence of input directory
if exist(dir,'dir')
    if isempty(varargin)
        bfname=char(['pt' mat2str(max(getlabs(dir)))]);
    else 
        bfname=varargin{1}; varargin=varargin(2:end);
    end
else 
    fprintf('Directory %s does not exist.\n',dir)
    return
end

% set up the name sytax of file
ffname=[dir '/' bfname '.mat']; 

% use load function 
try s=load(ffname,'p');
catch 
    fprintf('Point %s does not exist.\n',ffname); % error warning
    bfname=ptselect(dir); ffname=[dir '/' bfname '.mat']; 
    s=load(ffname,'p');
end

p=s.p; % pde object loaded successfully
m0=isfield(p,'mesh'); 
m1=isfield(p, 'pdeo');  
if m0==1; m0=isfield(p.mesh, 'geo'); end
if ((m0==0) && (m1==0)); try load(p.file.mname); p.mesh=m; catch; end; end 
if ~isempty(varargin); p=set_file_name(p,varargin{1}); end

p=setfemops(p); 
p.file.count=p.file.count+1;
