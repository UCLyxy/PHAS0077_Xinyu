function p=load_point(dir,varargin)
% load_point: load point from dir/fname.mat, for plotting
%
% USAGE:
%  p=load_point(dir)               - load file with largest label
%  p=load_point(dir,fname)
%  p=load_point(dir,fname,newname) - load and set problem name to newname
%
% Called by function:

% Check existence of directory
if exist(dir,'dir') 
    bfname=varargin{1}; 
    varargin=varargin(2:end);
else 
    fprintf('Directory %s does not exist.\n',dir); 
    return
end

ffname=[dir '/' bfname '.mat']; 
%load point from directory
try s=load(ffname,'p');

catch 
    fprintf('Point %s does not exist.\n',ffname); 
     bfname=ptselect(dir); 
     ffname=[dir '/' bfname '.mat']; 
     s=load(ffname,'p');
end

p=s.p; 
m0=isfield(p,'mesh'); 
m1=isfield(p, 'pdeo');

if m0==1
    m0=isfield(p.mesh, 'geo'); 
end

if ((m0==0) && (m1==0)); try load(p.file.mname); p.mesh=m; catch; end; end 
if ~isempty(varargin); p=set_file_name(p,varargin{1}); end
% set operators that were not saved in p
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
p.file.count=p.file.count+1;
if(p.sol.ptype==1) p.file.bcount=p.file.bcount+1; end
if(p.sol.ptype==2) p.file.fcount=p.file.fcount+1; end
