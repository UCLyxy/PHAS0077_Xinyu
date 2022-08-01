function plotsol(varargin) 
% plotsol: plot component of p.u in struct p; 
% plotsol(p,varargin) 

if ischar(varargin{1})
    dir=varargin{1};
    if nargin>1 && ischar(varargin{2}) % check if varargin{2} is a specific point in folder
    pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
    if exist(str,'file')==2
        p=load_point(dir,varargin{2}); anf=3;  %folder and point is given
    else  % check if user tried to load point which does not exist
         % (something like plotsol('h','pt4444')), or if user gives only folder and wants to
         % load the last point in the folder (call like plotsol('h','levn',3))
        if strcmp(pt(1:2),'pt') || strcmp(pt(1:3),'fpt') || strcmp(pt(1:3),'bpt')
            try p=load_point(dir,pt);  anf=3; catch; return; end
        else; p=load_point(dir);anf=2;
        end
    end
    else p=load_point(dir); anf=2; %user called something like plotsol('h',3,1,2) or plotsol('h')
    end % if nargin>1 && ischar(varargin{2})
    fprintf('lam=%g\n',get_para_lambda(p)); % show lambda in the work space
else p=varargin{1}; anf=2; pt=' '; % if first entry is a structure
end
%check existence of fields in p.plot and set fig. nr., component nr, and pstyle 
if isfield(p.plot,'pfig'); wnr=p.plot.pfig; else wnr=1; end
if isfield(p.plot,'pcmp'); cnr=p.plot.pcmp; else cnr=1; end
try; if cnr==0; cnr=1; end; end  % I do not why but sometime p.plot.pcmp=0 for saved sol.
if isfield(p.plot,'pstyle'); pstyle=p.plot.pstyle; else pstyle=1; end
% check if fig. nr., component nr, and pstyle are passed in form plotsol(p,11,1,2)
if nargin>=anf && isa(varargin{anf},'double'); wnr=varargin{anf}; end;
if nargin>=anf+1 && isa(varargin{anf},'double') && isa(varargin{anf+1},'double'); 
    cnr=varargin{anf+1};
end
if nargin>=anf+2 && isa(varargin{anf},'double') && isa(varargin{anf+1},'double')...
        && (isa(varargin{anf+2},'double')) 
    pstyle=varargin{anf+2};
end
% read out general options
try shsw=p.plot.shsw; catch shsw=1; end % shading-switch 
sub=0; 
for k=1:nargin
 if ischar(varargin{k})
   switch lower(varargin{k})
     case 'pfig'; wnr=varargin{k+1}; % figure number
     case 'pcmp'; cnr=varargin{k+1}; % component
     case 'pstyle'; pstyle=varargin{k+1}; % style
     case 'fs'; p.plot.fs=varargin{k+1}; % fontsize
     case 'cm'; p.plot.cm=varargin{k+1}; % colormap
     case 'axis'; p.plot.axis=varargin{k+1}; % axis (e.g.: 'tight', 'image')
     case 'sh'; shsw=varargin{k+1};  % shading
     case 'sub'; sub=varargin{k+1};  % subplot 
   end
 end
end
if pstyle==-1; userplot(p,wnr); return; end 
% set fig. nr. and u
figure(wnr); if sub~=0; subplot(sub(1),sub(2),sub(3:end)); end; cla 
try upde=p.mat.fill*p.u(1:p.nu); catch; upde=p.u; end 
n0=(cnr-1)*p.np+1; n1=cnr*p.np; u=upde(n0:n1); 
try; tname=[p.file.pname,num2str(p.file.count-1)]; catch tname=[]; end 
if p.sol.ptype==1 && pt(1)=='b'; tname=[p.file.bpname,num2str(p.file.bcount-1)]; end
if p.sol.ptype==2 && pt(1)=='f'; tname=[p.file.fpname,num2str(p.file.fcount-1)]; end
[po,~,~]=getpte(p);
switch size(po,1)
  case 1  %%%%%%%%%%%%%%%%%%% 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:nargin  % read out extra 1D-options    
      if ischar(varargin{k})
        switch lower(varargin{k})
          case 'lw'; p.plot.lw=varargin{k+1}; case 'cl'; p.plot.cl=varargin{k+1};
         end
      end
    end
    if isfield(p.plot,'lw')==0; p.plot.lw=1; end
    if isfield(p.plot,'cl')==0; p.plot.cl=zeros(length(cnr),3); end
    psc=cell(1,length(pstyle)); % plot style cell
    if isa(pstyle,'double')
      for i=1:length(pstyle)
        switch pstyle(i)
          case 1; psc{i}='-'; case 2; psc{i}='--'; case 3; psc{i}=':'; 
          case 4; psc{i}='*'; case 5; psc{i}='-*';
        end
      end
    else % convert pstyle to the cell {pstyle} 
      if isa(pstyle,'cell')==0; psc={pstyle}; else psc=pstyle; end
    end               
    if iscell(p.plot.cl); cl=zeros(length(p.plot.cl),3); % set colors
     for i=1:length(p.plot.cl); cl(i,:)=pde2pathcolors(p.plot.cl{i}); end
        p.plot.cl=cl;
    end
    if ischar(p.plot.cl); p.plot.cl=pde2pathcolors(p.plot.cl); end      
  %%%%%%%%%%%%%%%% begin filling %%%%%%%%%%%%%%%%%%%%%%%

    if length(cnr)>length(psc); dum=cell(1,length(cnr));
       for i=1:length(cnr); dum{i}=psc{1}; end; psc=dum;
    end 
  % if user types plotsol(p,11,[1,2],[1,2],'lw',3) 
  % instead of plotsol(p,11,[1,2],[1 2],'lw',[3,3]), we fill p.plot.lw
    if length(cnr)>length(p.plot.lw); dum=zeros(1,length(cnr));
      for i=1:length(cnr); dum(i)=p.plot.lw(1); end; p.plot.lw=dum;
    end 
  % if user types plotsol(p,11,[1,2],'cl','b3') instead of 
  % plotsol(p,11,[1,2],'lw',{'b3','b3'}), we fill p.plot.lw
    if length(cnr)>size(p.plot.cl,1); cl=zeros(length(cnr),3);
      for i=1:length(cnr); cl(i,:)=p.plot.cl(1,:); end; p.plot.cl=cl;
    end %%%%%%%% end filling, now plot solution(s)
   % clf(wnr); 
    for i=1:length(cnr);
       n0=(cnr(i)-1)*p.np+1; n1=cnr(i)*p.np; 
       p.pdeo.grid.plot(upde(n0:n1),psc{i},'linewidth',p.plot.lw(i),...
               'color',p.plot.cl(i,:));
    end
    drawnow; str=[]; % create string like u1, u2 
    if p.nc.neq==1; str='u'; % if we have one component system 
    else
      for i=1:length(cnr)-1; 
          str=[str,'\color[rgb]{',num2str(p.plot.cl(i,:)),'}u',num2str(cnr(i)),', ']; 
      end;
      if cnr==1; str=[str,'u',num2str(cnr(end))];  %use color black for u for 1 comp. plots
      else; str=[str,'\color[rgb]{',num2str(p.plot.cl(end,:)),'}u',num2str(cnr(end))]; % use different colors for u1 and u2         
      end
    end
    set(gca,'FontSize',p.plot.fs); title([str,'\color{black} at ',tname],'interpreter','tex');
    xlabel('x'); axis tight; box on; % set title, box, and label     

end 
drawnow 
if isfield(p.plot,'fs')==0; p.plot.fs=16; end
if p.nc.neq==1; t1=['u at ' tname]; % if we have a one component system
else t1=['u',num2str(cnr),' at ' tname]; 
end
% write title
set(gca,'fontsize',p.plot.fs); 
if size(po,1)~=1 % title writing for 1D solutions is done in the 1 D part above  
if size(po,1)==3 && pstyle==2; title({t1,t},'fontsize',p.plot.fs); % write point name and levels
else title(t1); % write point name
end
end
try; if p.plot.nola==1; nola; end; catch; end 
end % function plotsol ending

function anz=iszero(arg)% internal function, check if a field is zero (anz=1) or not (anz=0)
anz=0;
if isa(arg,'double');if max(size(arg))==1; if arg==0;anz=1; end; end; end
end