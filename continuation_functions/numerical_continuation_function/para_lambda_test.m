function p=para_lambda_test(p)
% para_lambda_test: test if an entry in p.usrlam or lammin or lammax has been passed; 
%  if so, locate 
%
%    p=para_lambda_test(p)
%
% (if more than one such entry skipped at once, only the 1st is detected)
%
% Called by function: find_bifurcation_point.m

if p.sw.para==3; return; end 
targets=[p.usrlam p.nc.lammin p.nc.lammax]; ind=[];
if(size(p.branch,2)>1)
    for j=1:length(targets)
        if(sign(p.branch(4,end)-targets(j)) ~= sign(p.branch(4,end-1)-targets(j)))
            ind=j; break
        end
    end
end
try
if(~isempty(ind))
    lamd=p.tau(p.nu+p.nc.nq+1); 
    ds=(targets(ind)-p.branch(4,end))/lamd;
    au1=transform_u_activeu(p,p.u,1)+ds*p.tau; u1=transform_activeu_u(p,au1,1); % predictor
    ntime=tic;  [u1,res,iter,Gu]=newton_corrector(p,u1); sol.meth='natural';
    p.time.newton=p.time.newton+toc(ntime); % newton-loop time, accumulated
    if(res<=p.nc.tol) %if convergence in Newton   
        obs=1; % since 2019 matlab complains about "if 0", hence stupid construction 
        if ~obs
        if(mod(p.file.count,p.file.smod)==0) % if last continuation point was saved, 
            p.file.count=p.file.count+1; % resave it with a label increased by 1
            p.branch(1,end)=p.branch(1,end)+1; 
            p.fuha.savefu(p); p.file.count = p.file.count-1;  
            p.branch(1,end)=p.branch(1,end)-1;            
        end
        end 
        q=p; q.u=u1; q.sol.ptype=-3;    
        if q.sw.spcalc>0 % calculate EVals
            sptime=tic; 
            [q.sol.ineg,q.sol.muv]=vspcalc(Gu,q);
            p.time.spec=p.time.spec+toc(sptime); % spectral-time, accumulated
        else 
            q.sol.ineg=-1;
        end
        if(q.sw.errcheck>0); q.sol.err=errcheck(q);end   
        brout=[default_branch_data(q); p.fuha.outfu(q,q.u)];  
        % userfu to append to bif-branches
        brplot=brout(length(bradat(q))+p.plot.bpcmp); %y-axis value in bif-figure        
        fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
            p.file.count,printcon(get_para_lambda(q)),printcon(brplot),res,iter, sol.meth,printcon(ds));

        if(q.sw.errcheck>0)
            fprintf('%5.2e ',q.sol.err);
        end

        if(q.sw.spcalc==1)
            fprintf(' %2i ', q.sol.ineg); 
        end

        npr=length(q.sw.bprint);
        for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(q))+p.sw.bprint(i)))); end
        fprintf('\n');
        % put on branch (the new point ends up between the last and the one but last)
        p.branch=[p.branch(:,1:end-1) brout p.branch(:,end)];  
        p.branch(1,end)=p.branch(1,end)+1; 
        q.branch=p.branch(:,1:end-1); q.file.count=p.file.count;
        p.fuha.savefu(q); % save the solution at the desired lambda to the file with label p.file.count
        p.file.count=p.file.count+1;
        
        %plot the point
        figure(p.plot.brfig); hold on;      % plot point with according symbol
        if q.sol.ineg<=0; plot(get_para_lambda(q),real(brplot),'*'); 
        else
            plot(get_para_lambda(q),real(brplot),'+'); 
        end         
    else % no convergence 
        fprintf('no convergence in newton_corrector for lambda = %g \n', targets(ind))
        %do nothing (temporary solution)
    end
end
catch 
end 

