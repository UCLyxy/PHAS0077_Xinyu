function stanheadfu(p)
% STANHEADFU: standard "header function" for runtime output called at start
% of numerical_continuation.m
%
%  stanheadfu(p)
%
% Output: "step lambda      y-axis    residual iter meth   ds"
%
% means: step # on branch, value of primary bif. param., y-axis value in bif diagram, 
%    residual, # newton iteration, continuation method (normal/arclength), stepsize
%
% Called by numerical_continuation.m

auxdictl=0;

if(isfield(p.plot,'auxdict')) 
    auxdictl=length(p.plot.auxdict); 
end

if(auxdictl>=p.nc.ilam(1)) 
    primaux=p.plot.auxdict{p.nc.ilam(1)};
else 
    primaux='lambda';
end

primaux(primaux=='\')=''; 
primaux=[primaux blanks(12-length(primaux))];

if(p.plot.bpcmp==0)
    y_axis='L2-norm';
elseif(p.plot.bpcmp<=auxdictl) 
    y_axis=p.plot.auxdict{p.plot.bpcmp};
else 
    y_axis='y-axis'; 
end

y_axis(y_axis=='\')='';
try
    y_axis=[y_axis blanks(8-length(y_axis))]; 
catch
end

fprintf(char(['step   ' primaux y_axis 'residual  iter meth   ds       ']));

if(p.sw.errcheck>0) 
    fprintf('err      '); 
end

if(p.sw.spcalc==1) 
    fprintf('#-EV '); 
end

npr=length(p.sw.bprint);

for i=1:npr 
    if(p.sw.bprint(i)>0 && p.sw.bprint(i)<=auxdictl)
        addaux=char([p.plot.auxdict{p.sw.bprint(i)}]);
        addaux(addaux=='\')=''; fprintf(addaux);
        fprintf(blanks(11-length(addaux)));
    else 
        fprintf('b(%i)       ',p.sw.bprint(i)); 
    end
end
% put anything else here 
fprintf('\n');

