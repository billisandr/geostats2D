%   DESCRIPTION:
% 
% 	The function 'aniso_dvf' estimates the anisotropy of a given random 
%   field. [R,phi,xi1,xi2,er1,er2]= ...
%   aniso_dvf(x,y,rf,model,objmod,ncpc,nrbins,phistep,phitol,flag)
%   returns the parameters of anisotropy (R,phi,xi1,xi2) of the random
%   field rf on the grid [x,y] and the errors of the variogram fitting and
%   the ellipses fitting (er1,er2).
% 
% 	The estimation of anisotropy is achieved by means of directional
% 	variograms with step=phistep and tolerance=phitol. Firstly, the
% 	experimental variogram is calculated. Then, the choosen theoretical
% 	model (in its isotropic form) is fitted to each directional variogram
% 	by minimizing the normalized error variance. Finally, an ellipsis is
% 	fitted to the estimated xis', which represents the anisotropy of the
% 	field. 
% 
%   'phi' represents the angle between the x axis of the coordinate
%   system and the first (not necessarily the major) principal axis
%   'R' represents the ratio of the correlation length along the first
%   principal axis over the correlation lengh along the orthogonal 
%   principal axis, i.e. R = xi1/xi2. 
% 
% 	INPUT VARIABLES:
%   
%   [x,y]:	  Numerical rectangular grid of the field 
%   rf:       Random field values on the grid
% 	model:	  Structure containing:
%             1. model.function: the desired model's name(string)
%             2. model.params0:the initial values of the model's parameters 
%             for minimization
%             3. model.paramslb:the lower bounds for the parameters
%             4. model.paramsub:the upper bounds for the parameters
% 	objmod:	  Objective function for the minimization (string)     
% 	ncpc:	  Percent of maximum distance of the grid taken into account to
% 	          the calculation of experimental variogram
%   nrbins:   Number of distance bins(integer)
%   phistep:  Anglular step of the direction variograms(degrees)
%   phitol:   Angular tolerance of the directional variograms(degrees)
% 	flag:	  If 1 the plots are drawn, differently they are suppressed.
%           
% 
% 	OUTPUT VARIABLES:
% 
% 	R:	  Anisotropy ratio
% 	phi:  Anisotropy orientiation angle (refers to xi1(degrees)
%   xi1:  Correlation length along the first principal axis 
%   xi2:  Correlation length along the orthogonal principal axis
%   er1:  Error of the directional variograms' fitting
%   er2:  Error of the anisotropy ellipses' fitting
%
% 	COMMENTS:
% 
% 	Five variogram models are currently supported:
% 
%   'Gexp':  Generalized Exponential
%   'Sphe':  Spherical
% 	'Gaus':  Gaussian		
% 	'Mate':  Matern	++'Mate1/3','Mate0.5','Mate1.0','Mate1.5','Mate2.0',
%                      'Mate2.5','Mate3.0','Mate3.5'
%   'Spar':  Spartan 3D ++'Spar1','Spar2','Spar3'
%
%    SPARTAN model is a new covariance function family introduced
%            by Hristopulos D. T. and Elogne S.N. (2007). Analytic
%            properties and covariance functions for a new class of
%            generalized gibbs random fields. IEEE Transactions on
%            Information Theory, 53(12):4667–4679
%   
%
%    Six objective functions are currently supported:
% 
%   'VarEr':   Variance of Error
%   'NEr':     Normalised Error
% 	'WEr_d':   Weighted Error (divide with number of pairs)		
% 	'WEr_m':   Weighted Error (multiply with number of pairs)
%   'NWEr_d':  Normalised and Weighted Error (divide with number of pairs)
%   'NWEr_m':  Normalised and Weighted Error (multiply with number of 
%              pairs)
%
%   EXAMPLE:
% 
%   To estimate the anisotropy parameters of a random field using Matern
%   covariance with smoothness index v=1.5 use the following:
%
%   %create random field
%   [x,y,rf]=randomfield('Mate',[1,3,1,36,0,2],60,60,1);
%   %Construct 'model' structure
%   model.function = 'Mate';
%   model.params0 = [0.5,2,0,1.5];
%   model.paramslb = [eps,eps,eps,1.5];
%   model.paramsub = [10,100,0.2,1.5];
%   %Estimate anisotropy
%   t6 = tic;
%   [R,phi,xi1,xi2,er1,er2]=aniso_dvf(x,y,rf,model,'NWEr_m',0.7,35,20,20,1)
%   T6 = toc(t6)
%   
%   See also: randomfield,variogramfit,expvar,ccvfun
%
%   Code written by V.Androulakis
%   Last update: 28 April,2017
%
%==========================================================================
% 
%   Copyright (C) 2016, 2017 V. Androulakis
%   This file is part of geostats2D.
% 
%   billisandr/geostats2D is free software: you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation, either version 3 of the
%   License, or any later version.
%
%   geostats2D is distributed in the hope that it will be 
%   useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
%   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%==========================================================================
                     

function [R,phi,xi1,xi2,er1,er2]=aniso_dvf(x,y,rf,model,objmod,ncpc,nrbins,phistep,phitol,flag)

iso = 0; %N_ms = 40;
[gexp, nr_pairs, c_centers] = expvar(x,y,rf,iso,ncpc,nrbins,phistep,phitol,flag);
nrphibins = size(gexp,2);

% >>> Theoretical Variogram Model Fitting <<<

% Minimization
npar = size(model.params0,2);
bmod(nrphibins,npar) = 0; fval(nrphibins,1) = 0; 
iso = 1;
niter = 0;

for j=1:nrphibins
    gexp1 = gexp(:,j);
    nr_pairs1 = nr_pairs(:,j);
    c1.c = {c_centers.c{1,1}(:,j),c_centers.c{1,2}(:,j)};
    c1.distc = c_centers.distc;
    c1.phic = c_centers.phic(1,j);
    [bmod(j,:),fval(j,1),tit]=variogramfit(gexp1,nr_pairs1,c1,iso,model,objmod,flag);
    
% %     %===============EXTRA STEP-REMOVE IF NOT NEEDED=====================%
% %     %2nd optimization
% %     model.params0 = bmod(j,:);
% %     [bmod(j,:),fval(j,1),tit]=variogramfit(gexp1,nr_pairs1,c1,iso,model,objmod,flag);
% %     %===============EXTRA STEP-REMOVE IF NOT NEEDED=====================%
    
    niter = niter+1 %#ok<*NOPRT>
end
    
er1 = sum(fval); %error of directional variograms' fitting

% >>> Elliptical Anisotropy  Model Definition <<<

% Initial Values
options = optimoptions('fmincon','Algorithm','interior-point');
ms = MultiStart('StartPointsToRun','bounds-ineqs','Display','none'); %construct a MultiStart object
% gs = GlobalSearch('StartPointsToRun','bounds-ineqs');%construct a Global Search object

phic = c_centers.phic;
th = [phic,phic+pi]';
r_thita = @ (A,th)(A(1)*A(2))./sqrt((A(2)*cos(th+A(3))).^2 + (A(1)*sin(th+A(3))).^2); %ellipsis function
fplot = 0:pi/1000:2*pi;
tit1 = sprintf('Elliptical Anisotropy Model of %s',tit); 

% Minimization and Ellipses Parameters Definition

%Remove outliers
range = repmat(bmod(:,2),2,1);
Nrg = size(bmod,1);
range2 = [bmod(Nrg-1:Nrg,2);range;bmod(1:2,2)];
d_out(2*Nrg,1)=0;
for i=1:2*Nrg
    d_out_m = mean([range2(i,1),range2(i+1,1),range2(i+3,1),range2(i+4,1)])
    d_out(i,1) = (range2(i+2,1)<1.3*d_out_m);%&&(range2(i+3,1)>0.7*d_out_m) ;
end    
d_out = logical(d_out);
% % q1 = quantile(range,0.25);%first quantile
% % q3 = quantile(range,0.75);%third quantile
% % iqr = q3 - q1; %interquantile range
% % d_out = range < q3 + 1.5*iqr; %outlier index
% % % nr_d_out = sum(d_out);

%A0 = [max(range),min(range),pi/4]; %a_el, b_el, phi
A0 = [max(range(d_out))/2,min(range(d_out))/2,pi/4]; %a_el, b_el, phi
bel_lb = [eps,eps,-pi/2]; bel_ub = [max(range(d_out))*2,max(range(d_out))*2,pi/2]; %lower and upper limits of fmincon

objfun_el = @(A_model)sum(((r_thita(A_model,th(d_out))-range(d_out)).^2));
%objfun_el = @(A_model)sum(((r_thita(A_model,th(d_out))-range(d_out)).^2)./...
%    (r_thita(A_model,th(d_out)).^2));

[A_mod, er2] = fmincon(objfun_el,A0,[],[],[],[],bel_lb,bel_ub,[],options);

% % problem = createOptimProblem('fmincon','objective',objfun_el,...
% %     'x0',A0,'lb',bel_lb,'ub',bel_ub,'nonlcon',[],'options',options);
% % 
% % % Run Solver
% % tic;
% % %[A_mod, er2] = run(gs,problem);
% % [A_mod, er2] = run(ms,problem,N_ms);
% % toc

% Definition of Anisotropy Angle and Ratio
fplot2 = -pi/2:pi/1000:pi/2;
r = r_thita(A_mod,fplot2);
phi_min = fplot2(1,r==min(r)); 
% phi_max = min(fplot(1,r==max(r)));

% phi = min(phi_min,phi_max);
% if phi == phi_min
%     xi1 = min(r);
%     xi2 = max(r);          % phi refers to the first principal axis (minor or major) 
% else
%     xi1 = max(r);
%     xi2 = min(r);
% end

phi = min(phi_min);    
xi1 = min(r);   % phi refers to the minor axis so as to have R<1
xi2 = max(r); 

R = xi1/xi2;
phi = phi*180/pi;

if flag==1
    %Plot ellipses
    figure;
    polar(fplot,r_thita(A_mod,fplot),'b--');
    %title(tit1);
    hold on
    polar(th(d_out),range(d_out),'r*')
end
      
end
