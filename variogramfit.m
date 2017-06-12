%   DESCRIPTION:
%
% 	The function 'variogramfit' estimates the parameters of a desired
% 	variogram model which give the best fitting of the theoretical model to
% 	the experimental variogram.The estimation is achieved by minimizing an
% 	specified error.
% 
% 	INPUT VARIABLES:
%   
%   gexp:	   Experimental variogram (vector of size [nrbins,1]) 
%   nr_pairs:  Number of pairs belong to each bin. (same size as 'gexp')
%   c_centers: Structure which contains the centers of the distance bins
%              (c_centers.distc), the centers of the angular bins 
%              (c_centers.phic) and the respective cartesian coordinates
%              (c_centers.c).
%   iso:       1 for isotropic, 0 for anisotropic and 00 for anisotropic2
% 	model:	   Structure containing:
%              1. model.function: the desired model's name(string)
%              2. model.params0:the initial values of the model's 
%              parameters for minimization
%              3. model.paramslb:the lower bounds for the parameters
%              4. model.paramsub:the upper bounds for the parameters
% 	objmod:	   Objective function for the minimization (string)
% 	flag:	   If flag=1 the plots are drawn, differently they are
% 	           suppressed.
%           
% 	OUTPUT VARIABLES:
% 
% 	bmodel:	Model's Parameters
% 	fval:   Objective function's value at the last iteration
%   tit:    Model's Name
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
%   To fit a Matern variogram model to the experimental
%   variogram of a constructed random field use the following:
% 
%   [x,y,rf]=randomfield('Mate',[8,3,1.2,36,0.2,1.5],60,60,1); 
%   [gexp, nr_pairs, c_centers] = expvar(x,y,rf,0,0.7,25,20,20,1);
%   model.function = 'Mate1.5'; model.params0 = [5,2,2,10,0.1];
%   model.paramslb = [eps,eps,eps,-90,eps]; 
%   model.paramsub = [20,20,20,90,0.4];
%   [bmodel,fval,tit]=...
%          variogramfit(gexp,nr_pairs,c_centers,0,model,'NWEr_m',1)
%   
%   See also: randomfield,ccvfun,expvar
%
%   Code written by V.Androulakis
%   Last update: 28 April,2017
%
%==========================================================================
% 
%   Copyright (C) 2016, 2017 V. Androulakis
%   This file is part of geostats2D.
% 
%   geostats2D is free software: you can redistribute it and/or
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
                     

function [bmodel,fval,tit]=...
    variogramfit(gexp,nr_pairs,c_centers,iso,model,objmod,flag)

% Number of model's parameters
modstr = {'Gexp';'Gaus';'Sphe';'Mate';'Mate1/3';'Mate0.5';'Mate1.0';...
    'Mate1.5';'Mate2.0';'Mate2.5';'Mate3.0';'Mate3.5';'Spar';'Spar1';'Spar2';'Spar3'};
npar0_is = [4;3;3;4;3;3;3;3;3;3;3;3;4;4;3;4];
npar0_an = npar0_is + 2;

if iso==1    
    npar = npar0_is(strcmp(model.function,modstr)==1,1);
elseif iso==0 || iso==00
    npar = npar0_an(strcmp(model.function,modstr)==1,1);
else
    error('Input iso must be 1,0 or 00.')
end

% Initialise Parameters
model1 = model.function;
c = c_centers.c; %bins' centers coordinations
distc = repmat(c_centers.distc,1,numel(c_centers.phic));
d_NaN = ~isnan(gexp);
gexpmax = max(max(gexp));

% Variogram Function
vrg = @(b,c) ccvfun(iso,model1,b,c,'vrg');
[~,~,~,~,~,~,tit,~] = ccvfun(iso,model1,zeros(1,npar),{0,0});

% Initial Parameters and Limits
bmod0 = model.params0;
bmod_lb = model.paramslb;
bmod_ub = model.paramsub;

% Objective function

switch objmod
    
    case 'VarErr' % Variance of error
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2));
        
    case 'NEr' % Normalised Error
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2)./...
            (vrg(bmod,c).^2)); 
        
    case 'WEr_d'
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2)./...
            gexp_nrp);
        
    case 'WEr_m'
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2).*...
            gexp_nrp);
        
    case 'NWEr_d'
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2)./...
            gexp_nrp./(vrg(bmod,c).^2));
        
    case 'NWEr_m'
        objfun0 = @(bmod,c,gexp,gexp_nrp)(((vrg(bmod,c)- gexp).^2).*...
            gexp_nrp./(vrg(bmod,c).^2)); 
        
    otherwise
        error('Error: Unregognised objective function')
                   
end

% >>>Error Minimization and Parameters Definition<<<

% Initial values and matrices preallocation
options = optimoptions('fmincon','Algorithm','interior-point');
% % opts = optimoptions('fmincon','Algorithm','interior-point',...
% %     'MaxIter',10000,'TolCon',0.0001,'TolFun',0.0001,'Display','iter'); %,'UseParallel','always');
% % options = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1000000,'MaxIter',300000);
% % ms = MultiStart('StartPointsToRun','bounds-ineqs','Display','none'); %construct a MultiStart object
% % N_ms = 40;
% % gs = GlobalSearch('StartPointsToRun','bounds-ineqs');%construct a Global Search object

% Minimization
c1 = {c{1,1}(d_NaN),c{1,2}(d_NaN)};
gexp1 = gexp(d_NaN);
gexp_nrp1 = nr_pairs(d_NaN);

objfun = @(bmod) sum(objfun0(bmod,c1,gexp1,gexp_nrp1));
[bmodel,fval] = fmincon(objfun,bmod0,[],[],[],[],bmod_lb,bmod_ub,[],options);
% % problem = createOptimProblem('fmincon','objective',objfun,...
% %     'x0',bmod0,'lb',bmod_lb,'ub',bmod_ub,'nonlcon',[],'options',options);
% % 
% % %Run Solver
% % tic;
% % [bmodel, fval] = run(ms,problem,N_ms);
% % %[bmodel, fval] = run(gs,problem);
% % toc


% Plots
if flag==1 && iso==1
    if isfield(c_centers,'phic')==1
        tit1 = sprintf('%s \nDirectional Experimental (Semi-)Variogram on {\\phi = %g}^{\\circ}',tit,c_centers.phic*180/pi);
    else
        tit1 = sprintf('%s \nExperimental (Semi-)Variogram',tit);
    end
    figure;
    scatter(distc(d_NaN),gexp(d_NaN),'r*');
    hold on;
    %plot(rc,gexp(d_NaN),'LineWidth',2);
    c2 = {0:0.5:max(distc(d_NaN)),0*(0:0.5:max(distc(d_NaN)))};
    gth = vrg(bmodel,c2);
    plot(c2{1,1},gth,'b--');
    axis([0 max(distc (d_NaN)) 0 gexpmax*1.1]);
    %title(tit1);
    xlabel('x');
    ylabel('\gamma(x)');

    
elseif flag==1 && (iso==0 || iso==00)
    
    %h = @(b,c) hfun(iso,b,c);
    nrphibins = size(gexp,2);
    
    for i=1:nrphibins
        if isfield(c_centers,'phic')==1
            tit1 = sprintf('%s \nDirectional Experimental (Semi-)Variogram on {\\phi = %g}^{\\circ}',tit,c_centers.phic(1,i)*180/pi);
        else
            tit1 = sprintf('%s \nExperimental (Semi-)Variogram',tit);
        end
        figure;
        %h1 = sqrt(c1{1,1}(i,1).^2 + c1{1,2}(i,1).^2);
        %h1 = h(bmodel,{c{1,1}(:,i),c{1,2}(:,i)});
        scatter(distc(d_NaN(:,i),i),gexp(d_NaN(:,i),i),'r*');
        hold on;
        c2 = {[0;c{1,1}(d_NaN(:,i),i)],[0;c{1,2}(d_NaN(:,i),i)]};
        %m = max(max(c{1,1}(:,i)),max(c{1,2}(:,i)));
        %c2 = {0:0.5:2*m, (0:0.5:2*m)};
        gth = vrg(bmodel,c2); %gth(1) = c0(bmodel);
        plot([0;distc(d_NaN(:,i),i)],gth,'b--');
        axis([0 max(distc(d_NaN(:,i),i)) 0 gexpmax*1.1]);
        %title(tit1);
        xlabel('r');
        ylabel('\gamma(r)');
    end
    
end
     
end
