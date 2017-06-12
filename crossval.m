%   DESCRIPTION:
%
% 	The function 'crossval' estimates the cross validation scores of a
% 	specified theoretical model compared to a given set of data.
%   [cv_scores, checks, checks_mat] =...
%            crossval(x,y,rf,iso,d_col,model)
%   returns the cross validation scores of the estimated values compared to
%   the known data. The estimator which is used is the ordinary kriging. 
%
% 	INPUT VARIABLES:
%   
%   [x,y]:	  Numerical rectangular grid of the field 
%   rf:       Random field values on the grid
%   iso:      If 1 the correlation models used are isotropic, if 0 they are
%             anisotropic and if 00 anisotropic2
%   d_col:    If 1 the same column neighboors are excluded from the
%             estimation
% 	model:	  Structure which contains:
%             1.model.function: the correlation function of the desired
%             model(anonymous function of 'beta' and 'c = {rx,ry}'), 
%             2.model.params: the model's parameters
%             3.model.r_ok: the radius of ordinary kriging-->[rx_ok,ry_ok]
%           
% 
% 	OUTPUT VARIABLES:
% 
% 	cv_scores:	  Table containing cross validation scores
%   checks:       Table of indices reffering to checka such as positive 
% 	              definiteness, symmetry and non-zero determinant of the 
%                 covariance matrices, etc.
%   checks_mat:   Structure containing the check matrices... 
%
% 	COMMENTS:
% 
% 	Five variogram models are currently supported:
% 
%   'Gexp':   Generalized Exponential
%   'Sphe':   Spherical
% 	'Gaus':   Gaussian		
% 	'Mate':   Matern ++'Mate1/3','Mate0.5','Mate1.0','Mate1.5','Mate2.0',
%                      'Mate2.5','Mate3.0','Mate3.5'
%   'Spar':   Spartan 3D ( ++'Spar1','Spar2','Spar3'
%   'Spar1':  Spartan with |eta1|<2 
%   'Spar2':  Spartan with eta1=2
%   'Spar3':  Spartan with eta1>2 )
%
%    SPARTAN model is a new covariance function family introduced
%            by Hristopulos D. T. and Elogne S.N. (2007). Analytic
%            properties and covariance functions for a new class of
%            generalized gibbs random fields. IEEE Transactions on
%            Information Theory, 53(12):4667–4679
% 
%   EXAMPLE:
% 
%   To estimate the cross validation scores of a specified model compared
%   to a given Matern random field use the following:
% 
%   [x,y,rf]=randomfield('Mate0.5',[8,1.5,3,20,0.2],50,50,1);
%   model.function = 'Mate0.5';
%   model.params = [8,1.5,3,20,0.2];
%   model.r_ok = [3,3];
%   [cv_scores, checks, checks_mat] = crossval(x,y,rf,0,1,model);
%   
%   See also: randomfield,ccvfun
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


function [cv_scores, checks, checks_mat] =...
       crossval(x,y,rf,iso,d_col,model)

[nx,ny] = size(rf); %size of grid 
N = nx*ny;
c = [x(:),y(:)]; %coordinations' matrix
rf1 = rf(:);

model1 = model.function;
bmodel = model.params; %model's parameters

% Ordinary Kriging Function
okr_func = @(la,z)sum(la.*z);

%Check number of models' given parameters
modstr = {'Gexp';'Gaus';'Sphe';'Mate';'Mate1/3';'Mate0.5';'Mate1.0';...
    'Mate1.5';'Mate2.0';'Mate2.5';'Mate3.0';'Mate3.5';'Spar';'Spar1';'Spar2';'Spar3'};
npar0_is = [4;3;3;4;3;3;3;3;3;3;3;3;4;4;3;4];
npar0_an = npar0_is + 2;

if iso==1    
    npar = npar0_is(strcmp(model1,modstr)==1,1);
elseif iso==0 || iso==00
    npar = npar0_an(strcmp(model1,modstr)==1,1);
else
    error('Error: Input iso must be 1,0 or 00.')
end

if length(bmodel)> npar
    warning('Size of input covariance parameters is bigger than the necessary')
elseif length(bmodel)<npar
    error('Error: Missing covariance parameters.')
end

%Covariance Function
covar = @(b,c) ccvfun(iso,model1,b,c,'cov');
[~,~,~,~,s2,c0,~,tit] = ccvfun(iso,model1,bmodel,{0,0});

% Radius of kriging neighborhood
if isfield(model,'r_ok')==1
    r_ok = model.r_ok;
    neighb = 1; 
else
    neighb = 0;
end

if neighb==0
    
    %Distance matrices
    u = c; %"unknown" point's coordinations
    s = c; %samples' coordinations
       
    rx_nb_nb = pdist2(s(:,1),s(:,1)); %neighboors' distance matrix on the x axis (columns)
    ry_nb_nb = pdist2(s(:,2),s(:,2)); %neighboors' distance matrix on the y axis (rows)
    c_nb_nb = {rx_nb_nb,ry_nb_nb}; %distances between neigboors-neighboors point on x and y axes
        
    rx_us(N,N)=0;ry_us(N,N)=0;
    for j=1:N
        rx_us(:,j) = pdist2(s(:,1),u(j,1)); %distance between samples-"unknown" point on the x axis (columns)
        ry_us(:,j) = pdist2(s(:,2),u(j,2)); %distance between samples-"unknown" point on the y axis (rows)
    end
    c_nb_pi = {rx_us,ry_us};%distances between neigbors-"unknown" point on x and y axes
    
    
     % Covariance matrices
     Cx_nb_nb = covar(bmodel,c_nb_nb); %covariance between neighboors
     Cx_nb_nb(rx_nb_nb==0 & ry_nb_nb==0) = s2 + c0; %correct values at zero distances
     Cx_nb_nb_exp = [[Cx_nb_nb ones(N,1)]; ones(1,N) 0]; %expand Cx_nb_nb
     Cx_nb_pi = covar(bmodel,c_nb_pi); %covariance between neighboors and "unknown" point
     Cx_nb_pi_exp = [Cx_nb_pi;ones(1,N)]; %expand Cx_nb_pi
     
     % Matrices and Cells preallocation
     posdef(N,1)=0; rcon(N,1)=0;Lagr(N,1)=0;
     wa(N,N)=0; rweig(N,1)=0; Z(N,1)=0; Z_error(N,1)=0;
     
     for i=1:N
         
         rf2 = rf1; %samples' values of random field
         
         Cx_nb_nb_exp2 = Cx_nb_nb_exp;
         Cx_nb_nb_exp2(i,:) = [];  Cx_nb_nb_exp2(:,i) = []; 
         Cx_nb_pi_exp2 = Cx_nb_pi_exp(:,i);
         Cx_nb_pi_exp2(i) = [];
         
         % Estimate weights
         wa(:,i) = Cx_nb_nb_exp2\Cx_nb_pi_exp2; %linear system solution
         
         % Checks
         [~,p] = chol(Cx_nb_nb_exp2(1:end-1,1:end-1));
         posdef(i,1) = (p==0); %check for positive definite covariance matrix
         %posdef = all(eig(Cx_nb_nb) > 0); %check for positive definite covariance matrix
         rcon(i,1) = rcond(Cx_nb_nb_exp2); %conditional number of expanded covariance matrix
         rweig(i,1) = isreal(wa); %check for real weights
         
         % Value estimation at the unknown point
         Lagr(i,1) = wa(end,i);%Lagrange multiplier
         la = wa(1:end-1,i); %neighbors weights
         rf2(i) = [];
         Z(i,1) = okr_func(la,rf2); %value estimation
         Z_error(i,1) = s2 + c0 - sum(la.*Cx_nb_pi_exp2(1:end-1)) - Lagr(i,1);  %minimum mean squared error (variance) of estimation
         
         clear Cx_nb_nb_exp2 Cx_nb_pi_exp2 la rf2
     end
         
         N_nb = N-1;
         r_ok = 'whole domain';
     
elseif neighb==1  

    % Matrices and Cells preallocation
    N_nb(N,1)=0;
    posdef(N,1)=0; rcon(N,1)=0;Lagr(N,1)=0;
    wa{N,1}=[]; rweig(N,1)=0; Z(N,1)=0; Z_error(N,1)=0;
    
    % Definition of Weights (lambda)
    
    for i=1:N
        % "Unknown" point and "Samples" definition
        u = c(i,:); %"unknown" point's coordinations
        if d_col == 1
            idx = c(:,1)~=c(i,1); %samples' indices - points of the same column are excluded
        else
            idx = ones(N,1);
            idx(i,1) = 0;  %samples' indices - only the examined point is excluded
            idx = logical(idx);
        end
        s = c(idx,:); %samples' coordinations
        rf2 = rf1(idx,:); %samples' values of random field
        
        %Distance matrices and neighboors' definition
        rx_us = pdist2(s(:,1),u(1,1)); %distance between samples-"unknown" point on the x axis (columns)
        ry_us = pdist2(s(:,2),u(1,2)); %distance between samples-"unknown" point on the y axis (rows)
        II = (rx_us <= r_ok(1,1))&(ry_us <= r_ok(1,2)); %neighboors' index
        N_nb(i,1) = sum(II); %number of neighbors
        c_nb_pi = {rx_us(II),ry_us(II)};%distance between neigbors-"unknown" point on x and y axes
        c_nb = s(II,:); %neighboors' coordinates
        nb_rf2 = rf2(II); %neighbors' values of random field
        
        if N_nb(i,1)==0
            
            wa{i,1} = 0;
            Z(i,1) = 10^6;
            Z_error(i,1) = 10^6;
            warning('No neighboors found')
            
        else
            
            rx_nb_nb = pdist2(c_nb(:,1),c_nb(:,1)); %neighboors' distance matrix on the x axis (columns)
            ry_nb_nb = pdist2(c_nb(:,2),c_nb(:,2)); %neighboors' distance matrix on the y axis (rows)
            c_nb_nb = {rx_nb_nb,ry_nb_nb};
            
            % Covariance matrices & Linear solution
            Cx_nb_nb = covar(bmodel,c_nb_nb); %covariance between neighboors
            Cx_nb_nb(rx_nb_nb==0 & ry_nb_nb==0) = s2 + c0; %correct values at zero distances
            Cx_nb_nb_exp = [[Cx_nb_nb ones(N_nb(i,1),1)]; ones(1,N_nb(i,1)) 0]; %expand cv_Cx_nb_nb
            Cx_nb_pi = covar(bmodel,c_nb_pi); %covariance between neighboors and "unknown" point
            Cx_nb_pi_exp = [Cx_nb_pi;1]; %expand cv_Cx_nb_pi
            
            % Estimste weights
            wa{i,1} = Cx_nb_nb_exp\Cx_nb_pi_exp; %linear system solution
            
            % Checks
            [~,p] = chol(Cx_nb_nb); posdef(i,1) = (p==0); %check for positive definite covariance matrix
            %posdef(i,1) = all(eig(Cx_nb_nb{i,1}) > 0); %check for positive definite covariance matrix
            rcon(i,1) = rcond(Cx_nb_nb_exp); %conditional number of expanded covariance matrix
            rweig(i,1) = isreal(wa{i,1}); %check for real weights
            
            % Value estimation at the unknown point
            Lagr(i,1) = wa{i,1}(end,1);%Lagrange multiplier
            la = wa{i,1}(1:end-1,1); %neighbors weights
            Z(i,1) = okr_func(la,nb_rf2); %value estimation
            Z_error(i,1) = s2 + c0 - sum(la.*Cx_nb_pi) - Lagr(i,1);  %minimum mean squared estimation error
            
        end
        
        clear u idx s rf2 rx_us ry_us c_us II c_nb_pi nb_rf2 rx_nb_nb ry_nb_nb c_nb r_nb_nb Cx_nb_nb Cx_nb_nb_exp  Cx_nb_pi Cx_nb_pi_exp la
    end
    
end

% Cross Validation Scores

S.Model = {tit};
S.MeanAbsErr = mean(abs(rf1 - Z));
S.MaxAbsErr = max(abs(rf1 - Z));
S.MSE = mean((rf1 - Z).^2);
S.RMSE = sqrt(mean((rf1 - Z).^2));
S.rpearson = corr(Z,rf1);
S.rspearman = corr(Z,rf1,'type','Spearman');
cv_scores = struct2table(S);

%Structure with matrices of checks 
checks_mat.Lagr = Lagr;
checks_mat.posdef = posdef;
checks_mat.rcon = rcon;
checks_mat.wa = wa;
checks_mat.N_nb = N_nb;
checks_mat.Z = Z;
checks_mat.Z_error = Z_error;
checks_mat.s2 = s2;
checks_mat.r_ok = r_ok;

% Indices of checks table
N_posdef = all(all(posdef==1));
N_rcond = min(min(rcon));
N_Lagr = all(all(Lagr<=0));
N_rweig = all(all(rweig==1));
N_Z_error = all(all(Z_error>0));
table_h2 = {'posdef','rcond','Lagr','rweights', 'est_error'};
checks = [N_posdef, N_rcond, N_Lagr, N_rweig, N_Z_error];
checks = array2table(checks,'VariableNames',table_h2);

end