%   DESCRIPTION:
%
% 	The function 'ordkrig' estimates the values at the unknown points
% 	of a random field from the known data and a specified theoretical
% 	model. 
%   [Z,Z_error,checks,checks_mat] = ordkrig(x,y,rf,xu,yu,iso,model)
%   returns the estimated values and the estimations' error. The estimator
%   which is used is the ordinary kriging.
%
% 	INPUT VARIABLES:
%   
%   [x,y]:	  Numerical rectangular grid of the field (known points)
%   rf:       Random field values at the known points (stochastic
%             component)
%   [xu,yu]:  Numerical rectangular grid of the field (unknown points)
%   iso:      If 1 the correlation models used are isotropic(s2,xi1,c0,v),
%             if 0 they are anisotropic (s2,xi1,xi2,phi,c0,v) and if 00
%             anisotropic2(s2,xi1,R,phi,c0,v)
% 	model:	  Structure which contains:
%             1.model.function: the correlation function of the desired
%             model(string), 
%             2.model.params: the model's parameters
%             3.model.r_ok: the radius of ordinary kriging-->[rx_ok,ry_ok]
%
% 	OUTPUT VARIABLES:
% 
%   Z:            Estimations of random field at unknown points with trend
%                 addition and inversed boxcox transformation
%   Z_error:      Estimations' variance
%   checks:       Table of indices reffering to checka such as positive 
% 	              definiteness, symmetry, non-zero determinant of the 
%                 covariance matrices,negative Lagrange multiplier,etc.
%   checks_mat:   Structure containing the checks matrices.
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
% 
%   EXAMPLE:
% 
%   To estimate the values at the unknown points of a random field from
%   given known data and a specified theoretical model use the following:
% 
%   [x,y,rf] = randomfield('Mate1.0',[8,1.5,3,20,0.2],50,50,1);
%   figure;hist(rf(:));
%   [~,idx] = datasample(rf(:),floor(0.4*2500),'Replace',false);
%   idx = sort(idx);
%   x1 = x(idx);
%   y1 = y(idx);
%   rf1 = rf(idx);
%   rf2 = rf; rf2(idx) = 0;
%   [yu,xu] = find(rf2);
%   model.function = 'Mate1.0';
%   model.params = [8,1.5,3,20,0.2];
%   [Z,Z_error, checks, checks_mat] = ordkrig(x1,y1,rf1,xu,yu,0,model);
%   cc = [x1',y1',rf1';xu,yu,Z];
%   for i = 1:2500
%     rf_est(cc(i,2),cc(i,1)) = cc(i,3);
%   end
%   figure;pcolor(rf_est);shading interp
%   figure;hist(rf_est(:));
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


function [Z,Z_error,checks,checks_mat] = ordkrig(x,y,rf,xu,yu,iso,model)

[nx,ny] = size(rf);  
N = nx*ny; %number of known points
[nxu,nyu] = size(xu);  
Nu = nxu*nyu; %number of estimation points
c = [x(:),y(:)]; %known data coordinations' matrix
cu = [xu(:),yu(:)]; %unknown points coordinations' matrix
rf1 = rf(:);

% Ordinary Kriging Function
okr_func = @(la,z)sum(la.*z);

%Inputs
model1 = model.function; %model's function
bmodel = model.params; %model's parameters

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
[~,~,~,~,s2,c0,~,~] = ccvfun(iso,model1,bmodel,{0,0});

% Radius of kriging neighborhood
if isfield(model,'r_ok')==1
    r_ok = model.r_ok;
    neighb = 1; 
else
    neighb = 0;
end

if neighb==0
    
    %Distance matrices
    u = cu; %"unknown" point's coordinations
    s = c; %samples' coordinations
    rf2 = rf1; %samples' values of random field
    
    rx_nb_nb = pdist2(s(:,1),s(:,1)); %neighboors' distance matrix on the x axis (columns)
    ry_nb_nb = pdist2(s(:,2),s(:,2)); %neighboors' distance matrix on the y axis (rows)
    c_nb_nb = {rx_nb_nb,ry_nb_nb}; %distances between neigboors-neighboors point on x and y axes
    %r_nb_nb = pdist2(s,s); %neighboors' distance matrix
    
    rx_us(N,Nu)=0;ry_us(N,Nu)=0;
    for j=1:Nu
        rx_us(:,j) = pdist2(s(:,1),u(j,1)); %distance between samples-"unknown" point on the x axis (columns)
        ry_us(:,j) = pdist2(s(:,2),u(j,2)); %distance between samples-"unknown" point on the y axis (rows)
    end
    c_nb_pi = {rx_us,ry_us};%distances between neigbors-"unknown" point on x and y axes
    
    
     % Covariance matrices
     Cx_nb_nb = covar(bmodel,c_nb_nb); %covariance between neighboors
     Cx_nb_nb(rx_nb_nb==0 & ry_nb_nb==0) = s2 + c0; %correct values at zero distances
     Cx_nb_nb_exp = [[Cx_nb_nb ones(N,1)]; ones(1,N) 0]; %expand Cx_nb_nb
     Cx_nb_pi = covar(bmodel,c_nb_pi); %covariance between neighboors and "unknown" point
     Cx_nb_pi_exp = [Cx_nb_pi;ones(1,Nu)]; %expand Cx_nb_pi
     
     % Estimate weights
     wa = Cx_nb_nb_exp\Cx_nb_pi_exp; %linear system solution
     
     % Checks
     [~,p] = chol(Cx_nb_nb);
     posdef = (p==0); %check for positive definite covariance matrix
     %posdef = all(eig(Cx_nb_nb) > 0); %check for positive definite covariance matrix
     rcon = rcond(Cx_nb_nb_exp); %conditional number of expanded covariance matrix
     rweig = isreal(wa); %check for real weights
     
     Lagr(Nu,1)=0; Z(Nu,1)=0; Z_error(Nu,1)=0;
     for j=1:Nu
         % Value estimation at the unknown point
         Lagr(j,1) = wa(end,j);%Lagrange multiplier
         la = wa(1:end-1,j); %neighbors weights
         Z(j,1) = okr_func(la,rf2); %value estimation
         Z_error(j,1) = s2 + c0 - sum(la.*Cx_nb_pi(:,j)) - Lagr(j,1);  %minimum mean squared error (variance) of estimation
     end  
    
     N_nb = N;
     r_ok = 'whole domain';
     
elseif neighb==1
    
    % Matrices and Cells preallocation
    N_nb(Nu,1)=0; wa{Nu,1}=[];
    posdef(Nu,1)=0; rcon(Nu,1)=0;
    Lagr(Nu,1)=0; rweig(Nu,1)=0;
    Z(Nu,1)=0; Z_error(Nu,1)=0;
    
    % Definition of Weights (lambda)
    
    for i=1:Nu
        
        % "Unknown" point and "Samples" definition
        u = cu(i,:); %"unknown" point's coordinations
        s = c; %samples' coordinations
        rf2 = rf1; %samples' values of random field
        
        %Distance matrices and neighboors' definition
        rx_us = pdist2(s(:,1),u(1,1)); %distance between samples-"unknown" point on the x axis (columns)
        ry_us = pdist2(s(:,2),u(1,2)); %distance between samples-"unknown" point on the y axis (rows)
        %r_us = pdist2(s,u); % 2D distance marix
        II = (rx_us <= r_ok(1,1))&(ry_us <= r_ok(1,2)); %neighboors' index
        %II = r_us <= r_ok; %neighboors' index
        N_nb(i,1) = sum(sum(II)); %number of neighbors
        c_nb_pi = {rx_us(II),ry_us(II)};%distances between neigbors-"unknown" point on x and y axes
        %c_nb_pi = r_us(II); %distances between neigbors-"unknown" point
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
            r_nb_nb = pdist2(c_nb,c_nb); %neighboors' distance matrix
            %h_nb_nb = h(bmodel,c_nb); % %neighboors' normalised distance matrix
            
            % Covariance matrices
            Cx_nb_nb = covar(bmodel,c_nb_nb); %covariance between neighboors
            Cx_nb_nb(r_nb_nb==0) = s2 + c0; %correct values at zero distances
            Cx_nb_nb_exp = [[Cx_nb_nb ones(N_nb(i,1),1)]; ones(1,N_nb(i,1)) 0]; %expand cv_Cx_nb_nb
            Cx_nb_pi = covar(bmodel,c_nb_pi); %covariance between neighboors and "unknown" point
            Cx_nb_pi_exp = [Cx_nb_pi;1]; %expand cv_Cx_nb_pi
            
            % Estimate weights
            wa{i,1} = Cx_nb_nb_exp\Cx_nb_pi_exp; %linear system solution
            
            % Checks
            [~,p] = chol(Cx_nb_nb);
            posdef(i,1) = (p==0); %check for positive definite covariance matrix
            %posdef(i,1) = all(eig(Cx_nb_nb) > 0); %check for positive definite covariance matrix
            rcon(i,1) = rcond(Cx_nb_nb_exp); %conditional number of expanded covariance matrix
            rweig(i,1) = isreal(wa{i,1}); %check for real weights
            
            % Value estimation at the unknown point
            Lagr(i,1) = wa{i,1}(end,1);%Lagrange multiplier
            la = wa{i,1}(1:end-1,1); %neighbors weights
            Z(i,1) = okr_func(la,nb_rf2); %value estimation
            Z_error(i,1) = s2 + c0 - sum(la.*Cx_nb_pi) - Lagr(i,1);  %minimum mean squared error (variance) of estimation
            
        end
        
        clear u s rf2 rx_us ry_us c_us II c_nb_pi nb_rf2 rx_nb_nb ry_nb_nb c_nb r_nb_nb Cx_nb_pi Cx_nb_pi_exp Cx_nb_pi Cx_nb_pi_exp la
    end

end

%Structure with matrices of checks 
checks_mat.Lagr = Lagr;
checks_mat.posdef = posdef;
checks_mat.rcon = rcon;
checks_mat.wa = wa;
checks_mat.N_nb = N_nb;
checks_mat.s2 = s2;
checks_mat.r_ok = r_ok;
if neighb==0
    checks_mat.rx = rx_nb_nb;
    checks_mat.ry = ry_nb_nb;
    checks_mat.Cx = Cx_nb_nb;
end

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