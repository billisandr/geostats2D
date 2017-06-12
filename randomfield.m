%   DESCRIPTION:
%
%	The function 'randomfield' constructs gaussian 2D random fields
%	on a grid. [x,y,rf,Cxx]=randomfield(covmod,beta,nx,ny,flag) returns a
%	2D random field rf on the grid [x,y].
%
%	The field is generated with multiplication of the covariance matrix
%	defined from the input arguments with a random vector (sampled from
%	gaussian distribution).
%   The anisotropy is inserted by calculating the covariance on the rotated
%   and rescaled coordination system.
%
%	INPUT VARIABLES:
%
%   covmod:	Covariance model (string)
%   beta:   Covariance model's parameters [s2,xi1,xi2,phi,c0,v or eta1]
%	nx:		Size of the grid in the x-direction (integer)
%	ny:		Size of the grid in the y-direction (integer)
%	flag:	If flag~=1 or flag=[] the plots are suppressed.
%           Default value is 0.
%
%
%	OUTPUT VARIABLES:
%
%	[x,y]:	Numerical rectangular grid (both of size [nx x ny])
%	rf:		Random field values on the grid (size [nx x ny])
%   Cxx:    Covariance matrix of the field (size [nx x ny])
%
%	COMMENTS:
%
%	Five covariance models are currently supported:
%
%   'Gexp': Generalized Exponential
%   'Sphe': Spherical
%	'Gaus':	Gaussian
%	'Mate':	Matern ++'Mate1/3','Mate0.5','Mate1.0','Mate1.5','Mate2.0',
%                    'Mate2.5','Mate3.0','Mate3.5'
%   'Spar': Spartan 3D ++'Spar1','Spar2','Spar3'
%
%    SPARTAN model is a new covariance function family introduced
%            by Hristopulos D. T. and Elogne S.N. (2007). Analytic
%            properties and covariance functions for a new class of
%            generalized gibbs random fields. IEEE Transactions on
%            Information Theory, 53(12):4667–4679
%
%   EXAMPLE1:
%
%   To generate and display a 2D isotropic random field with Matern
%   covariance, zero nugget effect and smoothness index v=2 use the
%   following:
%
%   [x,y,rf,Cxx]=randomfield('Mate',[1,4,4,0,0,2],60,60,1);
%
%   EXAMPLE2:
%
%   To generate and display a 2D anisotropic random field with Matern
%   covariance and smoothness index v=2 use the following:
%
%   [x,y,rf,Cxx]=randomfield('Mate',[1,4,1.5,30,0,2],60,60,1);
%
%   Code written by V.Androulakis and A.Pavlides 
%   Last update: 28 April, 2017
%
%==========================================================================
% 
%   Copyright (C) 2016, 2017 V. Androulakis and A. Pavlides
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


function [x,y,rf,Cxx]=randomfield(covmod,beta,nx,ny,flag)

% Check number of model's parameters
cmod = {'Gexp';'Gaus';'Sphe';'Mate';'Mate1/3';'Mate0.5';'Mate1.0';...
    'Mate2.0';'Mate2.5';'Mate3.0';'Mate3.5';'Spar';'Spar1';'Spar2';...
    'Spar3'};
npar0 = [6;5;5;6;5;5;5;5;5;5;5;5;6;6;5;6];
npar = npar0(strcmp(covmod,cmod)==1,1);

if size(beta,2)> npar
    warning('Size of input models parameters is bigger than the necessary')
elseif size(beta,2)<npar
    error('Error: Missing models parameters.')
end

% Initialize parameters

s2 = beta(1);

xi1 = beta(2);
xi2 = beta(3);
phi = beta(4)*pi/180;

c0 = beta(5);
if npar==6
    v = beta(6);
else
    v = [];
end

if nargin==4 || flag~=1
    flag = 0;
end

% Transformation Matrix (for anisotropic fields)
A = [cos(phi)/xi1, sin(phi)/xi1;
    -sin(phi)/xi2, cos(phi)/xi2];
% A = [cos(phi)/xi1, -sin(phi)/xi2;
%       sin(phi)/xi1, cos(phi)/xi2];

% Grid construction
rx = 1:nx;
ry = 1:ny;
%rx = linspace(1,20,nx);
%ry = linspace(1,20,ny);
[x,y] = meshgrid(rx,ry);

% Distance Matrices
X = reshape(x,1,nx*ny);
Y = reshape(y,1,nx*ny);
C = [X;Y];
Ctr = A*C; %transformed coordinations
xtr = Ctr(1,:);
ytr = Ctr(2,:);
distx = pdist2(xtr',xtr'); % distance matrix in the x-direction
disty = pdist2(ytr',ytr'); % distance matrix in the y-direction

% Anisotropic norm of distance
h = sqrt(distx.^2 + disty.^2);

switch covmod
    case 'Gexp'
        tit = sprintf('Gener. Exponential covariance model (v=%g)',v);
        Cxx = s2*(exp(-h.^v));
        
    case 'Gaus'
        tit = 'Gaussian covariance model';
        Cxx = s2*(exp(-h.^2));
        
    case 'Sphe'
        tit = 'Spherical covariance model';
        Cxx = s2*(1-1.5*h+0.5*(h.^3)).*(1-heaviside2(h-1));
        
    case 'Mate'
        tit = sprintf('Matern covariance model (v=%g)',v);
        Cxx = s2*(1/((2^(v-1))*gamma(v)))*(h.^v).*(besselk(v,h));
        
    case 'Mate1/3' % von Karman: v = 1/3
        tit = 'Matern model with v=1/3';
        Cxx  = s2*(((2^(2/3))/gamma(1/3))*(h.^(1/3)).*(besselk(1/3,h)));
        
    case 'Mate0.5' % Exponential: v = 1/2
        tit = 'Matern model with v=1/2';
        Cxx = s2*exp(-h);
        
    case 'Mate1.0'   % Whittle: v = 1.0
        tit = 'Matern model with v=1.0';
        Cxx = s2*(h.*besselk(1,h));
        
    case 'Mate1.5' % Modified Exponential: v = 1.5
        tit = 'Matern model with v=1.5';
        Cxx = s2*((1 + h).*exp(-h));
        
    case 'Mate2.0' % Matern with v = 2.0
        tit = 'Matern model with v=2.0';
        Cxx = s2*(0.5*(h.^2).*besselk(2,h));
        
    case 'Mate2.5' % Third-order Autoregressive: v = 2.5
        tit = 'Matern model with v=2.5';
        Cxx = s2*((1 + h + (1/3)*h.^2).*exp(-h));
        %Cxx =  (1/((2^(1.5))*gamma(2.5)))*(h.^2.5).*(besselk(2.5,h));
        
    case 'Mate3.0' % Matern with v = 3.0
        tit = 'Matern model with v=3.0';
        Cxx = s2*(((h.^3)/8).*besselk(3,h));
        
    case 'Mate3.5' % Matern with v = 3.5
        tit = 'Matern model with v=3.5';
        Cxx = s2*(((1/15)*(h.^3) + 0.4*(h.^2) + h + 1).*exp(-h));
        %Cxx = ((2^(-2.5))/gamma(3.5))*(h.^3.5).*besselk(3.5,h) ;
        
    case 'Spar'
        tit = sprintf('Spartan covariance model (eta1=%g)',v);
        eta0 = s2; eta1 = v;
        delta = sqrt(abs(eta1^2-4));
        if abs(eta1)< 2 %##3-dimensional covariance + use the form s2*correl(h) to avoide inaccuracies##
            beta1 = 0.5*(sqrt(2-eta1));
            beta2 = 0.5*(sqrt(2+eta1));
            s2 = eta0/(4*pi*sqrt(2+eta1));
            Cxx = s2*((exp(-h*beta1)).*(sin(h*beta2)./(h*beta2)));
        elseif  eta1==2
            s2 = eta0/(8*pi);
            Cxx = s2*exp(-h);
        elseif eta1>2
            omega1 = sqrt((eta1-delta)/2);
            omega2 = sqrt((eta1+delta)/2);
            s2 = eta0/(4*pi*delta); %(omega2-omega1) is missing knowingly
            Cxx = s2*((exp(-h*omega1) - exp(-h*omega2))./h);
        else
            disp('Permissibility conditions fall');
        end
        
    case 'Spar1'
        eta0 = s2; eta1 = v;
        beta1 = 0.5*(sqrt(2-eta1));
        beta2 = 0.5*(sqrt(2+eta1));
        s2 = eta0/(4*pi*sqrt(2+eta1));
        Cxx = s2*((exp(-h*beta1)).*(sin(h*beta2)./(h*beta2)));
        tit = 'Spartan model with {|\eta_1|} < 2';
    case 'Spar2'
        eta0 = s2;
        s2 = eta0/(8*pi);
        Cxx = s2*exp(-h);
        tit = 'Spartan model with {\eta_1} = 2';
    case 'Spar3'
        eta0 = s2; eta1 = v;
        delta = sqrt(abs(eta1^2-4));
        omega1 = sqrt((eta1-delta)/2);
        omega2 = sqrt((eta1+delta)/2);
        s2 = eta0*(omega2-omega1)/(4*pi*delta);
        Cxx =s2*((exp(-h*omega1) - exp(-h*omega2))./(h*(omega2-omega1)));
        tit = 'Spartan model with {\eta_1} > 2';
        
    otherwise
        error('Unrecognizable Model')
end

Cxx(distx==0 & disty==0) = c0 + s2;

n_posdef = all(eig(Cxx)>0); %check for positive definite Cxx

if n_posdef==0
    error('Non positive definite covariance matrix')
else
    B = sqrtm(Cxx);
    
    u1=randn(nx*ny,1); %random vector
    
    simRF=B*u1; %simulated random field
    rf = reshape(simRF,nx,ny);
end

if flag==1
    figure;
    surf(x,y,rf);
    shading interp;
    view(2);
    colorbar
    title(sprintf('Gaussian Random Field with %s',tit));
end
end


