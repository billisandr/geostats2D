%   DESCRIPTION:
%
%   The 'ccvfun' routine estimates correlation, covariance and 
%   (semi-)variogram of various models. Choice between isotropic (iso==1) 
%   and anisotropic(iso==0), and between models is enabled. 
% 
% 	INPUT VARIABLES:
%
%   iso:       If 1 isotropic is used, if 0 anisotropic with xi1,xi2,phi 
%              anisotropy parameters is used, and if 00 anisotropic with 
%              xi1,R, phi anisotropy parameters is used.
%	model:     Model (string).
%	bmodel:    Model's parameters [(s2,xi,c0,v)or(s2,xi1,xi2,phi,c0,v) or
%	           (eta0,xi,c0,eta1)or(eta0,xi1,xi2,phi,c0,eta1)]
%   c:         Distances: 1x2 cell (1x1--> rx, cell 1x2-->ry)
%   nout:      String which defines whether the first output(out1) is the
%              correlation(nout='cor'),the covariance(nout='cov') or the
%              variogram(nout='vrg'). Default is 'cor'. The last three
%              "out" variables are always h,s2,c0.
%
%	OUTPUT VARIABLES:
%
%	cor:              Correlation vector
%	covar:            Covariance vector
%   variogr:          Variogram vector
%   h:                Normalised distance lag
%   s2:               Variance
%   c0:               Nugget effect
%   tit0,tit:         Model's general and more specific name, respectively
%
%	COMMENTS:
%
%	Eight covariance models are currently supported:
%
%   'Gexp': Generalized Exponential
%   'Sphe': Spherical
%   'Cubi': Cubic
%	'Gaus':	Gaussian
%   'Quad': Rational Quadratic
%   'Card': Cardinal Sine
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
%
%   EXAMPLE:
%
%   To calculate the correlation, covariance and (semi-)variogram on given 
%   coordinates with Matern model, zero nugget effect and smoothness 
%   index v=2 use the following:
%   
%   iso=0; model='Mate'; bmodel=[1,4,1.5,30,0,2];  
%   x = [1:10]; y = [1:2:20]; c = {pdist2(x',x'),pdist2(y',y')};
%   [out1,out2,out3,out4,out5,out6,tit0,tit] = ccvfun(iso,model,bmodel,c)
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

function [out1,out2,out3,out4,out5,out6,tit0,tit] =...
            ccvfun(iso,model,bmodel,c,nout)

%Check number of models' given parameters
modstr = {'Gexp';'Gaus';'Sphe';'Cubi';'Quad';'Card';'Mate';'Mate1/3';'Mate0.5';'Mate1.0';...
    'Mate1.5';'Mate2.0';'Mate2.5';'Mate3.0';'Mate3.5';'Spar';'Spar1';'Spar2';'Spar3'};
npar0_is = [4;3;3;3;4;3;4;3;3;3;3;3;3;3;3;4;4;3;4];
npar0_an = npar0_is + 2;

if iso==1    
    npar = npar0_is(strcmp(model,modstr)==1,1);
elseif iso==0 || iso==00
    npar = npar0_an(strcmp(model,modstr)==1,1);
else
    error('Input iso must be 1,0 or 00.')
end
    
if length(bmodel)> npar
    warning('Size of input covariance parameters is bigger than the necessary')
elseif length(bmodel)<npar
    error('Missing covariance parameters.')
end


% Initialize parameters

if iso==1
    
    s2 = bmodel(1);
    xi1 = bmodel(2);
    xi2 = xi1;
    phi = 0;
    c0 = bmodel(3);
    if npar==4
        v = bmodel(4);
    else
        v = [];
    end
    
    x = c{1,1};
    if size(c,2)==1
        y = zeros(size(c{1,1}));
    else
        y = c{1,2};
    end
        
elseif iso==0
    
    s2 = bmodel(1);
    xi1 = bmodel(2);
    xi2 = bmodel(3);
    phi = bmodel(4)*pi/180;
%     if xi1>xi2
%         xi1 = bmodel(3);
%         xi2 = bmodel(2);
%         if phi>0
%             phi = (bmodel(4)-90)*pi/180;
%         elseif phi<0
%             phi = (bmodel(4)+90)*pi/180;
%         end
%     end
    c0 = bmodel(5);
    if npar==6
        v = bmodel(6);
    else
        v = [];
    end 
    
    x = c{1,1}; y = c{1,2};
    
elseif iso==00
    
    s2 = bmodel(1);
    xi1 = bmodel(2);
    R = bmodel(3);
    phi = bmodel(4)*pi/180;
%     if R>1
%         R = 1/R;
%         if phi>0
%             phi = (bmodel(4)-90)*pi/180;
%         elseif phi<0
%             phi = (bmodel(4)+90)*pi/180;
%         end
%     end
    xi2 = xi1/R;
    c0 = bmodel(5);
    if npar==6
        v = bmodel(6);
    else
        v = [];
    end  
    
    x = c{1,1}; y = c{1,2};
    
end

% Normalised Lag
A1 = (cos(phi)/xi1)^2 + (sin(phi)/xi2)^2;
A2 = (cos(phi)/xi2)^2 + (sin(phi)/xi1)^2;
A12 = (1/(xi1^2)- 1/(xi2^2))*sin(phi)*cos(phi);
h = sqrt(A1*(x.^2) + A2*(y.^2) + 2*A12*(x.*y));

switch model
    
    case 'Gexp'
        cor = exp(-h.^v);
        tit0 = 'Generalised Exponential model';
        if isa(bmodel,'sym')==1
            tit = tit0;
        else
            tit = sprintf('Generalised Exponential model (v = %.02g)',v);
        end
    case 'Gaus'
        cor = exp(-h.^2);
        tit0 = 'Gaussian model'; tit = tit0;
    case 'Sphe'
        cor = (1-1.5*h+0.5*(h.^3));
        cor(h>1) = 0;
        tit0 = 'Spherical model';tit = tit0;
    case 'Cubi'
        cor = 1 - 7*(h.^2) + 35/4*(h.^3) - 7/2*(h.^5) + 3/4*(h.^7);
        cor(h>1) = 0;
        tit0 = 'Cubic model';tit = tit0;
    case 'Quad'
        cor = (1 + h.^2).^(-v);
        tit0 = 'Rational Quadratic model';tit = tit0;
    case 'Card'
        cor = sin(h)./h;
        tit0 = 'Cardinal Sine model';tit = tit0;    
    case 'Mate'
        cor = (1/((2^(v-1))*gamma(v)))*(h.^v).*(besselk(v,h));
        tit0 = 'Generalised Matern model';
        if isa(bmodel,'sym')==1
            tit = tit0;
        else
             tit = sprintf('Matern model with v = %.02g',v); 
        end
    case 'Mate1/3' % von Karman: v = 1/3
        cor = ((2^(2/3))/gamma(1/3))*(h.^(1/3)).*(besselk(1/3,h));
        tit0 = 'Matern model with v = 1/3'; tit = tit0;
    case 'Mate0.5' % Exponential: v = 1/2
        cor = exp(-h);
        tit0 = 'Matern model with v = 0.5';tit = tit0;
    case 'Mate1.0'   % Whittle: v = 1.0
        cor = h.*besselk(1,h);
        tit0 = 'Matern model with v = 1.0';tit = tit0;
    case 'Mate1.5' % Modified Exponential: v = 1.5
        cor = (1 + h).*exp(-h);
        tit0 = 'Matern model with v = 1.5';tit = tit0;
    case 'Mate2.0' % Matern with v = 2.0
        cor = 0.5*(h.^2).*besselk(2,h);
        tit0 = 'Matern model with v = 2.0';tit = tit0;
    case 'Mate2.5' % Third-order Autoregressive: v = 2.5
        cor = (1 + h + (1/3)*h.^2).*exp(-h);
        %cor = (1/((2^(1.5))*gamma(2.5)))*(h.^2.5).*(besselk(2.5,h));
        tit0 = 'Matern model with v = 2.5';tit = tit0;
    case 'Mate3.0' % Matern with v = 3.0
        cor = ((h.^3)/8).*besselk(3,h);
        tit0 = 'Matern model with v = 3.0';tit = tit0;
    case 'Mate3.5' % Matern with v = 3.5
        %cor = ((2^(-2.5))/gamma(3.5))*(h.^3.5).*besselk(3.5,h) ;
        cor = ((1/15)*(h.^3) + 0.4*(h.^2) + h + 1).*exp(-h);
        tit0 = 'Matern model with v = 3.5';tit = tit0;
    case 'Spar'
        eta0 = s2; eta1 = v;
        delta = sqrt(abs(eta1^2-4));
        if abs(eta1)< 2
            beta1 = 0.5*(sqrt(2-eta1));
            beta2 = 0.5*(sqrt(2+eta1));
            s2 = eta0/(4*pi*sqrt(2+eta1));
            cor = (exp(-h*beta1)).*(sin(h*beta2)./(h*beta2));
        elseif  eta1==2
            s2 = eta0/(8*pi);
            cor = exp(-h);
        elseif eta1>2
            omega1 = sqrt((eta1-delta)/2);
            omega2 = sqrt((eta1+delta)/2);
            s2 = eta0*(omega2-omega1)/(4*pi*delta);
            cor = (exp(-h*omega1) - exp(-h*omega2))./(h*(omega2-omega1));
        else
            disp('Permissibility conditions fall');
        end
        tit0 = 'Spartan model';
        if isa(bmodel,'sym')==1
            tit = tit0;
        else
            tit = sprintf('Spartan model with {\\eta_1} = %.02g',v);
        end
        case 'Spar1'
        eta0 = s2; eta1 = v;
        beta1 = 0.5*(sqrt(2-eta1));
        beta2 = 0.5*(sqrt(2+eta1));
        s2 = eta0/(4*pi*sqrt(2+eta1));
        cor =(exp(-h*beta1)).*(sin(h*beta2)./(h*beta2));
        tit0 = 'Spartan model with {|\eta_1|} < 2';tit = tit0;
    case 'Spar2'
        eta0 = s2;
        s2 = eta0/(8*pi);
        cor = exp(-h);
        tit0 = 'Spartan model with {\eta_1} = 2';tit = tit0;
    case 'Spar3'
        eta0 = s2; eta1 = v;
        delta = sqrt(abs(eta1^2-4));
        omega1 = sqrt((eta1-delta)/2);
        omega2 = sqrt((eta1+delta)/2);
        s2 = eta0*(omega2-omega1)/(4*pi*delta); 
        cor =(exp(-h*omega1) - exp(-h*omega2))./(h*(omega2-omega1));
        tit0 = 'Spartan model with {\eta_1} > 2';tit = tit0;
    otherwise
        error('Unrecognizable Model')
end

cor(h==0) = 1;
correl = cor/(1 + c0/s2);
correl(h==0) = 1;
covar = s2*cor;
covar(h==0) = s2 + c0;
variogr = c0 + s2*(1 - cor);
variogr(h==0) = c0;

if nargin==4
    
    out1 = correl;
    out2 = covar;
    out3 = variogr;
    out4 = h;
    out5 = s2;
    out6 = c0;
    
elseif nargin==5
    
    if strcmp(nout,'cor')
        out1 = correl;
        out2 = covar;
        out3 = variogr;
        out4 = h;
        out5 = s2;
        out6 = c0;
    elseif strcmp(nout,'cov')
        out1 = covar;
        out2 = correl;
        out3 = variogr;
        out4 = h;
        out5 = s2;
        out6 = c0;
    elseif strcmp(nout,'vrg')
        out1 = variogr;
        out2 = correl;
        out3 = covar;
        out4 = h;
        out5 = s2;
        out6 = c0;
    elseif strcmp(nout,'h')
        out1 = h;
        out2 = s2;
        out3 = c0;
        out4 = correl;
        out5 = covar;
        out6 = variogr;
    elseif strcmp(nout,'s2')
        out1 = s2;
        out2 = h;
        out3 = c0;
        out4 = correl;
        out5 = covar;
        out6 = variogr;
    elseif strcmp(nout,'c0')
        out1 = c0;
        out2 = h;
        out3 = s2;
        out4 = correl;
        out5 = covar;
        out6 = variogr;
    else
        error('Error: Not valid input for nout')
    end
    
end


end