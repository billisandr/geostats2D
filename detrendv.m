%   DESCRIPTION:
%
%	The function 'detrendv' detrends 2D data. [fluc,trend,trfunc,a,scores,
%	dfreq,a_trends] = detrendv(x,y,rf,nfr,model,flag) calculates the linear
%	coefficients of the chosen trend model for the given data 'rf' on the
%	grid [x,y]. The coefficients of the model are estimated by using linear
%	regression. The definition of the best model (amongst the supported) by
%	means of AIC is also available.
%
%	INPUT VARIABLES:
%
%   [x,y]	      Numerical rectangular grid of the field (no matter if
%                 vector or matrix)
%   rf:           Random field's values on the grid(same size as x)
%   nfr:          Number of dominant frequencies required (for periodic
%                 trend model)
%   model:        Desired Model (string)
% 	flag:	      If 1 plots are drawn.
%
%
%	OUTPUT VARIABLES:
%
%	fluc:	 Fluctuations on the [x,y] grid
%	trend:	 Trend on the [x,y] grid
%   trfunc:  Trend models' function
%   a:       Coefficients of trend model
%   scores:  Table of scores of the examined models
%   dfreq:   Estimated dominant frequencies of the periodic trend model
%            (matrix of size [nfr x 2], 1st column corresponds to x axis
%            and 2nd to y axis)
%   a_trend:  Cell containing the linear coefficients for all the models
%
%	COMMENTS:
%
%	The following models are currently supported:
%
%   'mean':          Mean value (0 degree linear)
%   'linear':        Linear  (1st degree linear)
%	'quadratic':     Quadratic	(2nd degree linear)
%	'cubic':	     Cubic  (3rd degree linear)
%	'quartic':	     Quartic  (4th degree linear)
%   'periodic':      Periodic
%   'best':          Estimates all the above and chooses the best by means
%                    of corrected Akaike
%
%   EXAMPLE:
%
%   [x,y,rf] = randomfield('Mate',[1,4,4,0,0,2],60,60,1);
%   x1 = x(:,3:3:60);
%   y1 = y(:,3:3:60);
%   rf1 = rf(:,3:3:60);
%   %rf2 = rf1 + (3 + 0.5*x1 + 1.4*y1);%add linear trend
%   rf2 = rf1 + 4*sin(2*pi*0.3*x1) + 6*sin(2*pi*0.2*y1);%add periodic trend
%   [fluc,trend,trfunc,a,scores,dfreq,a_trend] =...
%                      detrendv(x1,y1,rf2,2,'linear',1);
%
%   See also: randomfield  
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


function [fluc,trend,trfunc,a,scores,dfreq,a_trends] =...
    detrendv(x,y,rf,nfr,model,flag)

x1 = x(:); y1 = y(:); rf1 = rf(:);
N = numel(rf);

%---Mean Value---

% Trend Function
Mx_func_m =  @(x,y,a)ones(size(x,1),1)*a;

% Multiple regression estimation of data trend
X_m = ones(N,1);
[a_m,~,fluc_m,~,~] = regress(rf1,X_m);
% [a_m,bint_m,fluc_m,rint_m,stats_m] = regress(rf1,X_m);
% a_m = mean(rf1);

%Trend of data
Mx_m = Mx_func_m(x1,y1,a_m);

fluc0{1,1}=fluc_m;trend0{1,1}=Mx_m;trfunc0{1,1}=Mx_func_m;a0{1,1}=a_m;

if flag==1
    figure;
    normplot(fluc_m(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(1,1),ks_p(1,1),ks_ksstat(1,1),ks_cv(1,1)] = kstest(fluc_m); %Kolmogorov-Smirnov test scores
end

%---Linear Model---

% Trend Function
Mx_func_l =  @(x,y,a)[ones(size(x,1),1), x, y]*a;

% Multiple regression estimation of data trend
X_l = [ones(N,1) x1 y1];
[a_l,~,fluc_l,~,~] = regress(rf1,X_l);
% [a_l,bint_l,fluc_l,rint_l,stats_l] = regress(rf1,X_l);

%Trend of data
Mx_l = Mx_func_l(x1,y1,a_l);

fluc0{2,1}=fluc_l;trend0{2,1}=Mx_l;trfunc0{2,1}=Mx_func_l;a0{2,1}=a_l;

if flag==1
    figure;
    normplot(fluc_l(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(2,1),ks_p(2,1),ks_ksstat(2,1),ks_cv(2,1)] = kstest(fluc_l);
end

%---Quadratic Model---

% Trend Function
Mx_func_q =  @(x,y,a_q)[ones(size(x,1),1),x, y, x.*y, x.^2, y.^2]*a_q;

% Multiple regression estimation of data trend
X_q = [ones(N,1) x1 y1 x1.*y1 x1.^2 y1.^2];
[a_q,~,fluc_q,~,~] = regress(rf1,X_q);
% [a_q,bint_q,fluc_q,rint_q,stats_q] = regress(rf1,X_q);

%Trend of data
Mx_q = Mx_func_q(x1,y1,a_q);

fluc0{3,1}=fluc_q;trend0{3,1}=Mx_q;trfunc0{3,1}=Mx_func_q;a0{3,1}=a_q;

if flag==1
    figure;
    normplot(fluc_q(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(3,1),ks_p(3,1),ks_ksstat(3,1),ks_cv(3,1)] = kstest(fluc_q);
end

%---Cubic Model---

% Trend Function
Mx_func_c =  @(x,y,a_c)[ones(size(x,1),1), x, y, x.*y, x.^2, y.^2, x.^2.*y, x.*y.^2, x.^3, y.^3]*a_c;

% Multiple regression estimation of data trend
X_c = [ones(N,1) x1 y1 x1.*y1 x1.^2 y1.^2 x1.^2.*y1 x1.*y1.^2 x1.^3 y1.^3];
[a_c,~,fluc_c,~,~] = regress(rf1,X_c);
% [a_c,bint_c,fluc_c,rint_c,stats_c] = regress(rf1,X_c);

%Trend of data
Mx_c = Mx_func_c(x1,y1,a_c);

fluc0{4,1}=fluc_c;trend0{4,1}=Mx_c;trfunc0{4,1}=Mx_func_c;a0{4,1}=a_c;

if flag==1
    figure;
    normplot(fluc_c(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(4,1),ks_p(4,1),ks_ksstat(4,1),ks_cv(4,1)] = kstest(fluc_c);
end

%---Quartic Model---

% Trend Function
Mx_func_q4 =  @(x,y,a_q4)[ones(size(x,1),1), x, y, x.*y, x.^2, y.^2,...
    x.^2.*y, x.*y.^2, x.^3, y.^3, x.^2.*y.^2,x.^3.*y, x.*y.^3, x.^4, y.^4]*a_q4;

% Multiple regression estimation of data trend
X_q4 = [ones(N,1) x1 y1 x1.*y1 x1.^2 y1.^2 x1.^2.*y1 x1.*y1.^2 x1.^3 y1.^3 ...
    x1.^2.*y1.^2 x1.^3.*y1  x1.*y1.^3  x1.^4  y1.^4];
[a_q4,~,fluc_q4,~,~] = regress(rf1,X_q4);
% [a_q4,bint_q4,fluc_q4,rint_q4,stats_q4] = regress(rf1,X_q4);

%Trend of data
Mx_q4 = Mx_func_q4(x1,y1,a_q4);

fluc0{5,1}=fluc_q4;trend0{5,1}=Mx_q4;trfunc0{5,1}=Mx_func_q4;a0{5,1}=a_q4;

if flag==1
    figure;
    normplot(fluc_q4(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(5,1),ks_p(5,1),ks_ksstat(5,1),ks_cv(5,1)] = kstest(fluc_q4);
end

%---Periodic Model---

% >>>Estimate dominant frequencies on x and y direction<<<

% FFT of data and Magnitude
[ny,nx] = size(rf);
Y = fft2(rf-mean(rf(:))); % FFT of 2D data
% if min(min(Y)) == 0
%     Magn = log(Y + 1e-7);
% else
%     Magn = log(Y);
% end
Magn = abs(Y).^2; % magnitude of 2D data
% Magn = log(abs(Y).^2); % magnitude of 2D data

% Magnitude of non-zero frequencies on x direction
MagnX = Magn(:,2:floor(nx/2));
MagnX = MagnX(:);
freqx0 = 0:1/nx:(1-1/nx);
freqx = 1/nx:1/nx:(1/2-1/nx);
freqx = repmat(freqx,ny,1);
freqx = freqx(:);

% Magnitude of non-zero frequencies on y direction
MagnY = Magn(2:floor(ny/2),:);
MagnY = MagnY(:);
freqy0 = 0:1/ny:(1-1/ny);
freqy = 1/ny:1/ny:(1/2-1/ny);
freqy = repmat(freqy',1,nx);
freqy = freqy(:);

% Plot Magnitude vs Frequencies on x and y directions
figure
subplot(2,1,1)
plot(freqx0,Magn)
%title('x dimension in the Frequency Domain')
xlabel('Frequency')
ylabel('Magnitude')
subplot(2,1,2)
plot(freqy0,Magn)
%title('y dimension in the Frequency Domain')
xlabel('Frequency')
ylabel('Magnitude')

% Dominant frequencies
[~,indx] = sort(MagnX,'descend');
[~,indy] = sort(MagnY,'descend');
dfreqx0 = freqx(indx); [dfreqx0,~,~] = unique(dfreqx0,'stable');
dfreqy0 = freqy(indy); [dfreqy0,~,~] = unique(dfreqy0,'stable');
dfreqx = dfreqx0(1:nfr);
dfreqy = dfreqy0(1:nfr);
% dfreqx = freqx(indx(1:nfr));
% dfreqy = freqy(indy(1:nfr));
dfreq = [dfreqx,dfreqy];

% Trend Function
syms x y
comp_str = []; % components' string
for i=1:nfr
    frx = dfreqx(i); fry = dfreqy(i);
    comp_str0 = [cos(2*pi*frx*x), sin(2*pi*frx*x),cos(2*pi*fry*y),sin(2*pi*fry*y)];
    comp_str = [comp_str comp_str0];  %#ok<AGROW>
    clear comp_str0
end
f = @(x,y)subs(comp_str);
Mx_func_p =  @(x,y,a_p)eval([ones(size(x,1),1) x y f(x,y)])*a_p;

X_p = [];
for i=1:nfr
    X_p0 = [cos(2*pi*dfreqx(i)*x1) sin(2*pi*dfreqx(i)*x1) cos(2*pi*dfreqy(i)*y1) sin(2*pi*dfreqy(i)*y1)];
    X_p = [X_p X_p0]; %#ok<AGROW>
    clear X_p0
end
X_p = [ones(N,1) x1 y1 X_p];
[a_p,~,fluc_p,~,~] = regress(rf1,X_p);
% [a_p,bint_p,fluc_p,rint_p,stats_p] = regress(rf1,X_p);
%a_p = X_p\rf(:);

%Trend and fluctuation of data
Mx_p = Mx_func_p(x1,y1,a_p);

fluc0{6,1}=fluc_p;trend0{6,1}=Mx_p;trfunc0{6,1}=Mx_func_p;a0{6,1}=a_p;

if flag==1
    figure;
    normplot(fluc_p(:))
    h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
    [ks_h(6,1),ks_p(6,1),ks_ksstat(6,1),ks_cv(6,1)] = kstest(fluc_p);
end

%Information Criteria & Scores table
S.TrendModel = {'Mean';'Linear';'Quadratic';'Cubic';'Quartic';'Periodic'};
for i=1:length(S.TrendModel)
    S.k(i,1) = length(a0{i,1});
    S.ks_h(i,1) = ks_h(i,1);
    S.ks_p(i,1) = ks_p(i,1);
    S.ks_ksstat(i,1) = ks_ksstat(i,1);
    S.ks_cv(i,1) = ks_cv(i,1);
    S.sq_error(i,1) = sum(fluc0{i,1}.^2);
    S.AIC(i,1) = log(S.sq_error(i,1)/N) + ((N+2*S.k(i,1))/N);
    S.BIC(i,1) = log(S.sq_error(i,1)/N) + ((S.k(i,1)*log(N))/N);
    S.AICc(i,1) = log(S.sq_error(i,1)/N) + ((N+S.k(i,1))/(N-S.k(i,1)-2));
end
a_trends = a0;
scores = struct2table(S);

switch model
    case 'best'
        %Best Model Definition (Squared Error or Acaike)
        % best_score = max(S.sq_error);
        % best_score = min(S.sq_error);
        best_score = min(S.AICc);
        % best_score = max(S.ks_p);
        idx = find (best_score==S.AICc);
        idx = idx(1);
        
        best_trendmodel = S.TrendModel{idx,1}   %#ok<NASGU,NOPRT>
        fluc = fluc0{idx,1};
        trend = trend0{idx,1};
        trfunc = trfunc0{idx,1};
        a = a0{idx,1};
        
    case 'mean'
        
        fluc = fluc0{1,1};
        trend = trend0{1,1};
        trfunc = trfunc0{1,1};
        a = a0{1,1};
        
    case 'linear'
        
        fluc = fluc0{2,1};
        trend = trend0{2,1};
        trfunc = trfunc0{2,1};
        a = a0{2,1};
        
    case 'quadratic'
        
        fluc = fluc0{3,1};
        trend = trend0{3,1};
        trfunc = trfunc0{3,1};
        a = a0{3,1};
        
    case 'cubic'
        
        fluc = fluc0{4,1};
        trend = trend0{4,1};
        trfunc = trfunc0{4,1};
        a = a0{4,1};
        
    case 'quartic'
        
        fluc = fluc0{5,1};
        trend = trend0{5,1};
        trfunc = trfunc0{5,1};
        a = a0{5,1};
        
    case 'periodic'
        
        fluc = fluc0{6,1};
        trend = trend0{6,1};
        trfunc = trfunc0{6,1};
        a = a0{6,1};
        
    otherwise
        
        error('Unrecognizable Model')
end



if flag==1
    
    %Plot data & trend
    X = x1; Y = y1;
    figure;
    scatter3(X,Y,rf1,'filled')
    hold on
    Xfit = min(X):1:max(X); nxfit = length(Xfit);
    Yfit = min(Y):1:max(Y); nyfit = length(Yfit);
    [XFIT,YFIT] = meshgrid(Xfit,Yfit);
    VFIT = trfunc(XFIT(:),YFIT(:),a);
    VFIT = reshape(VFIT,nyfit,nxfit);
    
    mesh(XFIT,YFIT,VFIT)
    colorbar
    xlabel('x')
    ylabel('y')
    %xlabel('Alongside Section')
    %ylabel('Depth')
    %title('Data & Trend Model')
    view(-120,9)
    shading interp
    colorbar
    
    %     figure;
    %     gridDelaunay = delaunay(X,Y);
    %     trimesh(gridDelaunay,X,Y,rf)
    %     hold on
    %     gridDelaunay = delaunay(X,Y);
    %     trimesh(gridDelaunay,X,Y,trend)
    %     shading flat
    %     az = 100;el = 30;
    %     view(az, el);
    %     xlabel X; ylabel Y;
    %     title('Data & Trend Model')
    %     clb = colorbar;
    %     clb.Label.String = 'Grey Value';
    
end

end