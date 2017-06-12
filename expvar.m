%   INSPIRATION:
%
%   The following code is inspired by, and makes use of, code published in
%   the file exchange by Wolfgang Schwanghart ("Experimental (Semi-)
%   Variogram" version 1.4)
%   (https://www.mathworks.com/matlabcentral/fileexchange/
%   20355-experimental--semi---variogram)
%
%   DESCRIPTION:
%
% 	The function 'expvar' estimates the experimental variogram of a given
% 	random field. [gexp, nr_pairs, c_centers] =...
%   expvar(x,y,rf,iso,ncpc,nrbins,phistep,phitol,flag) returns the values
%   of the experimental variogram ('gexp') calculated within bins whose
%   centers are defined by 'c_centers' and the number of the pairs which
%   belong to each bin ('nr_pairs').
% 
% 
% 	INPUT VARIABLES:
%   
%   [x,y]:      Numerical rectangular grid of the field (no matter if
%               matrix or vector)
%   rf:         Random field values on the grid (same size as input x)
%   iso:        If 1 the given random field is treated as isotropic. If   
%               0 then the directional experimental variogram is
%               calculated, by binning both distances and angles and
%               evaluating the experimental variogram within the
%               resulting cyclical sections.
% 	ncpc:	    Percent of maximum distance of the grid taken into
% 	            account to the calculation of experimental variogram
%   nrbins:     Number of distance bins(integer)
%   phistep:    Anglular step of the direction variograms(degrees)
%   phitol:     Angular tolerance of the directional variograms(degrees)
% 	flag:	    If 1 the plots are drawn.
%           
% 
% 	OUTPUT VARIABLES:
% 
% 	gexp:       Experimental variogram values
% 	nr_pairs:   Number of pairs belong to each bin
%   c_centers:  Coordinates of the bins' centers. Structure which contains 
%               the centers of the distance bins (c_centers.distc), the 
%               centers of the angular bins (c_centers.phic) and the 
%               respective cartesian coordinates (c_centers.c).
%
% 
%   EXAMPLE:
% 
%   To estimate the experimental variogram of an anisotropic random field
%   with Matern covariance, zero nugget effect and smoothness index v=2 use
%   the following:
% 
%   [x,y,rf] = randomfield('Mate',[1,3,1,36,0,2],60,60,1); 
%   [gexp, nr_pairs, c_centers] = expvar(x,y,rf,0,0.7,25,20,20,1)
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


function [gexp, nr_pairs, c_centers] =...
    expvar(x,y,rf,iso,ncpc,nrbins,phistep,phitol,flag)

[nx,ny] = size(rf); %size of grid
N = nx*ny;
c = [reshape(x,N,1),reshape(y,N,1)]; %coordinations' matrix
rf = reshape(rf,N,1);

%Distance Matrix
[k,l] = find(triu(true(N)));
d = hypot(c(k,1)-c(l,1),c(k,2)-c(l,2));
% d = triu(pdist2(c,c));
maxdist = max(max(d))*ncpc;

if iso==1
    
    %Distance binning
    disttol = maxdist/nrbins; %bin tolerance
    edges = linspace(0,maxdist,nrbins+1);edges(end) = inf; %distance bins
    distc = (edges(1:end-1)+ disttol/2)'; %centers of distance bins
    
    %Squared Value Differences
    I = d <= maxdist;
    distmat = [l(I) k(I) d(I)]; % distance matrix [point1 point2 distance]
    maxdist = max(distmat(:,3));
    gx = (rf(distmat(:,1))-rf(distmat(:,2))).^2; %squared differences
    
    
    % bin distance
    [~,ixedge] = histc(distmat(:,3),edges);
    % variogram anonymous function
    fvar = @(x) (1./(2*numel(x))) * sum(x); %we don't divide with 2 because we consider only the upper triangular matrix
    %calculate variogram
    gexp1 = accumarray(ixedge,gx,[numel(edges) 1],fvar,nan);
    gexp_nrp1 = accumarray(ixedge,ones(size(gx)),[numel(edges) 1],@sum,nan);
    gexp = gexp1(1:end-1,:);
    nr_pairs = gexp_nrp1(1:end-1,:);
    
    d_30 = nr_pairs <=30;
    gexp(d_30) = NaN;
    d_NaN = ~isnan(gexp);
    gexpmax = max(max(gexp));
    
    c_centers.c{1,1} = distc;% coordinates of  centers on x axis
    c_centers.c{1,2} = zeros(size(distc));% coordinates of  centers on y axis
    c_centers.distc = distc;
    c_centers.phic = 0;
    
    if flag == 1
        figure;
        plot(distc(d_NaN),gexp(d_NaN),'LineWidth',2);
        hold on;
        scatter(distc(d_NaN),gexp(d_NaN),'fill');
        axis([0 maxdist 0 gexpmax*1.1]);
        %title('Experimental (Semi-)Variogram');
        xlabel('x');
        ylabel('\gamma(x)');
    end
    
elseif iso==0
    
    %Distance and angle binning
    disttol = maxdist/nrbins; %bin tolerance
    phistep = phistep*pi/180; %angular step
    phitol = phitol*pi/180;  %angular tolerance
    
    edges = linspace(0,maxdist,nrbins+1);edges(end) = inf; %distance bins
    distc = (edges(1:end-1)+ disttol/2)'; %centers of distance bins
    
    phic = (0:phistep:pi-phistep); %centers of phi bins
    nrphibins = length(phic);
    
    % Squared Value Differences and Direction Matrix
    I = d <= maxdist;
    distmat = [l(I) k(I) d(I)]; % distance matrix [point1 point2 distance]
    maxdist = max(distmat(:,3));
    gx = (rf(distmat(:,1))-rf(distmat(:,2))).^2; %squared differences
    dirmat = atan2(c(distmat(:,2),2)-c(distmat(:,1),2),...
        c(distmat(:,2),1)-c(distmat(:,1),1)); % direction matrix
    
    %only the semicircle is necessary for the directions
    II = dirmat < 0;
    dirmat(II) = dirmat(II)+ pi;
    III = dirmat >= pi-phistep;
    dirmat(III) = 0;
    
    % Directional Experimental (Semi-)Variograms
    gexp(nrbins,nrphibins) = 0; %variogram values
    nr_pairs(nrbins,nrphibins) = 0; %number of pairs per lag
    
    for i=1:nrphibins
        dir_d = dirmat(:,:)>=phic(i)- phitol & dirmat(:,:)<=phic(i) + phitol;
        dir_distmat = distmat(:,3);
        dir_distmat(dir_d(:,:)==0) = [];
        dir_gx = gx;
        dir_gx(dir_d(:,:)==0) = [];
        % bin distance
        [~,ixedge] = histc(dir_distmat,edges);
        % variogram anonymous function
        fvar = @(x) (1./(2*numel(x))) * sum(x); %divide with 2 or not??? (because we consider only the upper triangular matrix)
        %calculate variogram
        gexp1 = accumarray(ixedge,dir_gx,[numel(edges) 1],fvar,nan);
        gexp_nrp1 = accumarray(ixedge,ones(size(dir_gx)),[numel(edges) 1],@sum,nan);
        gexp(:,i) = gexp1(1:end-1,:);
        nr_pairs(:,i) = gexp_nrp1(1:end-1,:);
        clear dir_d dir_distmat ixedge dir_gx gexp1 gexp_nrp1
    end
    
    d_30 = nr_pairs <=30;
    gexp(d_30) = NaN;
    % d_NaN = ~isnan(gexp);
    gexpmax = max(max(gexp));
    
    % Cartesian Coordinates of Bin Centers
    phirep = repmat(phic,nrbins,1);
    distrep = repmat(distc,1,nrphibins);
    [rx,ry] = pol2cart(phirep,distrep);
    [rxs,rys] = pol2cart(phirep+pi,distrep); %dimetrically opposite points
    
    c_centers.c{1,1} = rx;
    c_centers.c{1,2} = ry;
    c_centers.distc = distc;
    c_centers.phic = phic;
    
    if flag == 1
        % Plot experimental semivariogram
        cls = 13;
        %cls = 10*disttol;
        X = reshape( rx.' ,1,numel(rx));
        Y = reshape( ry.' ,1,numel(ry));
        XS = reshape( rxs.' ,1,numel(rxs));
        YS = reshape( rys.' ,1,numel(rys));
        gexp2 = reshape( gexp.' ,1,numel(gexp));
        figure;
        h_fake = polar(phirep,maxdist*ones(size(phirep)));
        hold on
        set(h_fake, 'Visible', 'Off');
        scatter(X, Y, cls, gexp2,'filled');
        scatter(XS, YS, cls, gexp2,'filled');
        %title('Directional Experimental (Semi-)Variograms (2D)')
        clb = colorbar;
        clb.Label.String = '(Semi-)Variogram,\gamma(x)';
        
        figure;
        gridDelaunay = delaunay(X,Y);
        trimesh(gridDelaunay,X,Y,gexp2)
        hold on
        gridDelaunay = delaunay(XS,YS);
        trimesh(gridDelaunay,XS,YS,gexp2)
        xlabel X; ylabel Y;
        %title('Directional Experimental (Semi-)Variogram (3D)')
        clb = colorbar;
        clb.Label.String = '(Semi-)Variogram,\gamma(x)';
        
        for i=1:nrphibins
            d_NaN = ~isnan(gexp(:,i));
            figure;
            plot(distc(d_NaN),gexp(d_NaN,i),'LineWidth',2);
            hold on;
            scatter(distc(d_NaN),gexp(d_NaN,i),'fill');
            axis([0 maxdist 0 gexpmax*1.1]);
            %title(sprintf('Directional Experimental (Semi-)Variogram on {\\phi = %g}^{\\circ}',phic(i)*180/pi));
            xlabel('x');
            ylabel('\gamma(x)');
        end
        
    elseif flag == 2
        
        % Plot experimental semivariogram
        cls = 13;
        %cls = 10*disttol;
        X = reshape( rx.' ,1,numel(rx));
        Y = reshape( ry.' ,1,numel(ry));
        XS = reshape( rxs.' ,1,numel(rxs));
        YS = reshape( rys.' ,1,numel(rys));
        gexp2 = reshape( gexp.' ,1,numel(gexp));
        figure;
        h_fake = polar(phirep,maxdist*ones(size(phirep)));
        hold on
        set(h_fake, 'Visible', 'Off');
        scatter(X, Y, cls, gexp2,'filled');
        scatter(XS, YS, cls, gexp2,'filled');
        %title('Directional Experimental (Semi-)Variograms (2D)')
        clb = colorbar;
        clb.Label.String = '(Semi-)Variogram,\gamma(x)';
        
        figure;
        gridDelaunay = delaunay(X,Y);
        trimesh(gridDelaunay,X,Y,gexp2)
        hold on
        gridDelaunay = delaunay(XS,YS);
        trimesh(gridDelaunay,XS,YS,gexp2)
        xlabel X; ylabel Y;
        %title('Directional Experimental (Semi-)Variogram (3D)')
        clb = colorbar;
        clb.Label.String = '(Semi-)Variogram,\gamma(x)';
        
    end
    
else
    
    error('Input iso must be 0 or 1.')
    
end


end