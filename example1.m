                       % Synthetic Data - Test   
                                       
clc;clear variables;close all;

%################## Preliminary Analysis ###########################
%===================================================================
% Construct Synthetic Data
[x,y,rf] = randomfield('Mate',[4,3,1.8,20,0.2,2],60,60,0); %random field
%Mx = 1.2 + 0.1*x + 0.01*y; %trend
data = rf; %synthetic data
nx = 60; ny = 60;

[row_all, col_all, v_all1] = find(data); 

% >>>Original Data Plots<<<

figure; %Simple plot plus colorbar
pcolor(data); shading interp %flat
%title('Synthetic data')
colorbar
set(gca,'XTick',[],'YTick',[])

data_tot_vec = data(:); %total image data in vector form


% Random Sampling (33% of original data)
spoints = 0.33*numel(data); %number of sample points 
[~,idx] = datasample(data_tot_vec,spoints,'Replace',false); %sampling
qpoints = data_tot_vec; qpoints(idx) = 0;
qpoints = reshape(qpoints,ny,nx); %"matrix form" of unknown points
sample = data - qpoints; %"matrix form" of known data
[qrow,qcol,qv1] = find(qpoints); %unknown points coordinations and values - VALIDATION SET
[row,col,v1] = find(sample); %known points coordinations and values - TRAINING SET
N = length(v1); %number of known points

% Plot Sample
figure; %Simple plot plus colorbar
pcolor(sample); shading interp %flat
%title('Random Sample')
colorbar
set(gca,'XTick',[],'YTick',[])

%>>>Data Histograms & Statistical Moments<<<

numbins = 15;
% Total image histogram
figure;                          
%subplot(1,2,1)                 
histfit(data_tot_vec,numbins)
alpha(0.5)
%title('Total image histogram')
% Sample data histogram
figure;%subplot(1,2,2)                 
histfit(v1,numbins)
alpha(0.5)
%title('Sample data histogram') 

% Total image data moments
totim_stats = [min(data(:)) max(data(:)) mean(data(:)) median(data(:)) var(data(:)) skewness(data(:)) kurtosis(data(:))];
% "Drill-holes" data moments
sample_stats = [min(v1) max(v1) mean(v1) median(v1) var(v1) skewness(v1) kurtosis(v1)];


%>>>Normality of data checking and Transformation<<<

% Histogram and Normal Probability Plot
figure; %same as the previous figure plus NPP                       
%subplot(1,2,1)                 
histfit(v1,numbins)
alpha(0.5)
%title('Sample data histogram')    
figure;%subplot(1,2,2)                 
nnp = normplot(v1); 
h_ch=get(gcf,'Children');h_str=get(h_ch(1),'Title');set(h_str,'String',''); % remove normplot title
[h_orig,kst_p_orig,ksstat_orig,cv_orig] = kstest(v1(:));

v = v1;
qv = qv1;
v_all = v_all1;

fluc = v(:); Mx = zeros(size(fluc)); qMx = zeros(size(qcol));
qfluc = qv-qMx; 
fluc_all = v_all-zeros(size(v_all)); 


%>>> Statistical Analysis of Residuals <<<

% Residuals/Fluctuations moments
fluc_stats = [min(fluc) max(fluc) mean(fluc) median(fluc) var(fluc) skewness(fluc) kurtosis(fluc)];


%Maximum distance
c = [col,row];
[k,l] = find(triu(true(N)));
d = hypot(c(k,1)-c(l,1),c(k,2)-c(l,2));
% d = triu(pdist2(c,c));
ncpc = 0.2;
maxdist = max(max(d))*ncpc;
clear k l d
%===================================================================



%% ################## Experimantal Variogram #######################
%===================================================================

% >>> Initialize Basic Parameters <<<
c = {col,row}; %known points' coordinations
qc = {qcol,qrow}; %"unknown" points' coordinations 
%ncpc = 0.2;
N = size(col,1);
Nu = size(qcol,1);
%maxdist = hypot(col(1,1)-col(N,1),row(1,1)-row(N,1))*ncpc;
nrbins = 45;
phistep = 15;
phitol = 20;
N_kr_mod = 3;

% >>> Experimental (Semi-)Variogram (anisotropic) <<<
x = col; y = row; rf = fluc; iso = 0; flag = 1;
[~,~,~] = expvar(x,y,rf,iso,ncpc,nrbins,4,phitol,2); %exper. variogr. of high analysis
[gexp, nr_pairs, c_centers] = expvar(x,y,rf,iso,ncpc,nrbins,phistep,phitol,flag);
gexpmax = max(max(gexp));
%===================================================================


%% ################## DirVar0 ######################################
%===================================================================

% >>> Fitting I and Parameters of Anisotropic Correlation Estimation (s2,xi1,xi2,phi,c0,v or eta1) <<<

% Desired Models
models = {'Gexp';'Gaus';'Sphe';'Mate';'Spar'};
n_models = length(models);

% Initial Values and Limits for optimization
b = [gexpmax,maxdist*1/3,0.5,10,gexpmax/100]; % s2,xi1,R,phi & c0
b_lb = [eps,eps,eps,-90,eps]; b_ub = [gexpmax*1.5,maxdist,30,90,gexpmax/5]; %lower and upper limits
bsp = [1000,maxdist*1/3,0.5,10,gexpmax/100]; % eta0, xi1,R,phi & c0
bsp_lb = [eps,eps,eps,-90,eps]; bsp_ub = [inf,maxdist,30,90,gexpmax/5]; %lower and upper limits

% Summary cells
model_par0 = {[b,1.5];b;b;[b,1.5];[bsp,1]}; %initial parameters values
model_par_lb = {[b_lb,eps];b_lb;b_lb;[b_lb,0.3];[bsp_lb,-2+eps]}; %lower bounds
model_par_ub = {[b_ub,2-eps];b_ub;b_ub;[b_ub,3.5];[bsp_ub,inf]}; %upper bounds

clear b b_lb b_ub bsp bsp_lb bsp_ub 

% Estimation of Parameters (s2,xi1,xi2,phi,c0,v or eta1)
bmodel{n_models,1}=[]; fval(n_models,1)=0; tit{n_models,1}=[];
iso = 00;objmod = 'NWEr_m';N_ms = 35;flag = 1;
for i=1:n_models
    model.function = models{i,1};
    model.params0 = model_par0{i,1};
    model.paramslb = model_par_lb{i,1};
    model.paramsub = model_par_ub{i,1};

   [bmodel{i,1},fval(i,1),tit{i,1}] =...
       variogramfit(gexp,nr_pairs,c_centers,iso,model,objmod,N_ms,flag);  
end

% Anisotropy Estimation (#Not Necessary#)
R0(n_models,1) = 0;phi0(n_models,1) = 0;xi10(n_models,1) = 0;xi20(n_models,1) = 0;
for i=1:n_models
    R0(i,1) = bmodel{i,1}(1,3);
    phi0(i,1) = bmodel{i,1}(1,4);
    xi10(i,1) = bmodel{i,1}(1,2);
    xi20(i,1) = bmodel{i,1}(1,2)/R0(i,1);
    if R0(i,1)>1
        R0(i,1) = 1/R0(i,1);bmodel{i,1}(1,3) = R0(i,1);
        xi101 = xi10(i,1);
        xi10(i,1) = xi20(i,1); bmodel{i,1}(1,2) = xi10(i,1);
        xi20(i,1) = xi101;
        if phi0(i,1)>0
            phi0(i,1) = phi0(i,1)-90;
        else
            phi0(i,1) = phi0(i,1)+90;
        end
        bmodel{i,1}(1,4) = phi0(i,1);
    end
end

R = mean(R0); phi = mean(phi0);xi1 = mean(xi10);xi2 = mean(xi20);


% >>> Cross Validation <<<

% Matrices and Cells preallocation
cv_scores{n_models,1}=[]; cv_checks{n_models,1}=[]; 
cv_matr{n_models,1}=[]; cv_Ss(n_models,6)= 0; 

% Inputs definition
iso = 00; d_col = 1; 

% Cross Validation
for i=1:n_models
    model.function = models{i,1};
    model.params = bmodel{i,1};
    model.r_ok = [6,6];
    
    [cv_scores{i,1}, cv_checks{i,1}, cv_matr{i,1}] =...
        crossval(x,y,rf,iso,d_col,model);
    cv_Ss(i,:) = table2array(cv_scores{i,1}(:,2:end));     
end

% Cross Validation Scores
table_h = {'MeanAbsErr','MaxAbsErr','MSE','RMSE','rpearson','rspearman','finalscore'};
% % table_r = {'Gexp';'Gaus';'Sphe';'Mate1/3';'Mate1.0';'Mate1.5';'Mate2.0';'Mate2.5';
% %     'Mate3.0';'Mate3.5';'Spar'};
table_r = {'Gexp';'Gaus';'Sphe';'Mate';'Spar'};

relMSE = cv_Ss(:,3)/min(cv_Ss(:,3)); %relative MSE
FinalScore = ((1./relMSE).^2).*cv_Ss(:,5).*cv_Ss(:,6); %finalscore
cv_Ss = [cv_Ss, FinalScore];
cv_Ssf = array2table(cv_Ss,'VariableNames',table_h,'RowNames',table_r);
clear relMSE FinalScore

%Trend addition and Boxcox Inversion
cv_St(n_models,6)= 0;
for i=1:n_models
    
    cv_matr{i,1}.Z_tr = cv_matr{i,1}.Z;

    cv_matr{i,1}.Z_ibt1 = cv_matr{i,1}.Z_tr;
    
    cv_matr{i,1}.Z_ibt = real(cv_matr{i,1}.Z_ibt1);
    cv_realZ(i,1) = isreal(cv_matr{i,1}.Z_ibt1); %#ok<SAGROW> 
    
    cv_St(i,:) = correlstats(v1,cv_matr{i,1}.Z_ibt);%total cv scores
          
end

% Total Cross Validation Scores
relMSE = cv_St(:,3)/min(cv_St(:,3)); %relative MSE
FinalScore = ((1./relMSE).^2).*cv_St(:,5).*cv_St(:,6); %finalscore
cv_St = [cv_St, FinalScore];
cv_Stf = array2table(cv_St,'VariableNames',table_h,'RowNames',table_r);
clear relMSE FinalScore

%Plots
for i=1:n_models
    
    % Stochastic Component's figures
    cc1 = [col,row,rf(:);qcol,qrow,zeros(Nu,1)];
    cc2 = [col,row,cv_matr{i,1}.Z(:);qcol,qrow,zeros(Nu,1)];
    Z1(max(cc1(:,2)),max(cc1(:,1)))=0; %#ok<SAGROW>
    Z2=Z1;
    for j = 1:size(cc1,1)
        Z1(cc1(j,2),cc1(j,1)) = cc1(j,3);
        Z2(cc2(j,2),cc2(j,1)) = cc2(j,3);
    end
    
    figure;pcolor(Z1);%title('Sample Stochastic Component'); 
    view(2);shading interp; colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;pcolor(Z2);%title(sprintf('Estimation of Sample Stochastic Component \n%s',tit{i,1}));
    view(2);shading interp;colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;scatter(rf(:),cv_matr{i,1}.Z(:),'filled','d');hold on;
    dvec1 = [rf(:);cv_matr{i,1}.Z(:)];
    plot([min(dvec1)-0.5,max(dvec1)+0.5],[min(dvec1)-0.5,max(dvec1)+0.5],'r');
    axis([min(dvec1)-0.5,max(dvec1)+0.5,min(dvec1)-0.5,max(dvec1)+0.5])
    %title('Scatter Plot')
    xlabel ('Observed Data'); ylabel ('Estimations');
    figure; h = histogram(rf(:),16,'EdgeColor',[0 0 1],'FaceAlpha',0.7);
    hold on
    histogram(cv_matr{i,1}.Z(:),'BinEdges',h.BinEdges,'EdgeColor',[0.2 1 0]','FaceAlpha',0.7)
    %title('Histograms of Sample Stochastic Component')
    legend({'Original', 'Estimated'});
    clear h

    
    % Total Data figures
    cc3 = [col,row,v1(:);qcol,qrow,zeros(Nu,1)];
    cc4 = [col,row,cv_matr{i,1}.Z_ibt(:);qcol,qrow,zeros(Nu,1)];
    Z3(max(cc3(:,2)),max(cc3(:,1)))=0; %#ok<SAGROW>
    Z4=Z3;
    for j = 1:size(cc3,1)
        Z3(cc3(j,2),cc3(j,1)) = cc3(j,3);
        Z4(cc4(j,2),cc4(j,1)) = cc4(j,3);
    end
%     figure;pcolor(Z3);view(2);shading interp; %title('Original Sample'); 
%     colorbar;set(gca,'XTick',[],'YTick',[]); 
    figure;pcolor(Z4);view(2);shading interp; %title(sprintf('Estimation of Original Sample \n%s',tit{i,1}));
    colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;scatter(v1(:),cv_matr{i,1}.Z_ibt(:),'filled','d');hold on;
    dvec2 = [v1(:);cv_matr{i,1}.Z_ibt(:)];
    plot([min(dvec2)-0.5,max(dvec2)+0.5],[min(dvec2)-0.5,max(dvec2)+0.5],'r');
    axis([min(dvec2)-0.5,max(dvec2)+0.5,min(dvec2)-0.5,max(dvec2)+0.5])
    %title('Scatter Plot')
    xlabel ('Observed Data')
    ylabel ('Estimations')
    figure;h = histogram(v1(:),16,'EdgeColor',[0 0 1],'FaceAlpha',0.7);
    hold on
    histogram(cv_matr{i,1}.Z_ibt(:),'BinEdges',h.BinEdges,'EdgeColor',[0.2 1 0],'FaceAlpha',0.7)
    %title('Histograms of Sample Data')
    legend({'Original', 'Estimated'});
    clear h
    
    clear Z1 Z2 Z3 Z4 dvec1 dvec2
            
end


% >>> Ordinary Kriging <<<

% Sort models based on cross validation scores 
[~,ind] = sort(table2array(cv_Stf(:,7)),'descend');

% Matrices and Cells preallocation
Z{N_kr_mod,1}=[]; Z_error{N_kr_mod,1}=[];kr_checks{N_kr_mod,1}=[]; 
kr_matr{N_kr_mod,1}=[];
kr_Ss(N_kr_mod,6)= 0; table_r2{N_kr_mod,1} = [];
CI1{N_kr_mod,1}=[];UNC{N_kr_mod,1}=[];

% Inputs definition
xu = qcol; yu = qrow; iso = 00; 

% Ordinary Kriging
for i=1:N_kr_mod
    model.function = models{ind(i),1};
    model.params = bmodel{ind(i),1};
    model.r_ok = [6,6];
    
    [Z{i,1},Z_error{i,1},kr_checks{i,1}, kr_matr{i,1}] =...
        ordkrig(x,y,rf,xu,yu,iso,model);
    kr_Ss(i,:) = correlstats(Z{i,1},qfluc); 
    table_r2{i,1} = table_r{ind(i),1};
    
    %Confidence Intervals (95%)
    CI1{i,1}.low = Z{i,1} - 1.96*sqrt(Z_error{i,1});
    CI1{i,1}.up = Z{i,1} + 1.96*sqrt(Z_error{i,1});
    CI1{i,1}.uncer = 1.96*sqrt(Z_error{i,1});
    UNC{i,1} = real(CI1{i,1}.uncer);
    realCI(i,1) = isreal(CI1{i,1}.uncer); %#ok<SAGROW>
end

% Kriging Scores
relMSE = kr_Ss(:,3)/min(kr_Ss(:,3)); %relative MSE
FinalScore = ((1./relMSE).^2).*kr_Ss(:,5).*kr_Ss(:,6); %finalscore
kr_Ss = [kr_Ss, FinalScore];
kr_Ssf = array2table(kr_Ss,'VariableNames',table_h,'RowNames',table_r2);
clear relMSE FinalScore

%Trend addition and Boxcox Inversion
kr_St(N_kr_mod,6)= 0; Z_tr{N_kr_mod,1} = 0;
Z_ibt1{N_kr_mod,1} = 0;Z_ibt{N_kr_mod,1} = 0;
for i=1:N_kr_mod
    
    Z_tr{i,1} = Z{i,1};

    Z_ibt1{i,1} = Z_tr{i,1};
    
    Z_ibt{i,1} = real(Z_ibt1{i,1});
    kr_realZ(i,1) = isreal(Z_ibt1{i,1});  %#ok<SAGROW>
    
    kr_St(i,:) = correlstats(qv,Z_ibt{i,1});%total kriging scores
          
end

% Total Kriging Scores
relMSE = kr_St(:,3)/min(kr_St(:,3)); %relative MSE
FinalScore = ((1./relMSE).^2).*kr_St(:,5).*kr_St(:,6); %finalscore
kr_St = [kr_St, FinalScore];
kr_Stf = array2table(kr_St,'VariableNames',table_h,'RowNames',table_r2);
clear relMSE FinalScore

%Plots
for i=1:N_kr_mod
    
    % Stochastic Component's figures
    cc1 = [col,row,rf(:);qcol,qrow,qfluc(:)];
    cc2 = [col,row,rf(:);qcol,qrow,Z{i,1}];
    Z1(max(cc1(:,2)),max(cc1(:,1)))=0; %#ok<SAGROW>
    Z2=Z1;
    for j = 1:size(cc1,1)
        Z1(cc1(j,2),cc1(j,1)) = cc1(j,3);
        Z2(cc2(j,2),cc2(j,1)) = cc2(j,3);
    end
    figure;pcolor(Z1);%title('Original Stochastic Component'); 
    view(2);shading interp; colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;pcolor(Z2); %title(sprintf('Estimation of Stochastic Component \n%s',tit{ind(i),1}));
    view(2);shading interp; colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;scatter(qfluc(:),Z{i,1}(:),'filled','d');hold on;
    dvec1 = [qfluc(:);Z{i,1}(:)];
    plot([min(dvec1)-0.5,max(dvec1)+0.5],[min(dvec1)-0.5,max(dvec1)+0.5],'r');
    axis([min(dvec1)-0.5,max(dvec1)+0.5,min(dvec1)-0.5,max(dvec1)+0.5])
    %title('Scatter Plot')
    xlabel ('Observed Data')
    ylabel ('Estimations')
    figure;h = histogram(qfluc(:),16,'FaceColor',[0 0 1],'FaceAlpha',0.7);
    hold on
    histogram(Z{i,1}(:),'BinEdges',h.BinEdges,'FaceColor',[0.2 1 0],'FaceAlpha',0.7)
    %title('Histograms of Stochastic Component')
    legend({'Original', 'Estimated'});
    clear h
    
    % Total Data figures
    cc3 = [col,row,v1(:);qcol,qrow,Z_ibt{i,1}];
    cc4 = [col,row,zeros(size(col,1),1); qcol,qrow,UNC{i,1}];
    Z3(max(cc3(:,2)),max(cc3(:,1)))=0; %#ok<SAGROW>
    Z4=Z3; 
    for j = 1:size(cc3,1)
        Z3(cc3(j,2),cc3(j,1)) = cc3(j,3);
        Z4(cc4(j,2),cc4(j,1)) = cc4(j,3);
    end
    figure;pcolor(Z3); view(2);shading interp; %title(sprintf('Estimation of Original Data \n%s',tit{ind(i),1}));
    colorbar;set(gca,'XTick',[],'YTick',[]);
    figure;scatter(qv1(:),Z_ibt{i,1}(:),'filled','d');hold on;
    dvec2 = [qv1(:);Z_ibt{i,1}(:)];
    plot([min(dvec2)-0.5,max(dvec2)+0.5],[min(dvec2)-0.5,max(dvec2)+0.5],'r');
    axis([min(dvec2)-0.5,max(dvec2)+0.5,min(dvec2)-0.5,max(dvec2)+0.5])
    %title('Scatter Plot')
    xlabel ('Observed Data')
    ylabel ('Estimations')
    figure;h = histogram(qv1(:),16,'FaceColor',[0 0 1],'FaceAlpha',0.7);
    hold on
    histogram(Z_ibt{i,1}(:),'BinEdges',h.BinEdges,'FaceColor',[0.2 1 0],'FaceAlpha',0.7)
    %title('Histograms of Data')
    legend({'Original', 'Estimated'});
    figure;pcolor(Z4); view(2);shading interp; %title('95% Confidence Interval'); 
    colorbar;set(gca,'XTick',[],'YTick',[]);
    clear h
    
    clear cc1 cc2 cc3 cc4 Z1 Z2 Z3 Z4 dvec1 dvec2
            
end
%===================================================================


