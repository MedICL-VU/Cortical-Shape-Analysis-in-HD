
%% setup
addpath(genpath('Tools/surfstat'));    % path to surfstat lib
cmeasure = 'CT';   % measure of interest

%% read data
% % load(sprintf('data/%s.mat',cmeasure));  % surface measure
% total_samples = size(demographics, 1);
% lh_data = cell(total_samples, 1);
% rh_data = cell(total_samples, 1);
% for i = 1:total_samples
% %     lh_data{i} = sprintf('%s/%s/%s-x-Reg/lh.mid.reg.%s.txt', data_p, demographics.project_id{i}, demographics.session_label{i}, cmeasure);
%     lh_data{i} = sprintf('%s/lh.mid.reg.%s.txt', demographics.fpath{i}, cmeasure);
%     rh_data{i} = sprintf('%s/rh.mid.reg.%s.txt', demographics.fpath{i}, cmeasure);
% end
% Y0=SurfStatReadData([lh_data, rh_data]);

load('data/y_ct');
load('data/demographics.mat');  % demographics
load('env/environment.mat');  % template surfaces, parcellation boundary, etc.



%% data filtering
subj = find(demographics.age_at_scan > 0 & demographics.cap_d_group>-1); % dx

Y=Y0(subj,:);

%% model initialization
% dx = cell(length(demographics.DEPRESS),1);
sex = cell(length(demographics.gender),1);
sex(demographics.gender == 1) = {'M'};
sex(demographics.gender == 0) = {'F'};

cag = demographics.cag;
% dx = demographics.diagnosis;
% project = demographics.Project;
cag = cag(subj);
age = demographics.age_at_scan(subj);
sex = sex(subj);
capd = demographics.cap_d_group(subj);
% dx  = dx(subj);
dx = cell(length(cag),1);
dx(demographics.cag <= 35) = {'HC'};
dx(demographics.cag > 35)  = {'HD'};



capg = cell(length(cag),1);
capg(capd == 0) = {'HC'};
capg(capd == 1) = {'LOW'};
capg(capd == 2) = {'MED'};
capg(capd == 3) = {'HIGH'};

scanner = demographics.scanid;

subject = demographics.external_id(subj);
% project = project(subj);

% last ='';
% nsubj = [0,0];
% for i=1:length(cag)
%     curr = subject{i};
%     
%     if ~strcmp(curr,last)
%         nsubj((cag(i)>35)+1) = nsubj((cag(i)>35)+1) + 1;
%     end
%     
%     last = curr;
% end


Age = term( age );
Sex = term( sex );
Patho = term(capg);
Subject = term(subject);
Scanner = term(scanner);
% Project = term( project );

%% model setup/fitting1
contrast = Patho.HC - Patho.HIGH;

M=Age+Sex+Patho+random(Subject)+I;

% contrast = age.*Depress.HC - age.*Depress.DP;
% M=1+Age+Sex+Depress+Age*Depress;

% Y = SurfStatSmooth(Y,surfwhite,10);
Y(:,sum(abs(Y))==0) = rand(size(Y(:,sum(abs(Y))==0)))*eps;   % prevent numerical unstability


slm = SurfStatLinMod(Y, M, surfwhite);               % fitting


% Omnibus test

slm = SurfStatT( slm, contrast );                   % t-stat


[ pval, peak, clus ] = SurfStatP( slm, mask, 0.01);  % cluster raw-p



%% residual: adjusted Y
X=Patho;
% X=Age*Depress;
G=term(1);
XG=X+G+X*G;
Yhat=slm.X*slm.coef;
XGmat=double(XG);
Yhatadj=XGmat*pinv(XGmat)*Yhat;
c=mean(Y)-mean(Yhatadj);
Yhatadj=Yhatadj+c;
Yadj=Yhatadj+Y-Yhat;
Yadj(:,~mask) = 0;

% %% visualization: mean
% meansubj = mean( double( Yadj ) );
% figure; SurfStatView( meansubj, surfinfl, 'adjusted mean' ); SurfStatColLim( [1 15] );   % 1-5: ct, 0-30: sd, 1-15: lgi
% 
% %% visualizaion: difference
% diff = mean(Yadj(contains(dx,'HC'),:))-mean(Yadj(contains(dx,'HD'),:));
% bdry = min(diff)-1;
% diff(~mask_b)=bdry;
% figure; SurfStatView( diff, surfinfl, 'adjusted difference (HC - HD)' );
% cmap = [0 0 0; jet];
% SurfStatColormap( jet );
% SurfStatColLim( [-.5 .5] );

EF =  slm.ef./slm.sd;
figure; SurfStatView( EF, surfinfl, 'Effect size' );

%% visualization: signficant regions after RFT
pval.mask = mask_b;
figure; SurfStatView( pval, surfinfl, 'p-val after multi-correction' );
% figure; SurfStatView( pval, surfwhite, 'p-val after multi-correction' );  % white surface

%% plot: box plot for significant regions
clusid=1;   % change this to see different clusters (also, use term( clus ))
cpval = clus.P(clusid);

cc = [1 0 0; 0 0 1];
figure('Position', [10 10 300 700]);
boxplot( mean(Yadj( :, pval.C == cpval ),2), dx, 'Notch','off', 'Colors', cc );
ylabel(sprintf('%s',cmeasure));
xlabel(sprintf('cluster p: %f',cpval));

h = findobj(gca,'tag','Median');
for j=1:length(h)
    set(h(j),'linestyle','-.','Color',[0 0 0],{'linew'},{2});
end
set(findobj(gca,'tag','Outliers'),'MarkerEdgeColor',[0 0 0]);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    c = cc(3-j,:);
    patch(get(h(j),'XData'),get(h(j),'YData'), c,'FaceAlpha',.51,'EdgeColor',cc(3-j,:));
end

%% plot: group by age interaction
clusid=1;

cpval = clus.P(clusid);
figure; SurfStatPlot( age, double(mean(Y(:, pval.C == cpval),2)), 1, dx ); ylabel(sprintf('%s',cmeasure));