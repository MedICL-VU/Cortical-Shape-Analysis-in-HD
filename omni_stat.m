close all
clear all

%% load data

addpath(genpath('Tools/surfstat'));    % path to surfstat lib

% Measure of interest ('lgi', 'sd', 'ct')
cmeasure = 'lgi';

load(sprintf('data/y_%s', cmeasure));
load('data/demographics_ticv.mat');  % demographics
load('env/environment.mat');         % template surfaces, parcellation boundary, etc.

%% Model initialization
valid = find(demographics.cap_e_group==0 | demographics.cap_e_group==1 | demographics.cap_e_group==2 | demographics.cap_e_group==3);
demographics = demographics(valid,:);
sex = cell(length(demographics.gender),1);
sex(demographics.gender == 1) = {'M'};
sex(demographics.gender == 0) = {'F'};

dx = cell(length(demographics.cap_e_group),1);
dx(demographics.cap_e_group==0) = {'CNTRL'};
dx(demographics.cap_e_group==1) = {'LOW'};
dx(demographics.cap_e_group==2) = {'MED'};
dx(demographics.cap_e_group==3) = {'HIGH'};

subj = find((demographics.cag > 0 & demographics.cag < 37) | demographics.cag >= 37);
Y0=Y0(subj,:);

age = demographics.age_at_scan(subj);
sex = sex(subj);
dx = dx(subj);
subject = demographics.external_id(subj);
entry = zeros(length(subject),1);
for i = 1:length(subject)
    entry(i) = min(age(strcmp(subject,subject(i))));
end
duration = age - entry;
scanner = demographics.scanner(subj);

Duration = term( duration );
Scanner = term( scanner );
Subject = term( subject );
Age = term( age );
Sex = term( sex );
Pathology = term( dx );
Cntrl = term ( Pathology.CNTRL );
CapL = term( Pathology.LOW );
CapM = term( Pathology.MED );
CapH = term( Pathology.HIGH );

%% model setup
contrast = Pathology.CNTRL - Pathology.MED;

%For CAP: 
M=Age+Duration+Pathology+Sex+random(Subject)+I;
fwhm = 1;
if strcmp(cmeasure, 'ct')
    fwhm = 6;
    Y0 = SurfStatSmooth(Y0, surfwhite, fwhm);
end

%% fitting
Y0(:,sum(abs(Y0))==0) = rand(size(Y0(:,sum(abs(Y0))==0)))*eps;   % prevent numerical unstability
slm = SurfStatLinMod(Y0,M,surfwhite);                            % fitting

%% correction
load(sprintf('data/chi^2_%s',cmeasure));
alpha = 0.01;
[pval, peak, clus] = rft_fwer(chis, 3, alpha, fwhm, mask_b, slm);

%% Visualization: signficant regions after RFT 
pval.mask = mask_b;
figure; SurfStatView(pval, surfinfl, 'p-val after multi-correction' );

%pval.mask(:) = true;
%figure; SurfStatView2( pval, surfinfl, 'p-val after multi-correction' );

