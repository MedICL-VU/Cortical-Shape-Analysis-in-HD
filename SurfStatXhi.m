function [ pval, peak, clus, clusid ] = SurfStatXhi( slm, mask, clusthresh )


[l,v]=size(slm.t);
if nargin<2 | isempty(mask)
    mask=logical(ones(1,v));
end
if nargin<3
    clusthresh=0.001;
end

if clusthresh<1
    thresh=stat_threshold(0,1,0,df,clusthresh,[],[],[],slm.k,[],[],0);
else
    thresh=clusthresh;
end

[resels,reselspvert,edg]=SurfStatResels(slm,mask);


if max(slm.t(1,mask))<thresh
    
    pval.P=stat_threshold(resels,v,1,df,[10 slm.t],[],[],[],slm.k,[],[],0);
    pval.P=pval.P((1:v)+1);
    peak=[];
    clus=[];
    clusid=[];
    
else
    % Peak and Cluster
    [peak,clus,clusid]=SurfStatPeakClus(slm,mask,thresh,reselspvert,edg);
    
    % Peak
    [pp,clpval]=stat_threshold(resels,v,1,df,...
        [10 peak.t' slm.t(1,:)],thresh,[10; clus.resels],[],slm.k,[],[],0);
    peak.P=pp((1:length(peak.t))+1)';
    pval.P=pp(length(peak.t)+(1:v)+1);
    
    % Cluster
    if slm.k>1
        j=(slm.k-1):-2:0;
        sphere=zeros(1,slm.k);
        sphere(j+1)=exp((j+1)*log(2)+(j/2)*log(pi)+gammaln((slm.k+1)/2)- ...
            gammaln(j+1)-gammaln((slm.k+1-j)/2));
        sphere=sphere.*(4*log(2)).^(-(0:(slm.k-1))/2)/ndf;
        [pp,clpval]=stat_threshold(conv(resels,sphere),Inf,1,df,...
            [],thresh,[10; clus.resels],[],[],[],[],0);
    end
    clus.P=clpval(2:length(clpval));
    pval.C=interp1([0; clus.clusid],[1; clus.P],clusid);        
end

tlim=stat_threshold(resels,v,1,df,[0.5 1],[],[],[],slm.k,[],[],0);
tlim=tlim(2);
pval.P=pval.P.*(slm.t(1,:)>tlim)+(slm.t(1,:)<=tlim);
pval.mask=mask;







end


