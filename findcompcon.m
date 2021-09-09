function BBW = findcompcon(bw, topo)

BBW = [];
BBW.VoxelIdxList = {};
BBW.NelemPerObject = [];

idx = find(bw>0);

topor = zeros(size(topo));
for i=1:length(idx)
    topor = topor + double(topo==idx(i));
end
topobin = sum(topor,2)>0;
n = sum(topobin)>0;

idtopo = find(topobin>0);
toponew = topo(idtopo,:);


n=0;
while ~isempty(idx)
    n = n+1;
    v = [idx(1)];
    idx = idx(2:end);
    w = v;
    
    lfv = 0;
    while length(v)>lfv
        lfv = length(v);
        
        u = [];
        for i=1:length(w)
            topobin = sum(toponew==w(i),2)>0;
            idtopo = find(topobin>0);
            for j=1:length(idtopo)
                for m=1:3
                    if sum(u==toponew(idtopo(j),m))==0 && sum(v==toponew(idtopo(j),m))==0
                        u = [u toponew(idtopo(j),m)];
                    end
                end
            end
        end
        
        idxnew = [];
        for i=1:length(idx)
            b=0;
            for j=1:length(u)
                if u(j)==idx(i)
                    b=b+1;
                end
            end
            if b==0
                idxnew = [idxnew idx(i)];
            end
        end
        idx = idxnew;
        v = [v u];
        w = u;
    end
    
    
    BBW.VoxelIdxList{n} = v;
    BBW.NelemPerObject = [BBW.NelemPerObject length(v)];
end
BBW.NumObjects = n;

end