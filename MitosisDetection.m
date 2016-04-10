function MitosisDetection(I,Iprev,L,LNew,LPrev,Kalmans,Params)
Params.K = 4;
%%
propsL = regionprops(L,I,'area','centroid','PixelValues');
propsLNew = regionprops(LNew,I,'area','centroid','PixelValues');
cents = cat(1,propsL.Centroid);
centsNew = cat(1,propsLNew.Centroid);
propsL = propsL(~isnan(cents(:,1)));
cents = cents(~isnan(cents(:,1)),:);
propsLNew = propsLNew(~isnan(centsNew(:,1)));
centsNew = centsNew(~isnan(centsNew(:,1)),:);
%%
glL = arrayfun(@(p) sum(p.PixelValues),propsL);
glLNew = arrayfun(@(p) sum(p.PixelValues),propsLNew);
gl = cat(1,glL,glLNew);
%%
[indNN,distNN] = knnsearch(cat(1,cents,centsNew),centsNew,'k',Params.K + 1);
indNN = indNN(:,2:end);
distNN = distNN(:,2:end);

%%
[X,Y] = meshgrid(1:Params.K,1:size(indNN,1));
nNew = numel(propsLNew);
cost = inf(nNew,Params.K);
linInd = sub2ind(size(cost),Y(:),X(:));
cost(linInd) = abs(gl(indNN(:))-glLNew(Y(:)));
[minCost,minIdx] = min(cost,[],2);
linInd = sub2ind(size(cost),1:size(cost,1),minIdx');
indNN(linInd)

%%
end