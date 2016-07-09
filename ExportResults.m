function ExportResults(dataPath,dataExt,dataExpCFP,dataExpYFP,resPath,resExt,resExpr,Name)

DataCFP = Load_Data(dataPath,dataExt,dataExpCFP);
DataYFP = Load_Data(dataPath,dataExt,dataExpYFP);
Res  = Load_Data(resPath,resExt,resExpr);
S = imread(Res.Frame_name{end});
maxID = max(S(:));
Area = nan(maxID,Res.Frame_Num);
CFP = nan(maxID,Res.Frame_Num);
YFP = nan(maxID,Res.Frame_Num);
X = nan(maxID,Res.Frame_Num);
Y = nan(maxID,Res.Frame_Num);

for t = 1:Res.Frame_Num
    S = imread(Res.Frame_name{t});
    uniqueS = unique(S(S(:)>0));
    ICFP = imread(DataCFP.Frame_name{t});
    IYFP= imread(DataYFP.Frame_name{t});
    regCFP = regionprops(S,ICFP,'Area','MeanIntensity','Centroid');
    regYFP = regionprops(S,IYFP,'MeanIntensity');
    cents =  cat(1,regCFP(uniqueS).Centroid);
    Area(uniqueS,t) = [regCFP(uniqueS).Area];
    CFP(uniqueS,t) = [regCFP(uniqueS).MeanIntensity];
    YFP(uniqueS,t) = [regYFP(uniqueS).MeanIntensity];
    X(uniqueS,t) = cents(:,1);
    Y(uniqueS,t) = cents(:,2);
end
save(fullfile(resPath,'..','Results.mat'),'Area','CFP','YFP','X','Y')
