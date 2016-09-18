function [varargout]=CheckJoseIlastik(txtfilePath,data_path,extention,expr,resPath,Name)
if ~exist('Name','var')||isempty(Name)||~ischar(Name)
    Name = '';
end
%%
manualTrack = readtable(txtfilePath,'Delimiter','\t','Format','%f%f%f%f%f%f%f');
manualTrackOrig = manualTrack;
varargout{1} = manualTrackOrig;
%%
%load(fullfile(resPath,'..','Link.mat'));
%%
Data = Load_Data(data_path,extention,expr) ;
fileinfo = hdf5info(resPath);
ResData = hdf5read(fileinfo.GroupHierarchy.Datasets(1));
ResData = squeeze(reshape(ResData,fileinfo.GroupHierarchy.Datasets.Dims));
%%
manualTrack = readtable(txtfilePath,'Delimiter','\t','Format','%f%f%f%f%f%f%f');
b = 10;
H = size(ResData,1);
W = size(ResData,2);
FrameNum = size(ResData,3);
%%manualTrack = readtable(txtfilePath,'Delimiter','\t','Format','%f%f%f%f%f%f%f');
manualTrackOrig = manualTrack;
existingCells = manualTrack.cell_id(manualTrack.timepoint==1);
existingSisters = [];
LinkMat = [];
start = [];
stop = [];
maxID = max(manualTrack.cell_id);
LinkStruct = struct([]);
%%
for t = 2:max(manualTrack.timepoint)
    manualTrack= manualTrack(manualTrack.cell_id>0,:);
    centroids = table2array(manualTrack(manualTrack.timepoint==t,2:3));
    currentCells =  manualTrack.cell_id(manualTrack.timepoint==t);
    futureCells =  manualTrack.cell_id(manualTrack.timepoint>t);
    [~,ia,ic] = unique(centroids,'rows');
    h = hist(ic,1:numel(ia));
    sisters = [];
    div = manualTrack.division(manualTrack.timepoint==t);
    for i = find(h>1)
        sisPair = currentCells(ic==i)';
        old = ismember(sisPair,existingCells);
        if all(old)
            continue;
        end
        currentCells(ismember(currentCells,sisPair(~old)))=0;
        sisters = cat(1,sisters,sisPair);
        div(ic==i) = 0;
    end
    sisters = cat(1,sisters,[currentCells(logical(div)),-currentCells(logical(div))]);
    if ~isempty(existingSisters)&&~isempty(sisters)

    splitSis = setdiff(existingSisters,sisters,'rows');
    elseif ~isempty(existingSisters)&&isempty(sisters)
        splitSis = existingSisters;
    else
        splitSis = [];
    end
    for i = 1:size(splitSis,1)

        splitSisi = splitSis(i,:);
        if splitSisi(2)<0
            motherCell = splitSisi(1);
            daughterCell = splitSisi(2);
        else
        existSis = ismember(splitSisi,existingCells);
        motherCell = splitSisi(existSis);
        daughterCell = splitSisi(~existSis);
        end
        maxID = maxID +1;
        futureCells(futureCells==motherCell) = maxID;
        currentCells(currentCells==motherCell) = maxID;
        motherCentroid = table2array(manualTrack(manualTrack.timepoint==(t-1)&manualTrack.cell_id==motherCell,2:3));
        daughter1Centroid = centroids(currentCells==maxID,:);

        if isempty(daughter1Centroid)
           maxID = maxID-1;
           stop = cat(1,stop,motherCell);
           continue
        end
        LinkStruct(end+1).t = t;
        LinkStruct(end).Mother.ID = motherCell;
        LinkStruct(end).Mother.Centroid = motherCentroid;
        LinkStruct(end).Daughter(1).ID = maxID;
        LinkStruct(end).Daughter(1).Centroid = daughter1Centroid;

        if daughterCell<0
            daughter2Centroid = [-1,-1];
        else
             daughter2Centroid = centroids(currentCells==daughterCell,:);
             if isempty(daughter2Centroid)
                 daughter2Centroid = [-1, -1];
             end
             LinkStruct(end).Daughter(2).ID = daughterCell;
             LinkStruct(end).Daughter(2).Centroid = daughter1Centroid;
        end

        LinkMat = cat(1,LinkMat,[t,motherCell,motherCentroid,maxID,daughter1Centroid,daughterCell,daughter2Centroid]); %
        start = cat(1,start,maxID,daughterCell(daughterCell>0));
        stop = cat(1,stop,motherCell);
    end

    for edgeIds =  manualTrack.cell_id(manualTrack.timepoint==t&...
            (manualTrack.centroid_col<b|manualTrack.centroid_row<b|manualTrack.centroid_col>(-b+H)...
            |manualTrack.centroid_row>(-b+W))&manualTrack.cell_id>0)'
            futureCells(futureCells==edgeIds) = 0;

    end
    existingCells = union(existingCells,currentCells(currentCells>0));
    manualTrack.cell_id(manualTrack.timepoint==t) = currentCells;
    manualTrack.cell_id(manualTrack.timepoint>t) = futureCells;
    if ~isempty(existingSisters)&&~isempty(sisters)
        existingSisters = union(existingSisters,sisters,'rows');
    elseif isempty(existingSisters)&&~isempty(sisters)
        existingSisters = sisters;
    end
    if ~isempty(existingSisters)&&~isempty(splitSis)
        existingSisters = setdiff(existingSisters,splitSis,'rows');
    end


end
%%
manualTrack = manualTrack(manualTrack.cell_id>0,:);
%%
varargout{1} = manualTrack;
varargout{2} = manualTrackOrig;
%%
Match = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),FrameNum));
GTgl = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),Data.Frame_Num));
GTvar = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),Data.Frame_Num));
AutoGL = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),FrameNum));
AutoVar = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),FrameNum));
%%
S = ResData(:,:,end);
maxIDAuto = max(S(:));
AutoGL2 = nan(maxIDAuto,min(max(manualTrack.timepoint),FrameNum));
AutoVar2 = nan(max(manualTrack.cell_id),min(max(manualTrack.timepoint),FrameNum));

for t = 1:Data.Frame_Num
    %%
     I = double(imread(Data.Frame_name{t}));
    if t<=FrameNum
    S = ResData(:,:,t);
     %for s = unique(S(S>0))'
     %   AutoGL2(s,t) = mean(I(S(:)==s));
     %
     %   AutoVar2(s,t) = var(I(S(:)==s));
     %end
    end

    sizeS = size(I);

    for i = find(manualTrack.timepoint==t)'
        col = manualTrack.centroid_col(i) + [-1,-1,-1,0,0,0,1,1,1]';
        row = manualTrack.centroid_row(i) + [-1,0,1,-1,0,1,-1,0,1]';
        col =min( max(col,1),sizeS(2));
        row =min( max(row,1),sizeS(1));
        cents = sub2ind(sizeS,col,row);
        if t<=FrameNum
        Match(manualTrack.cell_id(i),t)=median(S(cents));
        %AutoGL(manualTrack.cell_id(i),t) = mean(I(S(:)== Match(manualTrack.cell_id(i),t)));
        %AutoVar(manualTrack.cell_id(i),t)= var(I(S(:)== Match(manualTrack.cell_id(i),t)));
        end
        GTgl(manualTrack.cell_id(i),t)=mean(I(cents));
        GTvar(manualTrack.cell_id(i),t)=var(I(cents));
    end

    %cents = sub2ind(size(S),manualTrack.centroid_col(manualTrack.timepoint==t),manualTrack.centroid_row(manualTrack.timepoint==t));
    %Match(manualTrack.cell_id(manualTrack.timepoint==t),t)=S(cents);
    %%
end
MatchOrig = Match;

%%
%{
for l = numel(Link):-1:1
    for m = 1:numel(Link(l).Children)
      Match(Match== Link(l).Children(m)) = Link(l).Mother;
    end
end
%}
%%
varargout{3} = Match;
varargout{4} = MatchOrig;

%%
med_filt = vision.MedianFilter([1,7]);
MatchMed = step(med_filt,Match);
MatchMed(:,end-2:end) = Match(:,end-2:end);
%%
varargout{5} = MatchMed;
%%
uMatch = cell(size(Match,1),1);
for t =1:size(Match,2)
for m = 1:size(Match,1);
    rowM = MatchMed(m,1:t);
    uMatch{m} = unique(rowM(~isnan(rowM)));
    L(m,t) = numel(uMatch{m});

end
acc(t) = sum(L(:,t)==1&(sum(~isnan(Match),2)>=10))./sum(L(:,t)>0&(sum(~isnan(Match),2)>=10));
if t>1
    diffM = diff(MatchMed(:,1:t),1,2);
    acc2(t) = 1-sum(sum(diffM~=0&~isnan(diffM)))./sum(sum(~isnan(MatchMed(:,1:t-1))));
else
    acc2(t) = 1;
end
end
varargout{6} = diffM;
varargout{7} = L;
varargout{8} = GTgl;
varargout{9} = GTvar;
%varargout{10} = AutoGL;
%varargout{11} = AutoVar;
%varargout{12} = AutoGL2;
%varargout{12} = AutoVar2;
%%

figure; plot(acc);
xlabel('t'); ylabel('accuracy'); title(sprintf('%s: Accuracy of Full Track Measures',Name));
%%
figure; plot(acc2);
xlabel('t'); ylabel('accuracy'); title(sprintf('%s: Accuracy of Frame Pair Measure',Name))
%%
pathlength = hist(manualTrack.cell_id,[1:max(manualTrack.cell_id)]);

figure;
[pathLengthHist,histBins,histBinLoc]=histcounts(pathlength(pathlength>0),[1:5:(max(pathlength)+5)]);

Lvalid = L(pathlength>0,end);
for hb = unique(histBinLoc)

    trackErrByLength(hb) = sum(Lvalid(histBinLoc==hb)>1);
end

bar(histBins(1:end-1)+diff(histBins)/2-1,pathLengthHist,1); hold on;
bar(histBins(1:end-1)+diff(histBins)/2-1,trackErrByLength,0.8,'r'); hold off;
legend('Manual Track','Error')
xlabel('path lengths'); ylabel('freq'); title(sprintf('%s: Path Lebgths of Manual Tracks & Errors',Name))
%%
errPerFrame = sum(diffM~=0&~isnan(diffM),1);
cellNum = sum(~isnan(diffM),1);
figure; bar(cellNum); hold on; bar(errPerFrame,'r'); hold off;
legend('Manual Track','Error')
xlabel('frame #'); ylabel('freq'); title(sprintf('%s:Pairwise Manual Tracks & Errors: %d erros out of %d',Name,sum(errPerFrame),sum(cellNum)))


%%
%{
thr = 300;
for i = find(pathlength>thr)
    if ~all(isnan(GTgl(i,:)))
    figure;
    plot(GTgl(i,:)); hold on;
    plot(AutoGL(i,:));
    hold off;
    xlabel('frame #'); ylabel('gl'); title(sprintf('%s:Gray level of cell ID #%d\n Track Length:%d',Name,i,pathlength(i)))
    legend('Centroid GL','Mean GL');

    figure;
    subplot(2,1,1);
    plot(AutoGL(i,:));
    xlabel('frame #'); ylabel('mean gl'); title(sprintf('%s: Auto Algo - Gray level & Variance of cell ID #%d\n Track Length:%d',Name,i,pathlength(i)))
    subplot(2,1,2);
    plot(AutoVar(i,:));
    xlabel('frame #'); ylabel('var gl');
    end
end
%}
%%
totalLinks = 0;
correctLinks = nan(2,numel(LinkStruct));
correctLinksAll = zeros(size(MatchMed,2),2);
sisLinks = 0;
for i = 1:numel(LinkStruct)
     if LinkStruct(i).t>size(MatchMed,2)
        break;
    end
    %%

    time = LinkStruct(i).t;

    correctLinks(:,i)=-1;
    mother = LinkStruct(i).Mother.ID;
    timeVec = max(1,time-5):time-1;
    motherAuto = median(MatchMed(mother,timeVec));
    timeVec = time:min(size(MatchMed,2),time+4);
    daughterAuto = median(MatchMed([LinkStruct(i).Daughter.ID],timeVec),2);
    for d = 1:numel(daughterAuto)
      totalLinks = totalLinks+1;
      if daughterAuto(d)==motherAuto
          correctLinks(d,i) =1;
      else
           correctLinks(d,i) =0;
      end
    end
    if d==2&daughterAuto(1)==daughterAuto(2)
        sisLinks = sisLinks+1;
    end
    correctLinksAll(time,1) = correctLinksAll(time,1)+sum(correctLinks(:,i)>0);
    correctLinksAll(time,2) = correctLinksAll(time,2)+sum(correctLinks(:,i)>=0);
end

    fprintf('Sister Link Accuracy: %0.2f\n',sum(correctLinks(:)>0)./sum(correctLinks(:)>=0)*100)
figure; bar(correctLinksAll(:,2)); hold on; bar(correctLinksAll(:,2)-correctLinksAll(:,1),'r'); hold off;
legend('Manual Link','Error')
xlabel('frame #'); ylabel('freq'); title(sprintf('%s:Mitosis Links & Errors: %d erros out of %d',Name,sum(correctLinksAll(:,2))-sum(correctLinksAll(:,1)),sum(correctLinksAll(:,2))))

%%
%%
return
%%
visPath = fullfile(resPath,'..','Visualize');
visExp = strrep(resExp,'Seg','Vis');
visData = Load_Data(visPath,resExt,visExp);
%%
[cc,tt] = find(diffM~=0&~isnan(diffM))
ind = 0
%%
 ind = ind+1;
 cellid = cc(ind);
 tp = tt(ind);


 MatchMed(cellid,tp-2:tp+2)
 mt = manualTrack(manualTrack.timepoint==tp&manualTrack.cell_id==cellid,1:4)
 I1 = imread(Data.Frame_name{tp-2}); I2 = imread(Data.Frame_name{tp-1});
 I3 = imread(Data.Frame_name{tp});
 I4 = imread(Data.Frame_name{tp+1}); I5 = imread(Data.Frame_name{tp+2});


 p = prctile(double(I1(:)),[1,99]);

 figure(100);

 %imshow(I1,p); a(1) = gca; title(sprintf('Frame #%d',tp-2))
 %figure(101);
 subplot(1,12,[1:4]);
 imshow(I2,p); a(1) = gca; title(sprintf('Frame #%d',tp-1));hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp-1&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp-1&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 %figure(103);
 subplot(1,12,[9:12]);
 imshow(I4,p); a(2) = gca; title(sprintf('Frame #%d',tp+1));hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp+1&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp+1&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 %figure(104);
 %imshow(I5,p); a(5) = gca; title(sprintf('Frame #%d',tp+2))
 %figure(102);
 subplot(1,12,[5:8]);
 imshow(I3,p); a(3) = gca;hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 title(sprintf('Frame #%d',tp))

 I1 = imread(visData.Frame_name{tp-2}); I2 = imread(visData.Frame_name{tp-1});
 I3 = imread(visData.Frame_name{tp});
 I4 = imread(visData.Frame_name{tp+1}); I5 = imread(visData.Frame_name{tp+2});


 p = [];%prctile(double(I1(:)),[1,99]);

 figure(101);

 %imshow(I1,p); a(1) = gca; title(sprintf('Frame #%d',tp-2))
 %figure(101);
 subplot(1,12,[1:4]);
 imshow(I2,p); a(4) = gca; title(sprintf('Frame #%d',tp-1)) ; title(sprintf('Frame #%d',tp-1));hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp-1&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp-1&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 %figure(103);
 subplot(1,12,[9:12]);
 imshow(I4,p); a(5) = gca; title(sprintf('Frame #%d',tp+1)); title(sprintf('Frame #%d',tp+1));hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp+1&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp+1&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 %figure(104);
 %imshow(I5,p); a(5) = gca; title(sprintf('Frame #%d',tp+2))
 %figure(102);
 subplot(1,12,[5:8]);
 imshow(I3,p); a(6) = gca;hold on;scatter(manualTrack.centroid_row(manualTrack.timepoint==tp&manualTrack.cell_id==cellid),manualTrack.centroid_col(manualTrack.timepoint==tp&manualTrack.cell_id==cellid),'*g','linewidth',2); hold off;
 title(sprintf('Frame #%d',tp))
 linkaxes(a);


%%

%%

AutoGL = AutoGL(pathlength>0,:);
AutoVar = AutoVar(pathlength>0,:);
Step = 3;
TrueAutoGL = AutoGL(Lvalid==1,:);
Validpathlength = pathlength(pathlength>0);
ids = find(pathlength>0);
start = nan(size(TrueAutoGL,1),1);
for i = 1:size(TrueAutoGL,1)
    start(i) =  find(~isnan(TrueAutoGL(i,:)),1,'first');
end
[~,sidx] = sort(start);
TrueAutoGL = TrueAutoGL(sidx,:);
for i = 1:Step:sum(Lvalid==1)
    figure;
    s = min(i+2,size(TrueAutoGL,1));
    plot(TrueAutoGL(i:s,:)')
    xlabel('frame #'); ylabel('gl'); title(sprintf('%s:Gray level of correct tracks',Name));
    legend(num2str(sidx(i:s)))
end