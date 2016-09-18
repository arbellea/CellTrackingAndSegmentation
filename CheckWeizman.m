function [varargout]=CheckWeizman(cellData,LinkData,data_path,extention,expr,resPath,resExt,resExp,Name)
if ~exist('Name','var')||isempty(Name)||~ischar(Name)
    Name = '';
end

%%
load(fullfile(resPath,'..','Link.mat'));
%%
Data = Load_Data(data_path,extention,expr) ;
ResData = Load_Data(resPath,resExt,resExp) ;
%%
b = 15;
H = ResData.Height;
W = ResData.Width;


%%
maxID = size(cellData,3);
Match = nan(max(maxID),min(size(cellData,1),ResData.Frame_Num));
GTgl = nan(max(maxID),min(size(cellData,1),ResData.Frame_Num));
GTvar = nan(max(maxID),min(size(cellData,1),ResData.Frame_Num));
AutoGL =nan(max(maxID),min(size(cellData,1),ResData.Frame_Num));
AutoVar = nan(max(maxID),min(size(cellData,1),ResData.Frame_Num));
%%
S = imread(ResData.Frame_name{end});
maxIDAuto = max(S(:));
xyt = shiftdim(cellData,1);
xyt(:,squeeze(xyt(1,:,:))<b | squeeze(xyt(1,:,:)>(W-b))) = nan;
xyt(:,squeeze(xyt(2,:,:))<b | squeeze(xyt(2,:,:)>(H-b))) = nan;
for t = 1:Data.Frame_Num
    %%
     I = double(imread(Data.Frame_name{t}));
     xy = xyt(:,:,t);
     
     availableCells = find(~isnan(xy(1,:)));
    if t<=ResData.Frame_Num
    S = imread(ResData.Frame_name{t});
     %for s = unique(S(S>0))'
     %   AutoGL2(s,t) = mean(I(S(:)==s));
     %   
     %   AutoVar2(s,t) = var(I(S(:)==s));
     %end
    end
   
    sizeS = size(I);
    
    for i = availableCells
        xi = round(xy(1,i));
        yi = round(xy(2,i));
        xx = max(xi-1,1):min(xi+1,W);
        yy = max(yi-1,1):min(yi+1,H);
        s = S(yy,xx);
        
        
        if t<=ResData.Frame_Num
            Match(i,t) = mode(s(:));
            %AutoGL(manualTrack.cell_id(i),t) = mean(I(S(:)== Match(manualTrack.cell_id(i),t)));
            %AutoVar(manualTrack.cell_id(i),t)= var(I(S(:)== Match(manualTrack.cell_id(i),t)));
        end
        Ic = I(yy,xx);
        GTgl(i,t)=mean(Ic(:));
        GTvar(i,t)=var(Ic(:));
    end
   
    %cents = sub2ind(size(S),manualTrack.centroid_col(manualTrack.timepoint==t),manualTrack.centroid_row(manualTrack.timepoint==t));
    %Match(manualTrack.cell_id(manualTrack.timepoint==t),t)=S(cents);
    %%
end
MatchOrig = Match;

%%
for l = numel(Link):-1:1
    for m = 1:numel(Link(l).Children)
      Match(Match== Link(l).Children(m)) = Link(l).Mother;
    end
end
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
pathlength = sum(~isnan(xyt(1,:,:)),3);

figure; 
histBins = 1:5:(max(pathlength));
[pathLengthHist,histBinLoc]=histc(pathlength(pathlength>0),[histBins-2.5,inf]);

pathLengthHist = pathLengthHist(1:end-1);
Lvalid = L(pathlength>0,end);
trackErrByLength = zeros(size(pathLengthHist));

for hb = unique(histBinLoc)
    
    trackErrByLength(hb) = sum(Lvalid(histBinLoc==hb)>1);
end

bar(histBins,pathLengthHist,1); hold on;
bar(histBins,trackErrByLength,0.8,'r'); hold off;
legend('Manual Track','Error')
xlabel('path lengths'); ylabel('freq'); title(sprintf('%s: Path Lebgths of Manual Tracks & Errors, %3.2f %% Accuracy',Name,sum(trackErrByLength)./sum(pathLengthHist)))
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
correctLinks = nan(2,numel(LinkData));
correctLinksAll = zeros(size(MatchMed,2),2);
sisLinks = 0;
for i = 1:numel(LinkData)
     if LinkData(i).Time>size(MatchMed,2)
        break;
    end
    %%
    
    time = LinkData(i).Time;
   
    correctLinks(:,i)=-1;
    mother = LinkData(i).Mother;
    timeVec = max(1,time-5):time-1;
    motherAuto = median(MatchMed(mother,timeVec));
    timeVec = time:min(size(MatchMed,2),time+4);
    daughterAuto = median(MatchMed([LinkData(i).Children],timeVec),2);
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
if tp<3
    fprintf('cant show $d', tp)
    
else

 MatchMed(cellid,tp-2:tp+2)
 mt = xyt(:,cellid,t)
 %mt = manualTrack(manualTrack.timepoint==tp&manualTrack.cell_id==cellid,1:4)
 I1 = imread(Data.Frame_name{tp-2}); I2 = imread(Data.Frame_name{tp-1});
 I3 = imread(Data.Frame_name{tp});
 I4 = imread(Data.Frame_name{tp+1}); I5 = imread(Data.Frame_name{tp+2});
 
 
 p = prctile(double(I1(:)),[1,99]);
 
 figure(100);
 
 %imshow(I1,p); a(1) = gca; title(sprintf('Frame #%d',tp-2))
 %figure(101); 
 subplot(1,12,[1:4]);
 imshow(I2,p); a(1) = gca; title(sprintf('Frame #%d',tp-1));hold on;scatter(xyt(1,cellid,tp-1),xyt(2,cellid,tp-1),'*g','linewidth',2); hold off;
 %figure(103);
 subplot(1,12,[9:12]);
 imshow(I4,p); a(2) = gca; title(sprintf('Frame #%d',tp+1));hold on;scatter(xyt(1,cellid,tp+1),xyt(2,cellid,tp+1),'*g','linewidth',2); hold off;
 %figure(104); 
 %imshow(I5,p); a(5) = gca; title(sprintf('Frame #%d',tp+2))
 %figure(102);
 subplot(1,12,[5:8]);
 imshow(I3,p); a(3) = gca;hold on;scatter(xyt(1,cellid,tp),xyt(2,cellid,tp),'*g','linewidth',2); hold off;
 title(sprintf('Frame #%d',tp))
 
 I1 = imread(visData.Frame_name{tp-2}); I2 = imread(visData.Frame_name{tp-1});
 I3 = imread(visData.Frame_name{tp});
 I4 = imread(visData.Frame_name{tp+1}); I5 = imread(visData.Frame_name{tp+2});
 
 
 p = [];%prctile(double(I1(:)),[1,99]);
 
 figure(101);
 
 %imshow(I1,p); a(1) = gca; title(sprintf('Frame #%d',tp-2))
 %figure(101); 
 subplot(1,12,[1:4]);
 imshow(I2,p); a(4) = gca; title(sprintf('Frame #%d',tp-1)) ; title(sprintf('Frame #%d',tp-1));hold on;scatter(xyt(1,cellid,tp-1),xyt(2,cellid,tp-1),'*g','linewidth',2); hold off;
 %figure(103);
 subplot(1,12,[9:12]);
 imshow(I4,p); a(5) = gca; title(sprintf('Frame #%d',tp+1)); title(sprintf('Frame #%d',tp+1));hold on;scatter(xyt(1,cellid,tp+1),xyt(2,cellid,tp+1),'*g','linewidth',2); hold off;
 %figure(104); 
 %imshow(I5,p); a(5) = gca; title(sprintf('Frame #%d',tp+2))
 %figure(102);
 subplot(1,12,[5:8]);
 imshow(I3,p); a(6) = gca;hold on;scatter(xyt(1,cellid,tp),xyt(2,cellid,tp),'*g','linewidth',2); hold off;
 title(sprintf('Frame #%d',tp))
 linkaxes(a);
end
 
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
