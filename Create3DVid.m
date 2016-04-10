function Create3DVid(SaveDir)
%%
VisInfo = Load_Data(fullfile(SaveDir,'Visualize'),'*.tif');
SegInfo = Load_Data(fullfile(SaveDir,'Results'),'*.tif');
writerObj = VideoWriter(fullfile(SaveDir,'Results3D.avi'));
open(writerObj);
T = SegInfo.Frame_Num;
I = imread(SegInfo.Frame_name{T});
maxcell = max(I(:));
cents = nan(2,maxcell,T);
[Xv,Yv] = meshgrid(1:VisInfo.Height,1:VisInfo.Width);
Zv = ones(size(I));
for t = 1:T
    I = imread(SegInfo.Frame_name{t});
    Stats = regionprops(I,'Centroid');
    cents(:,1:length(Stats),t) = reshape([Stats.Centroid],2,[]);
    %All_Stats{t}= Stats;
end
%%
z = 1:T;
X = squeeze(cents(1,:,:))';
Y = squeeze(cents(2,:,:))';
figure(1); hold on;
c = varycolor(maxcell);
r = randperm(maxcell);
c= c(r,:);
set(gca,'YDir','Reverse')

view(gca,-32,36);
for t = 1:T
    hold on;
    for tt = 1:maxcell
        plot3(X(1:t,tt),Y(1:t,tt),z(1:t),'Color',c(tt,:),'LineWidth',2,'parent',gca);
    end
    
    IVis = imread(VisInfo.Frame_name{t});
    [Iind,map] = rgb2ind(IVis,1e3);
    surf(Xv,Yv,Zv*t,double(Iind),'LineStyle','none','parent',gca);colormap(map);alpha(0.8)
    xlim([1,size(I,2)]);
    ylim([1,size(I,1)]);
    zlim([1,T]);
    zlabel('Frame #','fontsize',12);
    hold off;
    pause(0.1);
    orig_mode = get(1, 'PaperPositionMode');
    set(1, 'PaperPositionMode', 'auto');
    img = hardcopy(1, '-Dzbuffer', '-r0');
    set(1, 'PaperPositionMode', orig_mode);
    frame = im2frame(img);
    writeVideo(writerObj,frame);
    cla;
    
    
end
close(writerObj);
