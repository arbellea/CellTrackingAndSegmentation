function output_3D_vid(indir,vidfile)
%%
 writerObj = VideoWriter([vidfile(1:end-4),'_tracks.avi']);
    writerObj.FrameRate=2;
    open(writerObj);
Vidobj = VideoReader(vidfile);
finfo = dir(fullfile(indir,'*.tif'));
name = arrayfun(@(f) f.name(5:end-4),finfo,'UniformOutput',false);
sorted_frames = sort(str2double(name));
l= length(sorted_frames);
fname = sprintf('Seg_%1.0f.tif',sorted_frames(end));
I = imread(fullfile(indir,fname));
maxcell = max(I(:));
cents = nan(2,maxcell,l);
[Xv,Yv] = meshgrid(1:size(I,2),1:size(I,1));
Zv = ones(size(I));
for i = 1:l
    fname = sprintf('Seg_%1.0f.tif',sorted_frames(i));
    I = imread(fullfile(indir,fname));

    Stats = regionprops(I,'Centroid');
    cents(:,1:length(Stats),i) = reshape([Stats.Centroid],2,[]);
    All_Stats{i}= Stats;
    
end
%%
z = 1:l;
X = squeeze(cents(1,:,:))';
Y = squeeze(cents(2,:,:))';
figure(1); hold on;
c = varycolor(maxcell);
r = randperm(maxcell);
c= c(r,:);
hs = [];
set(gca,'YDir','Reverse')

view(gca,-32,36);
for j = 1:l
hold on;
    for i = 1:maxcell
        plot3(X(1:j,i),Y(1:j,i),z(1:j),'Color',c(i,:),'LineWidth',2,'parent',gca);
    end

IVid = read(Vidobj,j);
[Iind,map] = rgb2ind(IVid,1e3);
hs = surf(Xv,Yv,Zv*j,double(Iind),'LineStyle','none','parent',gca);colormap(map);alpha(0.8)
xlim([1,size(I,2)]);
ylim([1,size(I,1)]);
zlim([1,l]);
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
