function [Tracking] = Initialize_Tracker(data,tagged_data,varargin)
if length(varargin)==1
    handles = varargin{1};
    GUI = true;
else
    GUI = false;
end
med_filt = vision.MedianFilter;
if strcmp(data.Type,'uint8')
    dens_x = 0:(2^8-1);
    Tracking.maxgray = 2^8-1;
else
    dens_x = 0:(2^16-1);
    Tracking.maxgray = 2^16-1;
end
DensCellPoints = [];
DensEdgePoints = [];
DensBGPoints = [];
Tracking.B  =0;
for t = 1:min(4,tagged_data.Frame_Num)

    whos Kalmans
    %% Pre-process Image
    if GUI
        
    set(handles.lbl_frame_num,'String',['Showing frame ',num2str(t),' out of: ',num2str(handles.stop_frame)]);
    end
    
    I = double(imread(data.Frame_name{t}));
    %per = prctile(I(:),[0.1,99.9]);
    %I = max(I,per(1));I = min(I,per(2));
    L = double(imread(tagged_data.Frame_name{t}));
    [I,B] =  CalcBGLighting2(I,L==0);
    if t==1
    
    prc = prctile(I(:),[sum(L(:)~=0)./numel(L),100-sum(L(:)~=0)./numel(L)]);
    Tracking.prc = prc;
    else
        prc = Tracking.prc;
    end
    I = min(max(I,prc(1)),prc(2));
    I = (I-min(I(:)))/(max(I(:))-min(I(:)))*Tracking.maxgray;
    
    LCells = L>0;
    LEdges = logical(imdilate(LCells,ones(3))- LCells);
    LBG = ~(LCells|LEdges);
    DensCellPoints = cat(1,DensCellPoints,I(LCells&I<Tracking.maxgray));
    DensBGPoints = cat(1,DensBGPoints,I(LBG));
    DensEdgePoints = cat(1,DensEdgePoints,I(LEdges));

    h_bg = hist(I(LBG),dens_x);
    [~,maxidx] = max(h_bg);
    
    %Ip = dens_x(maxidx)*ones(size(I)+20);
    %Ip(11:size(I,1)+10,11:size(I,2)+10) = I;
    %I = Ip;
    %Lp = zeros(size(L)+20);
    %Lp(11:size(L,1)+10,11:size(L,2)+10) = L;
    %L = Lp;
    uniqe_L = unique(L(L>0));
    %% Show Image
    if GUI
    axes(handles.axes1)
    imshow(I,[]);
    axes(handles.axes2)
    imshow(L,[]);
    end
    
    
    %% Create Kalman and predict
    if t==1
        states = Calculate_State(I,L);
        states = states;

        Kalmans = struct('ID',num2cell([states.ID]),'count',num2cell(1:size(states,2)),'enabled',true,'kalman',[],'num_props',0,...
            'prev_state',[],'children',[],'Contour',[],...
            'size',[],'HD',0.1,'new',true,'z_pred',[],'state_err',[],'U',[],'Mother',[],'Children',[],'cycle',0);
        Kalmans = arrayfun(@(s,k) fill_Kalman(s,k),states,Kalmans);
        Tracking.B = Tracking.B +B;
        Tracking.maxCellID = max(L(:));
        continue;
    end
    %% Predict and Correct
    prev_states = states;
    states =Calculate_State(I,L);
    Tracking.maxCellID = max(max(L(:)),Tracking.maxCellID);
    
    for n = 1:length(states)
        m = find([prev_states.ID]==states(n).ID);
        if isempty(m)
            mcount = max([Kalmans(:).count]);
            Kalmans(end+1).count = mcount+1;
            Kalmans(end) = fill_Kalman(states(n),Kalmans(end));
            Kalmans(end).ID = states(n).ID;
            Kalmans(end).HD = 0.1;
            Kalmans(end).num_props = 0;
            Kalmans(end).new = true;
            Kalmans(end).enabled = true;
            Kalmans(end).cycle = 0;
            continue;
        end
        states(n).kalman_state(6:10) =states(n).kalman_state(1:5)-prev_states(m).kalman_state(1:5);
        [z_pred, x_pred, P_pred] = predict(Kalmans(m).kalman);
        [z_corr,x_corr,P_corr] = correct(Kalmans(m).kalman,states(n).kalman_state);
        
        [p1y,p1x] = find(fullSingle(Kalmans(m).Contour)); [p2y,p2x]=find(fullSingle(states(n).Contour));
        p2x = p2x - z_corr(1)+Kalmans(m).prev_state(1) ;
        p2y = p2y - z_corr(2)+Kalmans(m).prev_state(2) ;
        p1 = [p1y,p1x];p2=[p2y,p2x];
        hd = HausdorffDist(p1,p2);
        if isempty(hd)
            hd = 0;
        end
        
        Kalmans(m).HD = Kalmans(m).HD*(t-2)/(t-1)+hd/(t-1) ;
        Kalmans(m).prev_state = z_corr';
        Kalmans(m).Contour = states(n).Contour;
        Kalmans(m).BW = states(n).BW;
        Kalmans(m).cycle = Kalmans(m).cycle+1;
    end
Tracking.B = Tracking.B +B;
    
end
u = (4/(3*min(numel(DensBGPoints)+numel(DensCellPoints))))^(1./5)*std(cat(1,DensCellPoints,DensBGPoints));
dens_cells = FastKDE(DensCellPoints,dens_x,u);
            dens_BG = FastKDE(DensBGPoints,dens_x,u);
%dens_cells = ksdensity(DensCellPoints,dens_x);
%dens_BG = ksdensity(DensBGPoints,dens_x);
%dens_edge = ksdensity(DensEdgePoints,dens_x);
mucells = mean(DensCellPoints);
mubg = mean(DensBGPoints);
%dens_BG(dens_BG<=dens_cells&dens_x<=mubg) = 1000*dens_cells(dens_BG<=dens_cells&dens_x<=mubg);
%dens_BG(dens_x<=mubg) = dens_BG(round(mubg));
%dens_cells(dens_cells==0&dens_x>=mucells) = 100*eps;
Tracking.dens_cells = dens_cells;
Tracking.dens_BG = dens_BG;
%Tracking.dens_edge = dens_edge;
Tracking.priorBG = sum(L(:)==0)./length(L(:));
Tracking.priorCell = 1-Tracking.priorBG;
Tracking.dens_x = dens_x;
Tracking.B = Tracking.B./t;

Tracking.Kalmans = Kalmans;
Tracking.current_t = t+1;
Tracking.med_filt = med_filt;


