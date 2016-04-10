function [save_dir_name,Kalmans] = Track_and_Segment(data,Tracking,Params,varargin)
rng(1);
GUI = false;
if length(varargin)==2
    GUI = true;
end
if length(varargin)>0&&isstruct(varargin{1})
    taggedData = varargin{1};
end

Link = struct();
Save_images = Params.Flags.WriteVideo;
SaveCheckPoints = Params.Flags.SaveCheckPoints;
LoadCheckPoints = Params.Flags.LoadCheckPoints;
segParams = Params.parameters;
save_debug = Params.Flags.SaveDebug;
deleteIfErr = Params.Flags.deleteIfErr;
t = datestr(now);
t(t==' ')='_';
t(t==':')='-';
if SaveCheckPoints||LoadCheckPoints
    t = 'CheckPoints';
end
if Save_images
    save_dir_name =fullfile('Outputs',sprintf('Results_%s_%s',Params.General.Name,t));
    save_dir_vis = fullfile(save_dir_name,'Visualize');
    save_dir_res = fullfile(save_dir_name,'Results');
    
    mkdir(save_dir_vis);
    mkdir(save_dir_res);
    mkdir(save_dir_name);
    
end
 save_dir_checkp = fullfile(save_dir_name,'CheckPoints');
if SaveCheckPoints
   
    mkdir(save_dir_checkp);
end
if any(save_debug)
    mkdir(fullfile(save_dir_name,'Debug'))
end
min_cell_size = 1;
Height = data.Height+20;
Width = data.Width+20;
Kalmans = Tracking.Kalmans;
stopFrame = Tracking.stop_frame;
try
    
    if LoadCheckPoints&& exist(save_dir_checkp,'dir')
        file_list = dir(fullfile(save_dir_checkp,'*.mat'));
        fileNames = {file_list.name};
        expr = 'CheckPoint_(\d+).mat';
        tokens = cellfun(@(str) regexp(str,expr,'tokens'),fileNames,'uniformoutput',false);
        emptyTokens = cellfun(@isempty, tokens);
        tokens = tokens(~emptyTokens);
        token = cellfun(@(token) str2double(token{1}{1}),tokens);
        if ~isempty(token) 
            if islogical(LoadCheckPoints)
            t = max(token);
            elseif isscalar(LoadCheckPoints)
                tokensSort = sort(token);
                t =tokensSort(find((tokensSort-LoadCheckPoints)<=0,1,'last'));
                if isempty(t)
                    t = tokensSort(1);
                end
            end
            stopFrame = Tracking.stop_frame;
            Tracking = load(fullfile(save_dir_checkp,sprintf('CheckPoint_%d.mat',t)));
            Kalmans = Tracking.Kalmans;
            if isfield(Tracking,'Link')
            Link = Tracking.Link;
            end
            Tracking.current_t = t;
        end
        
    end
    for t =Tracking.current_t:stopFrame;
        if SaveCheckPoints>0&& mod(t,SaveCheckPoints)==0
            
            Tracking.Link = Link;
            save(fullfile(save_dir_name,'CheckPoints',sprintf('CheckPoint_%d.mat',t)),'-struct','Tracking');
        end
        tStartFrame = tic;
        if GUI
            while get(handles.tb_pause_cont,'Value')
                pause(0.1)
            end
        end
        [~,mbgidx] = max(Tracking.dens_BG);
        mBG1 = Tracking.dens_x(mbgidx);
        I = double(imread(data.Frame_name{t}));
        I = I-Tracking.B;
        prc = Tracking.prc;
        I = min(max(I,prc(1)),prc(2));
        I = (I-min(I(:)))/(max(I(:))-min(I(:)))*Tracking.maxgray;
        I = step(Tracking.med_filt,I);
        %Ip = mBG1*ones(size(I)+20);
        %Ip(11:size(I,1)+10,11:size(I,2)+10) = I;
        %I = Ip;
        figure(1); imshow(I,[]);
        I_prev = double(imread(data.Frame_name{t-1}));
        I_prev = I_prev-Tracking.B;
        I_prev = min(max(I_prev,prc(1)),prc(2));
        I_prev = (I_prev-min(I_prev(:)))/(max(I_prev(:))-min(I_prev(:)))*Tracking.maxgray;
        I_prev = step(Tracking.med_filt,I_prev);
        
        %Ip = mBG1*ones(size(I_prev)+20);
        %Ip(11:size(I_prev,1)+10,11:size(I_prev,2)+10) = I_prev;
        %I_prev = Ip;
        %L = []; Debug.Input.Trakcing = Tracking;Debug.Input.Kalmans = Kalmans;Debug.Input.I = I;Debug.Input.I_prev = I_prev;Debug.Input.segParams = segParams;
        tSeg = tic;
        fprintf('Start Segmentation of frame %d...\n',t);
        [~,L,L_New_Cells,~,Kalmans,z_pred,cog_diff,Debug] = Fuzzy_Segmentation(Tracking,Kalmans,I,I_prev,segParams,any(save_debug));
        %Debug.Output.L = L; Debug.Output.L_New_Cells = L_New_Cells; Debug.Output.Kalmans = Kalmans; Debug.Output.z_pred = z_pred;Debug.Output.Debug = DebugTmp;
        timeSeg = toc(tSeg);
        fprintf('Done Segmentation of frame %d in %f seconds...\n',t,timeSeg);
        
        
        
        Obj_num = length(Kalmans);
        
        
        disabeledKalmans = Kalmans(~[Kalmans.enabled]);
        
        Kalmans = Kalmans([Kalmans.enabled]);
        tCalc = tic;
        states = Calculate_State(I,L,Kalmans);
        timeCalc = toc(tCalc);
        fprintf('Done Calculate State of frame %d in %0.3f seconds...\n',t,tCalc);
        tUp = tic;
        for n = 1:length(Kalmans);
            
            if isempty(states(n).kalman_state)
                Kalmans(n).state = Kalmans(n).z_pred';
                Kalmans(n).state(6:end)=0;
                Kalmans(n).HD = 10;
                Kalmans(n).num_props=endiKalmans(n).num_props+1;
                continue
            end
            %states(n).kalman_state(1:2) = states(n).kalman_state(1:2)-cog_diff;
            Kalmans(n).state =states(n).kalman_state';
            Kalmans(n).state_err(end+1,:) = Kalmans(n).z_pred-states(n).kalman_state;
            Kalmans(n).Contour = states(n).Contour;
            Kalmans(n).BW = states(n).BW;
            Kalmans(n).HD = states(n).HD;
            
            p=1;
            Kalmans(n).state(6:end) = p*(Kalmans(n).state(1:5)-Kalmans(n).prev_state(1:5))+(1-p)*(Kalmans(n).prev_state(6:end));
            Kalmans(n).state(6:7) = Kalmans(n).state(6:7)-cog_diff';
            if Kalmans(n).num_props
                Kalmans(n).state(6:end) = 0;
            end
            Kalmans(n).num_props = 0;
            Kalmans(n).size = states(n).size;
            Kalmans(n).cycle = Kalmans(n).cycle+1; 
        end
        timeUp = toc(tUp);
        fprintf('Done State update of frame %d in %f seconds...\n',t,timeUp);
        
        if ~isfield(Kalmans,'state')
            Kalmans(1).state =[];
        end
        
        if ~isempty(L_New_Cells)
            tnCalc = tic;
            states =Calculate_State(I,L_New_Cells);
            timeNCalc = toc(tnCalc);
            fprintf('Done New State Calculation of frame %d in %f seconds...\n',t,timeNCalc);
            tnUp = tic;
            uniqe_L = unique(L_New_Cells(L_New_Cells>0));
            states = states(uniqe_L);
            m =length(Kalmans);
            for n = 1:size(states,2);
                mu = states(n).kalman_state(1:2);
                if mu(1)<=1||mu(2)<=1||mu(1)>(Width)||mu(2)>(Height)
                    continue;
                end
                m = m+1;
                Kalmans(m).new = true;
                Kalmans(m).count =m;
                Kalmans(m).ID =Tracking.maxCellID +1;
                Tracking.maxCellID  = Tracking.maxCellID +1;
                Kalmans(m).enabled = true;
                Kalmans(m).kalman = Create_New_Kalman(states(n).kalman_state,8^2*eye(10),2*eye(10));
                Kalmans(m).num_props = 0;
                Kalmans(m).Children =[];
                Kalmans(m).prev_state = states(n).kalman_state;
                Kalmans(m).state = states(n).kalman_state;
                Kalmans(m).HD =100;
                Kalmans(m).size = states(n).size;
                Kalmans(m).z_pred = zeros(1,10);
                Kalmans(m).state_err = [];
                Kalmans(m).BW = states(n).BW;
                Kalmans(m).Contour = states(n).Contour;
                Kalmans(m).Mother = [];
                Kalmans(m).cycle = 1;
                L(fullSingle(Kalmans(m).BW)) = Kalmans(m).ID;
            end
            timeNUp = toc(tnUp);
            fprintf('Done New State Update of frame %d in %f seconds...\n',t,timeNUp);
        end
        if ~isempty(disabeledKalmans)
            tmitLink = tic;
            for n = 1:length(disabeledKalmans)
                if ~isempty(disabeledKalmans(n).Children);
                    motherID = disabeledKalmans(n).ID;
                    children = zeros(length(disabeledKalmans(n).Children),1);
                    for child = 1:length(disabeledKalmans(n).Children)
                        c = round(disabeledKalmans(n).Children(child).Centroid);
                        Kalmans([Kalmans.ID]==L(disabeledKalmans(n).Children(child).PixelIdxList(1))).Mother=disabeledKalmans(n).ID;
                        children(child) =L(disabeledKalmans(n).Children(child).PixelIdxList(1)) ;
                    end
                    if ~isfield(Link,'Mother');
                    Link(1).Mother = motherID;
                    Link(1).Children = children;
                    Link(1).Time = t;
                    else
                    Link(end+1).Mother = motherID;
                    Link(end).Children = children;
                    Link(end).Time = t;
                    end
                end
            end
            timeMitLink = toc(tmitLink);
            fprintf('Done Mitosis Link of frame %d in %f seconds...\n',t,timeMitLink);
        end
        tKUp = tic;
        for n = 1:length(Kalmans)
            if  Kalmans(n).enabled
                [z_corr(:,n),~,~] = correct(Kalmans(n).kalman,Kalmans(n).state);
                Kalmans(n).prev_state = z_corr(:,n);
                mes_mus(n,:)=Kalmans(n).state(1:2)';%.*[Width,Height];
                if mes_mus(n,1)<=1||mes_mus(n,2)<=1||mes_mus(n,1)>=(Width)||mes_mus(n,2)>=(Height)||Kalmans(n).size<min_cell_size % mes_mus(n,1)<=10||mes_mus(n,2)<=10||mes_mus(n,1)>=(Width-10)||mes_mus(n,2)>=(Height-10)||Kalmans(n).size<min_cell_size
                    Kalmans(n).enabled = false;
                    
                end
            end
            
            
        end
        z_mat= cell2mat(z_pred')';
        pred_mus = z_mat(1:2,:)';
        cor_mus = z_corr(1:2,:)';
        timeKUp = toc(tKUp);
        fprintf('Done Kalman Update of frame %d in %f seconds...\n',t,timeKUp);
        [Kalmans([Kalmans.new]).new] = deal(false);
        tSave = tic;
        if Save_images
            
            [~,fname,ext]=fileparts(data.Frame_name{t});
            frame_name = sprintf('Seg_%s%s',fname,ext);
            frame_name = fullfile(save_dir_res,frame_name);
            imwrite(uint16(L(1:end,1:end)),frame_name);%imwrite(uint16(L(11:end-10,11:end-10)),frame_name);
            
        end
        timeSave = toc(tSave);
        
        fprintf('Done Save frame %d in %f seconds...\n',t,timeSave);
        if GUI
            set(handles.lbl_frame_num,'String',['Showing frame ',num2str(t),' out of: ',num2str(handles.stop_frame)]);
            axes(handles.axes1)
            imshow(img);
            axes(handles.axes2) %#ok<*LAXES>
            h=imagesc(L);colormap(jet);colorbar; hold on;
            scatter(cor_mus([Kalmans.enabled],1),cor_mus([Kalmans.enabled],2),'g','filled')
            scatter(mes_mus([Kalmans.enabled],1),mes_mus([Kalmans.enabled],2),'r','filled')
            scatter(pred_mus(:,1),pred_mus(:,2),'b','filled')
            hold off
            set(h,'ButtonDownFcn',@clicked);
            pause(0.001);
        end
        
        if ismember(t,save_debug)
            Kalmans_debug= Kalmans; %#ok<*NASGU>
            Seg_Debug = Debug;
            save(fullfile(save_dir_name,'Debug',sprintf('Kalmans_debug_frame_%d.mat',t)),'Kalmans_debug','Seg_Debug','Tracking','-v7.3');
        end
        %I = I(11:end-10,11:end-10);
        %L = L(11:end-10,11:end-10);
        tKDE = tic;
            DensCellPoints =[];
            DensBGPoints = [];
            DensEdgePoints = [];
        if mod(t,1)==0&&exist('DensCellPoints','var')&&exist('DensBGPoints','var')
            LCells = L>0;
            %LEdges = logical(imdilate(LCells,ones(3))- LCells);
            LBG = ~(LCells);
            DensCellPoints = cat(1,DensCellPoints,I(LCells&I<Tracking.maxgray));
            DensBGPoints = cat(1,DensBGPoints,I(LBG));
            %DensEdgePoints = cat(1,DensEdgePoints,I(LEdges));
            u = (4/(3*min(numel(DensBGPoints)+numel(DensCellPoints))))^(1./5)*max(std(DensCellPoints),std(DensBGPoints));
            dens_cells = FastKDE(DensCellPoints,Tracking.dens_x,u);
            dens_BG = FastKDE(DensBGPoints,Tracking.dens_x,u);
            zz = find(dens_cells==0&dens_cells==0);
            %dens_edge = FastKDE(DensEdgePoints,Tracking.dens_x,u);
            % dens_cells = ksdensity(DensCellPoints,Tracking.dens_x,'bandwidth',u);
            % dens_BG = ksdensity(DensBGPoints,Tracking.dens_x,'bandwidth',u);
            % dens_edge = ksdensity(DensEdgePoints,Tracking.dens_x,'bandwidth',u);
            zd = knnsearch([DensCellPoints;DensBGPoints],zz'); 
            dens_cells(zz(zd<=numel(DensCellPoints))) = eps;
            dens_BG(zz(zd>numel(DensCellPoints))) = eps;
            mucells = mean(DensCellPoints);
            mubg = mean(DensBGPoints);
            %dens_BG(Tracking.dens_x<=mubg) = dens_BG(round(mubg));
            %dens_cells(dens_cells==0&Tracking.dens_x>=mucells) = 100*eps;
            %Tracking.dens_edge = dens_edge;
            Tracking.dens_cells = dens_cells;
            Tracking.dens_BG = dens_BG;
            Tracking.priorBG = sum(L(:)==0)./length(L(:));
            Tracking.priorCell = 1-Tracking.priorBG;
        elseif mod(t,1)>=0
            LCells = L>0;
            LEdges = logical(imdilate(LCells,ones(3))- LCells);
            LBG = ~(LCells|LEdges);
            DensCellPoints = I(LCells);
            DensBGPoints = I(LBG);
            %DensEdgePoints = I(LEdges);
        end
        timeKDE = toc(tKDE);
        fprintf('Done KDE Update %d in %f seconds...\n',t,timeKDE);
        %[~,B] = CalcBGLighting(I(11:end-10,11:end-10),L(11:end-10,11:end-10)==0);
        %Tracking.B = B;
        tEndFrame = toc(tStartFrame);
        Kalmans = Kalmans([Kalmans.enabled]);
        Tracking.Kalmans = Kalmans;
        whos Kalmans
        fprintf('Elapsed time for Frame %d is %2.5f\n',t,tEndFrame)
        
    end
    if Save_images
        save(fullfile(save_dir_name,'Link.mat'),'Link');
    end
catch err
    if ismember(t,save_debug)
        %Kalmans_debug= Kalmans; %#ok<*NASGU>
        %Seg_Debug = [];
        mkdir([save_dir_name,'_Debug']);
        save(fullfile([save_dir_name,'_Debug'],sprintf('Kalmans_debug_frame_%d.mat',t)),'-v7.3');
    end
    if deleteIfErr&&isdir(save_dir_name)
        try
            rmdir(save_dir_name,'s');
        catch err1
            disp('Could not delete output directory')
        end
        
    end
    rethrow(err);
end

end

function Position = clicked(hObject,~)
return;
end