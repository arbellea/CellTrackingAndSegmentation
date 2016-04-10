function [Ureturn,L,L_New_Cell,z,Kalmans,z_pred,DEBUG] = Fuzzy_SegmentationPar(Tracking,Kalmans,I,I_prev,params,save_debug) 
[Height,Width]=size(I);
[X,Y] = meshgrid(1:Width,1:Height);

XY = [Y(:),X(:)]';
cog_diff = (flipud(XY)*I(:)./sum(I(:))-flipud(XY)*I_prev(:)./sum(I_prev(:)))';
if isfield(params,'patchSize')
    patchSize = params.patchSize;
else
    patchSize = 100;
end
if isfield(params,'minErr')
    min_err = params.minErr;
else
    min_err = 15;
end
if isfield(params,'maxItr')
    max_itr = params.maxItr;
else
    max_itr =5;
end

if isfield(params,'minCellSize')
    min_cell_size = params.minCellSize;
else
    min_cell_size = 100;
end

if isfield(params,'solidityThr')
    Solidity_Thr = params.solidityThr;
else
    Solidity_Thr = 0.8;
end

if isfield(params,'cellPrior')
alpha = params.cellPrior;
else
alpha = 0.5;
end

q = 1;
enabled = [Kalmans.enabled]'&([Kalmans.size]>=ceil(min_cell_size/10))';
enKalmans = Kalmans(enabled);
if save_debug
    DEBUG = cell(max_itr);
else
    DEBUG = [];
end
try
    [z_pred,~, P_pred] = arrayfun(@(x) predict(x.kalman),enKalmans,'UniformOutput',false);
    z_pred_mat = cell2mat(shiftdim(z_pred,1));
    
    new_enabled = ~(z_pred_mat(:,1)<3|z_pred_mat(:,2)<3|z_pred_mat(:,1)>Width-3|z_pred_mat(:,2)>Height-3);
    enabled(enabled) = new_enabled;
    z_pred = z_pred(new_enabled);
    
    z_pred = cellfun(@(z) z+[cog_diff,0,0,0,0,0,0,0,0],z_pred,'UniformOutput',false);
    mu = cellfun(@(z) flipud(z(1:2)'), z_pred,'UniformOutput',false);
    P_pred = P_pred(new_enabled);
    enKalmans = enKalmans(new_enabled);
    [enKalmans.z_pred] = deal(z_pred{:});
    cell_size = arrayfun(@(x) x.size,enKalmans,'UniformOutput',false);
    BWs = arrayfun(@(x) x.BW,enKalmans,'UniformOutput',false);
    HDs = arrayfun(@(x) max(x.HD,0.01),enKalmans,'UniformOutput',false);
    last_mu = arrayfun(@(x) (flipud(x.prev_state(1:2))),enKalmans,'UniformOutput',false);
    hg = fspecial('gauss');
    
    I = imfilter(double(I), hg, 'replicate');
    gradI = calc_Grad_I(I);
    invG = max(1./(1+(gradI./(10*std(gradI(:))))).^q,0.1);
    
    
    %%
    
    %%
    
    
    itr = 0;
    err = min_err+1;
    PBG = Tracking.dens_BG(round(I)+1);
    P_Cell = Tracking.dens_cells(round(I)+1);
    nBG = alpha*P_Cell./(alpha*P_Cell+(1-alpha)*PBG);
%%
    NKalmans = sum(enabled);
    U = cell(1,NKalmans);
    %BWs_moved =cell(1,NKalmans);
    %SDF_moved = cell(1,NKalmans);
    %Phi_moved = cell(1,NKalmans);
    %Phi_norm =cell(1,NKalmans);
    %Phi_cropped =cell(1,NKalmans);
    %gradI_cropped = cell(1,NKalmans);
    %invG_cropped = cell(1,NKalmans);
    %nBG_cropped = cell(1,NKalmans);
    %mu_cropped = cell(1,NKalmans);
    %Speed_cropped = cell(1,NKalmans);
    %D_cropped =cell(1,NKalmans);
    %P_cropped = cell(1,NKalmans);
    P = cell(1,NKalmans);
    %U_cropped =cell(1,NKalmans);
    valid = true(1,NKalmans);
    while err>=min_err&&itr<max_itr
        
        
        itr=itr+1;
        valid = true(1,NKalmans);
        parfor i = 1:NKalmans
            [Phi_moved{i},~,~,~,~,valid(i)] = movePhi(fullSingle(BWs{i}),mu{i},last_mu{i},HDs{i},patchSize);
        end
            if any(~valid)
                
                HDs = HDs(valid);
                Phi_moved = Phi_moved(valid);
                %BWs_moved = BWs_moved(valid);
                BWs = BWs(valid);
                %BW_cropped = BW_cropped(valid);
                %SDF_moved_cropped = SDF_moved_cropped(valid);
                %Phi_moved_cropped = Phi_moved_cropped(valid);
                mu = mu(valid);
                last_mu = last_mu(valid);
                cell_size = cell_size(valid);
                enKalmans = enKalmans(valid);
                NKalmans = length(enKalmans);
                P_pred = P_pred(valid);
                enabled(enabled) = logical(valid);
                if exist('U','var')
                    U = U(valid);
                end
                if exist('P','var')
                    P = P(valid);
                end
                
            end
            
        Phimat = cat(3,Phi_moved{:});
        Sphi = sum(Phimat,3) + (1-max(Phimat,[],3));
        Sphi(Sphi(:)==0)=1;
        parfor i = 1:NKalmans
            Phi_norm = Phi_moved{i}./Sphi;% Phi_norm =  cellfun(@(phi) phi./Sphi,Phi_moved,'uniformoutput',0);
            Phi_cropped = CropImage(Phi_norm,mu{i},patchSize,patchSize); %Phi_cropped = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_norm,mu,'uniformoutput',0);
            invG_cropped = CropImage(invG,mu{i},patchSize,patchSize);%invG_cropped = cellfun(@(m) CropImage(invG,m,patchSize,patchSize),mu,'uniformoutput',0);
            nBG_cropped = CropImage(nBG,mu{i},patchSize,patchSize);%nBG_cropped = cellfun(@(m) CropImage(nBG,m,patchSize,patchSize),mu,'uniformoutput',0);
        if itr==1
            mu_cropped = FindMostLiklymu(nBG_cropped,Phi_cropped); %mu_cropped = cellfun(@FindMostLiklymu,nBG_cropped,Phi_cropped,'uniformoutput',0);
        else
            U_cropped = CropImage(U{i},mu{i},patchSize,patchSize);%U_cropped = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),U,mu,'uniformoutput',0);
            mu_cropped = FindMuCOM(U_cropped);%mu_cropped = cellfun(@FindMuCOM,U_cropped,'uniformoutput',0);
        end
            Speed_cropped = max((nBG_cropped).*(invG_cropped),1e-8); %Speed_cropped = cellfun(@(nbg,invg) max((nbg).*(invg),1e-8),nBG_cropped,invG_cropped,'uniformoutput',0);
            Speed_cropped = fixSpeedIm(Speed_cropped,mu_cropped,P_pred{i});%Speed_cropped = cellfun(@fixSpeedIm,Speed_cropped,mu_cropped,P_pred,'uniformoutput',0);
            D_cropped = CreateDiffDist(Speed_cropped,mu_cropped);%D_cropped = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped,'uniformoutput',0);
            P_cropped = Phi_cropped./(D_cropped+1);%P_cropped = cellfun(@(phi,d) phi./(+1),Phi_cropped,D_cropped,'uniformoutput',0);
            P{i} =PadImage(P_cropped,mu{i},patchSize,patchSize,Height,Width,0);%P = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),P_cropped,mu,'uniformoutput',0);
        end
        Pmat = cat(3,P{:});
        S = sum(Pmat,3)+1-nBG;
        S(S==0)=eps;
        UBG = (1-nBG+eps)./S;
        [U,~] = cellfun(@(p,m) NormalizeAndCropU(p,S,m,patchSize),P,mu,'uniformoutput',0);
        mu = cellfun(@(u) XY*u(:)./max(sum(u(:)),eps),U,'UniformOutput',false);
        prev_size = cell_size;
        cell_size = num2cell(cellfun(@(u) round(sum(u(:))),U));
        errs = sqrt(sum((cell2mat(cell_size)-cell2mat(prev_size)).^2,1));
        [err,erridx] = max(errs);
        if any([cell_size{:}]==0)
            remove_cell = [cell_size{:}]>0;
            BWs = BWs(remove_cell);
            P_pred = P_pred(remove_cell);
            HDs = HDs(remove_cell);
            %BWs_moved = BWs_moved(remove_cell);
            mu = mu(remove_cell);
            P = P(remove_cell);
            U = U(remove_cell);
            last_mu = last_mu(remove_cell);
            [enKalmans(~remove_cell).enabled] = deal(false);
            enabled(enabled) = logical(remove_cell);
	        enKalmans= enKalmans(remove_cell);
            NKalmans = length(enKalmans);
            cell_size = cell_size(remove_cell);
        end
        
    end
    
    Us = cellfun(@(u) sparseSingle(u),U,'UniformOutput',false);
    [enKalmans(:).U] = deal(Us{:});
    
    Kalmans(enabled) = enKalmans;
    enabledcell = num2cell(logical(enabled));
    [Kalmans(:).enabled] = deal(enabledcell{:});
    FullU = U;
    FullU{end+1} = UBG;
    
    Ureturn = cat(3,FullU{:});
    
    FullP = P;
    FullP{end+1}=PBG;
    
    
    z = cell2mat(shiftdim(FullP,-1));
    
    [~,L] = max(Ureturn,[],3);
    L(L>length(P))=0;
    labels = [enKalmans.ID];
    L(L>0) = labels(L(L>0));
    NotBG = false(size(I)); 
    nBGinvG = nBG.*invG;
    NotBG(L==0) = logical(nBGinvG(L==0)>0.5);
    NotBG = imerode(NotBG,ones(3));
    NotBGFull = NotBG;
    NotBG_Candidates = [];
    for l = unique(L(L>0))'
        BW = L==l;
        NotBG_Candidates_l = regionprops(BW,'Solidity','Area','PixelIdxList','Centroid');
        NotBG_Candidates_l =  NotBG_Candidates_l([NotBG_Candidates_l.Solidity]>Solidity_Thr&[NotBG_Candidates_l.Area]>min_cell_size);
        if length(NotBG_Candidates_l)<=1
            continue
        end
        flags = false(length(NotBG_Candidates_l),1);
        for i = 1:length(NotBG_Candidates_l)
            [idxy,idxx]= ind2sub([Height,Width],NotBG_Candidates_l(i).PixelIdxList);
            if any(idxx==11)||any(idxy==11)||any(idxx==Width-10)||any(idxy==Height-10)
                flags(i) = true;
                NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_l(i));
                NotBGFull(NotBG_Candidates_l(i).PixelIdxList) = 1;
                L(NotBG_Candidates_l(i).PixelIdxList) = 0;
            end
        end
        
        if sum(~flags)>1
            NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_l(~flags));
            NotBGFull(L==l) = 1;
            L(L==l) = 0;
            Kalmans([Kalmans(:).ID]==l).enabled = false;
        end
    end
    
    NotBG_Candidates_bg = regionprops(NotBG,'Solidity','Area','PixelIdxList','Centroid');
    NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_bg);
    New_Cells = NotBG_Candidates([NotBG_Candidates.Solidity]>Solidity_Thr&[NotBG_Candidates.Area]>min_cell_size);
    if ~isempty(New_Cells)
        BW_New_Cell = arrayfun(@(x) Create_New_Cell_L(NotBGFull,x.PixelIdxList),New_Cells,'UniformOutput',false);
        L_New_Cell = any(cell2mat(shiftdim(BW_New_Cell,-2)),3);
        L_New_Cell = bwlabel(L_New_Cell);
    else
        L_New_Cell=[];
    end
catch err
    save('Debug.mat','DEBUG','-v7.3');
    rethrow(err);
end

end

function BW = Create_New_Cell_L(BW_Full,Idx)
BW = zeros(size(BW_Full));
BW(Idx(:)) = BW_Full(Idx(:));
end
function Cropped = CropImage(I,c,h,w,varargin)
y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),size(I,1));
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),size(I,2));
Cropped = I(y1:y2,x1:x2);
end

function Cropped_mu = CropMu(c,h,w) %#ok<*DEFNU>
y1 = max(round(c(1)-h/2),1);
x1 = max(round(c(2)-w/2),1);
Cropped_mu = round(c-[y1-1;x1-1]);
end

function PadI = PadImage(I,c,h,w,H,W,padv)
PadI = padv*ones(H,W);

y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),H);
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),W);

if (y2-y1+1)==size(I,1)&&(x2-x1+1)==size(I,2)
    PadI(y1:y2,x1:x2) = I;
end

end

function gradI = calc_Grad_I(I)

hy = [-1;1]%fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradI = Ix.^2 + Iy.^2;
end

function mu = FindMostLiklymu(PG,Phi)
%[X,Y] = meshgrid(1:size(PG,2),1:size(PG,1));
if size(PG)==size(Phi)
    [muy,mux]=find(PG.*Phi==max(PG(:).*Phi(:)),1);
    %mu = round([Y(:),X(:)]'*(PG(:).*Phi(:))./sum(PG(:).*Phi(:)));
    mu = [muy;mux];
else
    [muy,mux]=find(PG==max(PG(:)),1);
    mu = [muy;mux];
end

end

function mu = FindMuCOM(P)
[X,Y] = meshgrid(1:size(P,2),1:size(P,1));
mu = round([Y(:),X(:)]'*P(:)./sum(P(:)));
end


function Cropped = CropImage_prev(I,c,c_n,h,w,fillval)
yn1 = max(round(c_n(1)-h/2),1);
yn2 = min(round(c_n(1)+h/2),size(I,1));
xn1 = max(round(c_n(2)-w/2),1);
xn2 = min(round(c_n(2)+w/2),size(I,2));

Cropped = fillval.*ones(yn2-yn1+1,xn2-xn1+1);

y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),size(I,1));
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),size(I,2));

Y1 = 1; Y2 = size(Cropped,1);
X1 = 1; X2 = size(Cropped,2);
if (y2-y1)<(yn2-yn1)
    if y1==1
        Y1 = (yn2-yn1)-(y2-y1)+1;
    else
        Y2 =(y2-y1)+1;
    end
elseif (y2-y1)>(yn2-yn1)
    if yn1==1
        y1 = y2-Y2+1;
    else
        y2 = y1 + Y2-1;
    end
end

if (x2-x1)<(xn2-xn1)
    if x1==1
        X1 = (xn2-xn1)-(x2-x1)+1;
    else
        X2 = (x2-x1)+1;
    end
elseif (x2-x1)>(xn2-xn1)
    if xn1==1
        x1 = x2-X2+1;
    else
        x2 = x1 + X2-1;
    end
end
Cropped(Y1:Y2,X1:X2) = I(y1:y2,x1:x2);
end

function [DiffD,Dg,De] = CreateDiffDist(SE,mu)
[X,Y] = meshgrid(1:size(SE,2),1:size(SE,1));
De = max(sqrt((Y-mu(1)).^2+(X-mu(2)).^2),eps);
if any(isnan(SE(:)))||any(isnan(mu))||any(mu<1)||mu(1)>size(SE,1)||mu(2)>size(SE,2)
    error('Fast Marching will crash!')
end
Dg = max(msfm2d(SE,mu,true,true),eps);
DiffD = max(Dg-De,0);
end

function [Speed_G] = fixSpeedIm(Speed,mu,C)
[X,Y] = meshgrid(1:size(Speed,2),1:size(Speed,1));
D = (X(:)-mu(2)).^2./C(1,1)+(Y(:)-mu(1)).^2./C(2,2);
P = reshape(exp(-0.5*D),size(Speed));
Speed_G = max(Speed,P);
end

function [Phi_smoothed] = smoothPhi(Phi,C)
w = max(C(1,1),C(1,2))*5;
[X,Y] = meshgrid(1:w,1:w);
mu = round(size(X)./2);
D = (X(:)-mu(2)).^2./C(1,1)+(Y(:)-mu(1)).^2./C(2,2);
P = reshape(exp(-0.5*D),size(X));
P = P./sum(P(:));
Phi_smoothed = imfilter(Phi,P,'same','conv');
end

function [bwout] = moveBW(bw,m,lastm)

bwout = zeros(size(bw));
[y,x] = find(bw);
y = round(y-lastm(1)+m(1));
x = round(x-lastm(2)+m(2));
invalid = (y<1|y>size(bw,1))|(x<1|x>size(bw,2));
y(invalid) = [];
x(invalid) = [];
idx = sub2ind(size(bw),y,x);
bwout(idx) = 1;
end

function [Phi_moved,BWs_moved,BW_cropped,SDF_moved_cropped,Phi_moved_cropped,valid] = movePhi(BWs,mu,last_mu,HD,patchSize)
BWs_moved=  moveBW(BWs,mu,last_mu);
[H,W] = size(BWs_moved);
BW_cropped= CropImage(BWs_moved,mu,patchSize,patchSize);
SDF_moved_cropped = SDF(BW_cropped);
Phi_moved_cropped = 1./(1+exp(-2*SDF_moved_cropped./(HD.*pi).*sqrt(3)));
Phi_moved = PadImage(Phi_moved_cropped,mu,patchSize,patchSize,H,W,0);
valid = any(BW_cropped(:));
end

function [U,U_cropped] = NormalizeAndCropU(P,S,mu,patchSize)
U = P./S;
U_cropped= CropImage(U,mu,patchSize,patchSize);
end
