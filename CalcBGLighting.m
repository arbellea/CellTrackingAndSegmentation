function varargout = CalcBGLighting(I,varargin)
[X,Y] = meshgrid(1:size(I,1),1:size(I,2));
%X = X./max(X(:));
%Y = Y./max(Y(:));
if length(varargin)>=1;
    M = varargin{1};
    if isempty(M)
         M = true(size(I));
    end
else
    M = true(size(I));
end
if length(varargin)>=2
    p = 1:varargin{2};
else
    p = 1:2;
end
bias = 1;

if p>0
    XX = bsxfun(@(x,q) x.^q,X(:),p);
    YY = bsxfun(@(x,q) x.^q,Y(:),p);
    
    XXM = bsxfun(@(x,q) x.^q,X(M(:)),p);
    YYM = bsxfun(@(x,q) x.^q,Y(M(:)),p);
else
    XX = [];
    YY = [];
    XXM = [];
    YYM = [];
end
    
W = [XXM,X(M(:)).*Y(M(:)),YYM,ones(sum(M(:)),bias)];
W_all = [XX,X(:).*Y(:),YY,ones(numel(I),bias)];
T = W'*W;
A = inv(T)*W'*I(M(:));
Bvec = W_all*A;
B = (reshape(Bvec,size(I)));
ICorrected = I-B;
%ICorrected = (ICorrected-min(ICorrected(:)))./(max(ICorrected(:))-min(ICorrected(:)))*255;
varargout{1} = ICorrected;
varargout{2} = B;
varargout{3} = A;

