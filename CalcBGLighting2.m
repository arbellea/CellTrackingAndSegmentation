function varargout = CalcBGLighting2(I,varargin)
h = fspecial('gaussian',210,70);
B = imfilter(I,h,'replicate');
ICorrected = I-B;

if length(varargin)>=1;
    M = varargin{1};
    if isempty(M)
         M = true(size(I));
    end
else
    M = true(size(I));
end
%ICorrected(M) = I(M)-B(M);
%ICorrected(M) = I(M)-B(M);

varargout{1} = ICorrected;
varargout{2} = B;


