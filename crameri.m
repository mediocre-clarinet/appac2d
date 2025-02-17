function cmap = crameri(ColormapName,varargin)
filename = mfilename('fullpath');
filepath = fileparts( filename );
S = load([filepath '/colormaps/' ColormapName '.mat']);

N = 256;
if nargin == 2
    N = varargin{1};
end

position = linspace(1,256,N).';
intBelow = floor(position);
intAbove = intBelow + 1;
intAbove(intAbove>256) = 256;
drgb = S.(ColormapName)(intAbove,:) - S.(ColormapName)(intBelow,:);
cmap = S.(ColormapName)(intBelow,:) + drgb.*(position - intBelow);