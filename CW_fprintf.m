function CW_fprintf(FID, Str, varargin)

if nargin == 2
  fprintf(FID, Str);
  fprintf(Str);
elseif nargin > 2
  Out = sprintf(Str, varargin{:});
  fprintf(FID, Out);
  fprintf(Out);
else
  error('not enough inputs');
end

return;