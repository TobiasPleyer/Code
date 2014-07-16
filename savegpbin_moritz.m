%% Write gnuplot 3D binary file from a cartesian matrix.
%%	z = rows(z) x columns(z) = length(x) x length(y)
%%
%% Gnuplot binary file data organization:
%% <N> <y1> .. <yN>  <x1> <row of z> <x2> <row of z> ...
%% where all numbers are floats.
%%
%% Author: Petr Mikulik
%% Version: June 2004 (docs improved)
%% Original version: October 2001
%% License: Public domain

function savegpbin_moritz ( x, y, z, filename )
  % check input arguments
  err = 0;
  if (nargin ~= 4) err = 1;
  elseif (~isvector(x) || ~isvector(y)) 
    error('x and y must be vectors\n');
  elseif (any(size(z) - [length(x) length(y)]))
    % Cannot do: (size(z) ~= [length(x) length(y)])  because result is matrix!
    error('z must be matrix of size length(x) x length(y)\n');
  end

  if (err) 
    a=[	'\tUSAGE: savegpbin( x, y, z, filename )\n', ...
	'\tWhere x, y are vectors and z is cartesian matrix length(x) x length(y)\n', ...
	'\t(z rotated by +90 deg is cartesian in x x y). For example:\n', ...
	'\tx=1:5; y=1:8; [xx,yy]=meshgrid(y,x); z=xx.*yy; savegpbin(x,y,z,''a.dat'')\n' ...
	];
    fprintf(a);
    return
  end
  
  if any(size(x)-[1,length(x)]) x=x'; end
  if any(size(y)-[1,length(y)]) y=y'; end

  % working code
  f = fopen(filename, 'wb');
  [tmp, nx] = size(x);
  fwrite(f, nx, 'float');
  fwrite(f, x, 'float');
  fwrite(f, [y; z], 'float');
  fclose(f);

  [nzx, nzy] = size(z);
%  sajz = 4 * ( 1 + nzx*nzy + nzx + nzy );
%  sprintf('File size should be %i B\n',sajz)

% endfunction

% eof savegpbin.m
