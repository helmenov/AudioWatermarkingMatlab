function str = dss_dec(file, p)

% Decoding echo hiding
%   extract embedded data from echo hided WAV file
% Type
%   str = echo_dec;
% on the command line of MATLAB or Octave 2.1.x.
%   It reads a host signal 'echoenc.wav' and extract embedded data to str.

%  Copyright (C) 2007 by Akira Nishimura
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin < 1
  file = 'dssenc.wav';
end

if nargin < 2
  p.trate = 4;         % tip rate
  p.drate = 8192;      % data rate
  p.cf = 0;            % carrier freq.
%  p.lev = -30;         % relative to signal level
  p.key = 1374;        % key of pseudo-random sequence; seed of randomize
end

rand('seed', p.key);           % set seed of random sequence
a = rand(p.drate/p.trate, 1);
pnseq = (a > 0.5)*2-1;        % pseudo-random sequence consisting of +1 and -1 
pnsig = reshape(ones(p.trate, 1)*pnseq', p.drate, 1);

[sig, sf] = wavread(file);
len = length(sig);
nframe = floor(len/p.drate);

xsig = reshape(sig(1:p.drate*nframe,1), p.drate, nframe);

for k=1:nframe
  x = xsig(:,k).*pnsig;
  sp = fft(x);
  f = round(p.cf/sf*p.drate);
  if p.cf == 0
    ddat(k) = num2str((sign(sp(f+1))+1)/2);
  else
    rad = angle(sp(f+1));
    ddat(k) = num2str(rad < pi/2 | rad > pi*3/2);
  end
end
nx = floor(nframe/7);
x = reshape(ddat(1:7*nx), 7, nx)';
str = char(bin2dec(x))';
str(1:min(length(str),70))

