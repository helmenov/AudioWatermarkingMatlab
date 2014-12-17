function str = echo_dec(file, p)

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
  file = 'echoenc.wav';
end
if nargin < 2
  p.dt = [100 150];  % deley time of bit0 and bit1 
  p.flen = 8192;     % frame length;  bit rate is sf/p.flen.
%  p.tr = 100;        % transition length
%  p.a = 0.2;         % relative amplitude of an echo
end

[sig, sf] = wavread(file);
len = length(sig);
nframe = floor(len/p.flen);

xsig = reshape(sig(1:p.flen*nframe,1), p.flen, nframe);

for k=1:nframe
  tmp1 = fft(xsig(:,k));
  tmp2 = ifft(log(tmp1));
  ps = abs(fft(abs(tmp2))).^2;
  ac = real(ifft(ps));
  if ac(p.dt(1)+1) >= ac(p.dt(2)+1)
    ddat(k) = '0';
  else
    ddat(k) = '1';
  end
end
nx = floor(nframe/7);
x = reshape(ddat(1:7*nx), 7, nx)';
str = char(bin2dec(x))';
str(1:min(length(str),70))