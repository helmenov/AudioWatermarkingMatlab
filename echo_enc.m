function [out, encbit] = echo_enc(file, embstr, p);

% Echo hiding 
%   embed embstr string to two different delay times
% Type
%   echo_enc;
% on the command line of MATLAB or Octave 2.1.x.
%    Host signal is 'wmtest.wav'.
%    Stego signal is written in 'echoenc.wav'. 

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
  file = 'wmtest.wav';
end

if nargin < 2
  embstr = 'JASJ';
end

if nargin < 3
  p.dt = [100 150];  % deley time of bit0 and bit1 
  p.flen = 8192;     % frame length;  bit rate is sf/p.flen.
  p.tr = 200;        % transition length
  p.a = 0.6;         % relative amplitude of an echo 
end

[sig, sf] = wavread(file);
[len, ch] = size(sig);

embnum = double(embstr);  % convert char to int
bitstr = dec2bin(embnum);  % convert int(dec) to int(bin)
[n, cn] = size(bitstr);      
b = reshape(bitstr', cn*n, 1);
rept = floor((len/p.flen)/(n*cn));
encbit = reshape(str2num(b)*ones(1, rept), rept*n*cn, 1);  % repeat
                                                           % embedding
                                                           % data
dtsig = (cos(pi*(1:p.tr)'/p.tr)+1)/2;  % transition envelope
utsig = flipud(dtsig);

dsig0 = [zeros(p.dt(1), ch); sig]*p.a;   % delayed signal of bit0 
dsig1 = [zeros(p.dt(2), ch); sig]*p.a;   % delayed signal of bit1

mlen = p.flen*rept*n*cn;
bsig = reshape(ones(p.flen, 1)*encbit', mlen, 1);
dseg = diff(encbit);
useg = find(dseg > 0);
dseg = find(dseg < 0);

for k=1:length(useg)    % transition from bit0 to bit1
  bsig(useg(k)*p.flen-p.tr/2+1:useg(k)*p.flen+p.tr/2) = utsig;
end
for k=1:length(dseg)    % transition from bit1 to bit0
  bsig(dseg(k)*p.flen-p.tr/2+1:dseg(k)*p.flen+p.tr/2) = dtsig;
end
bsig = bsig*ones(1,ch); 

% echo hiding
out = sig(1:mlen,:)+dsig0(1:mlen,:).*abs(bsig-1)+dsig1(1:mlen,:).*bsig;

if len > mlen
  out = [out; sig(mlen+1:len,:)];  % padding original signal
end

if nargout == 0
  wavwrite(out, sf, 'echoenc.wav');
end
