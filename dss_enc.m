function [out, encbit] = dss_enc(file, embstr, p);

% Data Hiding by Direct Spread Spectram 
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
  p.trate = 4;         % tip rate
  p.drate = 8192;      % data rate
  p.cf = 0;            % carrier freq.
%  p.cf = 44100/4;     % carrier freq.
  p.lev = -20;         % embedding level relative to the signal level
  p.key = 1374;        % key of pseudo-random sequence; seed of randomize
  p.seg = 512;         % length of level calculation
end

[sig, sf] = wavread(file);
[len, ch] = size(sig);

embnum = double(embstr);  % convert char to int
bitstr = dec2bin(embnum);  % convert int(dec) to int(bin)
[n, cn] = size(bitstr);      
b = reshape(bitstr', cn*n, 1);
rept = ceil((len/p.drate)/(n*cn));
encbit = reshape(str2num(b)*ones(1, rept), rept*n*cn, 1)*2-1;  % repeat
                                                           % embedding
                                                           % data
rand('seed', p.key);           % set seed of random sequence
a = rand(p.drate/p.trate, 1);
pnseq = (a > 0.5)*2-1;        % pseudo-random sequence consisting of +1 and -1 
pnsig = reshape(ones(p.trate, 1)*pnseq', p.drate, 1);

fnum = floor(len/p.drate);     % number of data frame
if p.cf == 0
  car = ones(p.drate*fnum, 1);   % carrier signal
else
  car = sin(2*pi*p.cf/sf*(1:p.drate*fnum)');
end

for k=1:ch
  fsig = reshape(sig(1:fnum*p.drate,k), p.drate, fnum);
  csig = reshape(car, p.drate, fnum);
  xsig = zeros(p.drate, fnum);
  for kk=1:fnum
    % power level of the signal
    slev = 10*log10(mean(reshape(fsig(:,kk), p.seg, p.drate/p.seg).^2)+eps);
    % power level of the carrier
    clev = 10*log10(mean(reshape(csig(:,kk), p.seg, p.drate/p.seg).^2)+eps);
    % calculate embedding power
    lev = max(slev-clev+p.lev, -60);
    % reshaping
    xpnsig = reshape(pnsig, p.seg, p.drate/p.seg).*(ones(p.seg,1)*10.^(lev/20))*encbit(kk);
    % add DSS signal of BPSK 
    xsig(:,kk) = fsig(:,kk) + csig(:,kk).*reshape(xpnsig, p.drate, 1);
  end
  out(:,k) = reshape(xsig, p.drate*fnum, 1);
end

if len > p.drate*fnum   % padding original signal
  out = [out; sig(p.drate*fnum+1:len,:)];
end

if nargout == 0
  wavwrite(out, sf, 'dssenc.wav');
end
