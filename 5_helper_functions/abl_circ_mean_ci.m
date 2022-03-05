function [mu,cl,ci1,ci2] = abl_circ_mean_ci(alpha,w,dim,xi,d,deg)

%%%% abl edited from circ_mean_ci from CircStat2012a (on github)

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 5 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin< 4 || isempty(xi)
   xi = .05; 
end

if nargin<6 || isempty(deg)
    deg = false;
end


r = circ_r(alpha,w,d,dim);
% else r = circ_r(alpha);
% end
n = sum(w,dim);
R = n.*r;
ch2 = chi2inv((1-xi),1);
rr = sum(w.*exp(1i*alpha),dim);

if r < .9
    ta = sqrt((2*n*(2*R^2-n*ch2))/(4*n-ch2));  % equ. 26.24
elseif r >= .9
    ta = sqrt(n^2-(n^2-R^2)*exp(ch2/n));      % equ. 26.25
else
    disp(r)
end
t = acos(ta./R);
mu = angle(rr);
ci1 = angle(rr-t); ci2 = angle(rr+t);
cl = angle(t);

if deg
   mu = mod(rad2deg(mu)+360,360);
   ci1 = mod(rad2deg(ci1)+360,360);
   ci2 = mod(rad2deg(ci2)+360,360);
   cl = mod(rad2deg(cl)+360,360);
end