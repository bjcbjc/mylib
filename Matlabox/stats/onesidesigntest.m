function [p, h, stats] = onesidesigntest(x,y,side,varargin)
%SIGNTEST Sign test for zero median.
%   P = SIGNTEST(X) performs a two-sided sign test of the hypothesis that
%   the data in the vector X come from a distribution whose median is zero,
%   and returns the p-value from the test.  P is the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis ("median is zero") is true.  Small values of P cast doubt
%   on the validity of the null hypothesis.  The data are assumed to come
%   from an arbitrary continuous distribution.  SIGNTEST omits values where
%   X is zero or NaN.
%
%   P = SIGNTEST(X,M) performs a two-sided test of the hypothesis that the
%   data in the vector X come from a distribution whose median is M.  M
%   must be a scalar.  SIGNTEST omits values where X is M or NaN.
%
%   P = SIGNTEST(X,Y) performs a paired, two-sided test of the hypothesis
%   that the difference between the matched samples in the vectors X and
%   Y comes from a distribution whose median is zero.  The differences X-Y
%   are assumed to come from an arbitrary continuous distribution.  X and Y
%   must be the same length.  The two-sided p-value is computed by doubling
%   the most significant one-sided value.  SIGNTEST omits values where X-Y
%   is 0 or NaN.
%
%   [P,H] = SIGNTEST(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("median is zero") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = SIGNTEST(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = SIGNTEST(...,'method',METHOD) computes the p-value using an
%   exact algorithm if METHOD is 'exact', or a normal approximation if
%   METHOD is 'approximate'.  The default is to use an exact method for
%   small samples.
%
%   [P,H,STATS] = SIGNTEST(...) returns STATS, a structure with two fields
%   'sign' and 'zval'. The field 'sign' contains the value of the sign
%   statistic. If P is calculated using a normal approximation, then the
%   field 'zval' contains the value of the normal (Z) statistic, otherwise, 
%   'zval' contains NaN.
%
%   See also SIGNRANK, RANKSUM, TTEST, ZTEST.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.

%   Copyright 1993-2009 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:26:37 $

% Check most of the inputs now
alpha = 0.05;
if nargin>3 && isnumeric(varargin{1})
   % Grandfathered syntax:  signtest(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
oknames = {'alpha' 'method'};
dflts   = {alpha   ''};
[eid,emsg,alpha,method] = internal.stats.getargs(oknames,dflts,varargin{:});
if ~isempty(eid)
   error(sprintf('stats:signtest:%s',eid),emsg);
end

if ~isscalar(alpha)
   error('stats:signtest:BadAlpha','SIGNRANK requires a scalar ALPHA value.');
end
if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error('stats:signtest:BadAlpha','SIGNRANK requires 0 < ALPHA < 1.');
end

if isempty(y)
    y = zeros(size(x));
elseif isscalar(y)
    y = repmat(y, size(x));
end

if ~isvector(x) || ~isvector(y)
    error(message('stats:signtest:InvalidData'));
elseif numel(x) ~= numel(y)
    error(message('stats:signtest:InputSizeMismatch'));
end

diffxy = x(:) - y(:);
nodiff = (diffxy == 0);
diffxy(nodiff | isnan(diffxy)) = [];
n = numel(diffxy);

if n == 0 % this means the two vectors are identical
    p = 1;
    h = 0;
    stats.sign = 0;
    stats.zval = NaN;
    return
end

% Now deal with the method argument
if isempty(method)
   if n<100
      method = 'exact';
   else
      method = 'approximate';
   end
elseif ischar(method)
   okmethods = {'exact' 'approximate'};
   j = strmatch(lower(method),okmethods);
   if isempty(j)
      error(message('stats:signrank:BadMethod'));
   end
   method = okmethods{j};
else
   error(message('stats:signrank:BadMethod'));
end

npos = length(find(diffxy>0));
nneg = n-npos;

%sgn = min(nneg,npos);
if strcmp(side, '-')
    sgn = nneg;
else
    sgn = npos;
end

%if isequal(method,'exact')
    %p = min(1, binocdf(sgn,n,0.5));  % p>1 means center value double-counted
    p = min(1, sum(binocdf(sgn:n, n, 0.5)));
    stats.zval = NaN;
%else
    % Do a continuity correction, keeping in mind the right direction
%     z = (npos-nneg - sign(npos-nneg))/sqrt(n);
%     p = 2*normcdf(-abs(z),0,1);
%     stats.zval = z;
%     
%end

 h = (p<=alpha);
 stats.sign = sgn;
    

