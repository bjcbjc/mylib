function [p, h, stats] = signrankplus(x,y,varargin)
%SIGNRANK Wilcoxon signed rank test for zero median.
%   P = SIGNRANK(X) performs a two-sided signed rank test of the hypothesis
%   that the data in the vector X come from a distribution whose median
%   (and mean, if it exists) is zero, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("median is zero") is true.
%   Small values of P cast doubt on the validity of the null hypothesis.
%   The data are assumed to come from a continuous distribution, symmetric
%   about its median.
%
%   P = SIGNRANK(X,M) performs a two-sided test of the hypothesis that the
%   data in the vector X come from a distribution whose median is M.  M
%   must be a scalar.
%
%   P = SIGNRANK(X,Y) performs a paired, two-sided test of the hypothesis
%   that the difference between the matched samples in the vectors X and Y
%   comes from a distribution whose median is zero.  The differences X-Y
%   are assumed to come from a continuous distribution, symmetric about its
%   median.  X and Y must be the same length.  The two-sided p-value is
%   computed by doubling the most significant one-sided value.
%
%   SIGNRANK treats NaNs in X or Y as missing values, and removes them.
%
%   [P,H] = SIGNRANK(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("median is zero") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = SIGNRANK(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = SIGNRANK(...,'method',METHOD) computes the p-value using an
%   exact algorithm if METHOD is 'exact', or a normal approximation if
%   METHOD is 'approximate'.  The default is to use an exact method for
%   small samples.
%
%   [P,H,STATS] = SIGNRANK(...) returns STATS, a structure with one or two
%   fields.  The field 'signedrank' contains the value of the signed rank
%   statistic.  If P is calculated using a normal approximation, then the
%   field 'zval' contains the value of the normal (Z) statistic.
%
%   See also SIGNTEST, RANKSUM, TTEST, ZTEST.

%   For the two-sample case, SIGNRANK uses a tolerance based on the
%   values EPSD=EPS(X)+EPS(Y). Any pair of values of D=X-Y that differ by
%   no more than the sum of their two EPSD values are treated as ties.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.


%   Copyright 1993-2011 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2011/05/09 01:26:51 $

% Check most of the inputs now
alpha = 0.05;
epsdiff = [];
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  signrank(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
oknames = {'alpha' 'method' 'epsdiff'};
dflts   = {alpha   '' epsdiff};
[alpha,method, epsdiff] = internal.stats.parseArgs(oknames,dflts,varargin{:});

if ~isscalar(alpha)
   error(message('stats:signrank:BadAlpha'));
end
if ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error(message('stats:signrank:BadAlpha'));
end

onesample = false;
if nargin < 2 || isempty(y)
    y = zeros(size(x));
    onesample = true;
elseif isscalar(y)
    y = repmat(y, size(x));
end

if ~isvector(x) || ~isvector(y)
    error(message('stats:signrank:InvalidData'));
elseif numel(x) ~= numel(y)
    error(message('stats:signrank:InputSizeMismatch'));
end

diffxy = x(:) - y(:);
if isempty(epsdiff)
    if onesample
        epsdiff = zeros(size(x(:)));
    else
        epsdiff = eps(x(:)) + eps(y(:));
    end
else
    epsdiff = repmat(epsdiff, size(x(:)), 1);
end
% Remove missing data
t = isnan(diffxy);
diffxy(t) = [];
epsdiff(t) = [];
if isempty(diffxy)
   %error(message('stats:signrank:NotEnoughData'));
   p = NaN;
   h = NaN;
   stats = [];
   return
end

t = (abs(diffxy) < epsdiff);
diffxy(t) = [];
epsdiff(t) = [];

n = length(diffxy);

if (n == 0)         % degenerate case, all ties
    p = 1;
    if (nargout > 1)
        h = 0;
        if (nargout > 2)
            stats.signedrank = 0;
        end
    end
    return
end

% Now deal with the method argument
if isempty(method)
   if n<=15
      method = 'exact';
   else
      method = 'approximate';
   end
elseif isequal(lower(method),'oldexact')
    % OK
else
   method = internal.stats.getParamVal(method,{'exact' 'approximate'},'''method''');
end

% Find negative differences and ranks of absolute differences
neg = (diffxy<0);
[tierank, tieadj] = tiedrank(abs(diffxy),0,0,epsdiff);

% Compute signed rank statistic (most extreme version)
w = sum(tierank(neg));
if w < n*(n+1)/2-w
    stats.signdirection = -1;
else
    stats.signdirection = +1;
end
w = min(w, n*(n+1)/2-w);

if isequal(method,'approximate')
    z = (w-n*(n+1)/4) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
    p = 2*normcdf(z,0,1);
    if (nargout > 2)
        stats.zval = z;
    end
elseif isequal(method,'oldexact')
    % Enumerates all possibilities and does not adjust for ties
    allposs = (ff2n(n))';
    idx = (1:n)';
    idx = idx(:,ones(2.^n,1));
    pranks = sum(allposs.*idx,1);
    tail = 2*length(find(pranks <= w)); % two side.

    % Avoid p>1 if w is in the middle and is double-counted
    p = min(1, tail./(2.^n));
else % isequal(method,'exact')
    p = statsrexact(tierank,w);
    p = min(1, 2*p);   % two-sided, don't double-count the middle value
end

if nargout > 1
    h = (p<=alpha);
    if (nargout > 2)
        stats.signedrank = w;
    end
end

end


function [pval,P] = statsrexact(v,w)
%STATSREXACT Compute exact tail probability for signed rank statistic.
%   [PVAL,ALLP]=STATSREXACT(V,W) computes the tail probability PVAL
%   for the statistic W with the vector V of ranks.  ALLP is a matrix
%   containing the probabilities (col. 2) for each W value (col. 1).
%
%   Private function used by the SIGNRANK function.

%   Copyright 2003-2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:30:30 $

n = length(v);
v = sort(v(:)'); % make sure it's a row

% For convenience we can just compute the lower tail.  If w is
% in the upper tail, compute its equivalent lower tail value.
maxw = n*(n+1)/2;
folded = (w>maxw/2);
if folded
   w = maxw-w;
end

% We would like to use the elements of w and v as indexes into
% arrays that enumerate possible values.  If there are ties causing
% non-integer ranks, multiply by 2 to force everything to integer.
doubled = any(v~=floor(v));
if doubled
   v = round(2*v);
   w = round(2*w);
end

C = zeros(w+1,1);  % C(w+1) will be the number of combinations adding
                   % to w at each step
C(1) = 1;          % just one combination includes nothing
top = 1;           % top entry currently in use in C vector

% Look at all combinations of ranks that could contribute
% to the observed value of W
for vj=v(v<=w)

   % C now enumerates combinations not including v(j).  Now update the
   % elements that could include v(j).
   newtop = min(top+vj,w+1);
   hi = min(vj,w+1)+1:newtop;
   lo = 1:length(hi);

   C(hi) = C(hi) + C(lo);

   top = newtop;
end

% Convert to probabilities
C = C / (2^n);

% Get tail probability
pval = sum(C);

if nargout>1
   allw = 0:w;
   if doubled
      allw = allw/2;
   end
   if folded
      allw = n*(n+1)/2 - allw;
   end
   
   P = [allw(:), C(:)];
end
end
