function [p,anovatab,stats] = anova1_figurecontrol(x,group,displayopt,extra)
%ANOVA1 One-way analysis of variance (ANOVA).
%   ANOVA1 performs a one-way ANOVA for comparing the means of two or more 
%   groups of data. It returns the p-value for the null hypothesis that the
%   means of the groups are equal.
%
%   P = ANOVA1(X,GROUP,DISPLAYOPT)
%   If X is a matrix, ANOVA1 treats each column as a separate group, and
%     determines whether the population means of the columns are equal.
%     This form of ANOVA1 is appropriate when each group has the same
%     number of elements (balanced ANOVA).  GROUP can be a character
%     array or a cell array of strings, with one row per column of
%     X, containing the group names.  Enter an empty array ([]) or
%     omit this argument if you do not want to specify group names.
%   If X is a vector, GROUP must be a categorical variable, vector,
%     string array, or cell array of strings with one group name for
%     each element of X.  X values corresponding to the same value of
%     GROUP are placed in the same group.
%
%   DISPLAYOPT can be 'on' (the default) to display figures
%   containing a standard one-way anova table and a boxplot, or
%   'off' to omit these displays.  Note that the notches in the
%   boxplot provide a test of group medians (see HELP BOXPLOT),
%   and this is not the same as the F test for different means
%   in the anova table.
%
%   [P,ANOVATAB] = ANOVA1(...) returns the ANOVA table values as the
%   cell array ANOVATAB.
%
%   [P,ANOVATAB,STATS] = ANOVA1(...) returns an additional structure
%   of statistics useful for performing a multiple comparison of means
%   with the MULTCOMPARE function.
%
%   See also ANOVA2, ANOVAN, BOXPLOT, MANOVA1, MULTCOMPARE.

%   Reference: Robert V. Hogg, and Johannes Ledolter, Engineering Statistics
%   Macmillan 1987 pp. 205-206.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.6 $  $Date: 2011/08/29 20:38:09 $



narginchk(1,4);

classical = 1;
nargs = nargin;
if (nargin>0 && strcmp(x,'kruskalwallis'))
   % Called via kruskalwallis function, adjust inputs
   classical = 0;
   if (nargin >= 2), x = group; group = []; end
   if (nargin >= 3), group = displayopt; displayopt = []; end
   if (nargin >= 4), displayopt = extra; end
   nargs = nargs-1;
end

if (nargs < 2), group = []; end
if (nargs < 3), displayopt = 'on'; end
% Note: for backwards compatibility, accept 'nodisplay' for 'off'
if isempty(displayopt)
    willdisplay = true;    
elseif isnumeric(displayopt)
    if displayopt > 0
        willdisplay = true;
    else
        willdisplay = false;
    end
else
    willdisplay = ~(strcmp(displayopt,'nodisplay') | strcmp(displayopt,'n') ...
                | strcmp(displayopt,'off'));
end
% Convert group to cell array from character array, make it a column
if (ischar(group) && ~isempty(group)), group = cellstr(group); end
if (size(group, 1) == 1), group = group'; end

% If X is a matrix with NaNs, convert to vector form.
if ~isvector(x)
   if (any(isnan(x(:))))
      [n,m] = size(x);
      x = x(:);
      gi = reshape(repmat((1:m), n, 1), n*m, 1);
      if isempty(group)     % no group names
         group = gi;
      elseif (size(group,1) == m)
         group = group(gi,:);
      else
         error(message('stats:anova1:InputSizeMismatch'));
      end
   end
elseif ~isempty(group) && (size(group,1) ~= length(x))
   error(message('stats:anova1:InputSizeMismatch'));
end

% If X is a matrix and GROUP is strings, use GROUPs as names
if (iscell(group) && (length(x) < numel(x)) ...
                  && (size(x,2) == size(group,1)))
   named = 1;
   gnames = group;
   grouped = 0;
else
   named = 0;
   gnames = [];
   grouped = ~isempty(group);
end

if (grouped)
   % Single data vector and a separate grouping variable
   x = x(:);
   lx = length(x);
   if (lx ~= numel(x))
      error(message('stats:anova1:VectorRequired'))
   end
   nonan = ~isnan(x);
   x = x(nonan);

   % Convert group to indices 1,...,g and separate names
   group = group(nonan,:);
   [groupnum, gnames] = grp2idx(group);
   named = 1;

   % Remove NaN values
   nonan = ~isnan(groupnum);
   if (~all(nonan))
      groupnum = groupnum(nonan);
      x = x(nonan);
   end

   lx = length(x);
   xorig = x;                    % use uncentered version to make M
   groupnum = groupnum(:);
   maxi = size(gnames, 1);
   if isa(x,'single')
      xm = zeros(1,maxi,'single');
   else
      xm = zeros(1,maxi);
   end
   countx = xm;
   if (classical)
      mu = mean(x);
      x = x - mu;                % center to improve accuracy
      xr = x;
   else
      [xr,tieadj] = tiedrank(x);
   end
   
   for j = 1:maxi
      % Get group sizes and means
      k = find(groupnum == j);
      lk = length(k);
      countx(j) = lk;
      xm(j) = mean(xr(k));       % column means
   end

   gm = mean(xr);                      % grand mean
   df1 = sum(countx>0) - 1;            % Column degrees of freedom
   df2 = lx - df1 - 1;                 % Error degrees of freedom
   xc = xm - gm;                       % centered
   xc(countx==0) = 0;
   RSS = dot(countx, xc.^2);           % Regression Sum of Squares
else
   % Data in matrix form, no separate grouping variable
   [r,c] = size(x);
   lx = r * c;
   if (classical)
      xr = x;
      mu = mean(xr(:));
      xr = xr - mu;           % center to improve accuracy
   else
      [xr,tieadj] = tiedrank(x(:));
      xr = reshape(xr, size(x));
   end
   countx = repmat(r, 1, c);
   xorig = x;                 % save uncentered version for boxplot
   xm = mean(xr);             % column means
   gm = mean(xm);             % grand mean
   df1 = c-1;                 % Column degrees of freedom
   df2 = c*(r-1);             % Error degrees of freedom
   RSS = r*(xm - gm)*(xm-gm)';        % Regression Sum of Squares
end

TSS = (xr(:) - gm)'*(xr(:) - gm);  % Total Sum of Squares
SSE = TSS - RSS;                   % Error Sum of Squares

if (df2 > 0)
   mse = SSE/df2;
else
   mse = NaN;
end

if (classical)
   if (SSE~=0)
      F = (RSS/df1) / mse;
      p = fpval(F,df1,df2);        % Probability of F given equal means.
   elseif (RSS==0)                 % Constant Matrix case.
      F = NaN;
      p = NaN;
   else                            % Perfect fit case.
      F = Inf;
      p = 0;
   end
else
   F = (12 * RSS) / (lx * (lx+1));
   if (tieadj > 0)
      F = F / (1 - 2 * tieadj/(lx^3-lx));
   end
   p = chi2pval(F,df1);
end


Table=zeros(3,5);               %Formatting for ANOVA Table printout
Table(:,1)=[ RSS SSE TSS]';
Table(:,2)=[df1 df2 df1+df2]';
Table(:,3)=[ RSS/df1 mse Inf ]';
Table(:,4)=[ F Inf Inf ]';
Table(:,5)=[ p Inf Inf ]';

colheads = {getString(message('stats:anova1:ColHeadSource')), getString(message('stats:anova1:ColHeadSS')), getString(message('stats:anova1:ColHeadDf')), getString(message('stats:anova1:ColHeadMS')), getString(message('stats:anova1:ColHeadF')), getString(message('stats:anova1:ColHeadProbGtF'))};

if (~classical)
   colheads{5} = getString(message('stats:anova1:ColHeadChisq'));
   colheads{6} = getString(message('stats:anova1:ColHeadProbGtChisq'));
end
rowheads = {getString(message('stats:anova1:RowHeadColumns')), getString(message('stats:anova1:RowHeadError')), getString(message('stats:anova1:RowHeadTotal'))};
if (grouped)
   rowheads{1} = getString(message('stats:anova1:RowHeadGroups'));
end

% Create cell array version of table
atab = num2cell(Table);
for i=1:size(atab,1)
   for j=1:size(atab,2)
      if (isinf(atab{i,j}))
         atab{i,j} = [];
      end
   end
end
atab = [rowheads' atab];
atab = [colheads; atab];
if (nargout > 1)
   anovatab = atab;
end

% Create output stats structure if requested, used by MULTCOMPARE
if (nargout > 2)
   if ~isempty(gnames)
      stats.gnames = gnames;
   else
      stats.gnames = strjust(num2str((1:length(xm))'),'left');
   end
   stats.n = countx;
   if (classical)
      stats.source = 'anova1';
      stats.means = xm + mu;
      stats.df = df2;
      stats.s = sqrt(mse);
   else
      stats.source = 'kruskalwallis';
      stats.meanranks = xm;
      stats.sumt = 2 * tieadj;
   end
end


if (~willdisplay), return; end

% digits = [-1 -1 0 -1 2 4];
% if (classical)
%    wtitle = getString(message('stats:anova1:OnewayANOVA'));
%    ttitle = getString(message('stats:anova1:ANOVATable'));
% else
%    wtitle = getString(message('stats:anova1:KruskalWallisOnewayANOVA'));
%    ttitle = getString(message('stats:anova1:KruskalWallisANOVATable'));
% end
% tblfig = statdisptable(atab, wtitle, ttitle, '', digits);
% set(tblfig,'tag','table');

%f1 = figure('pos',get(gcf,'pos') + [0,-200,0,0],'tag','boxplot');
if ~isempty(displayopt)
    %figure(displayopt);
    set(0,'CurrentFigure',displayopt);
end
%ax = axes('Parent',f1);
%ax = get(f1, 'currentaxes');
ax = gca;
if (~grouped)
   boxplot(ax,xorig,'notch','on');
else
   boxplot(ax,xorig,groupnum,'notch','on');
   h = get(ax,'XLabel');
   set(h,'String',getString(message('stats:anova1:GroupNumber')));
end

% If there are group names, use them
if ~isempty(gnames)
   h = get(ax,'XLabel');
   if (named)
      set(h,'String','');
   end
   set(ax, 'xtick', (1:size(gnames,1)), 'xticklabel', gnames);
end
end


function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2010/11/08 02:37:46 $

if nargin < 3, 
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
end

function p = chi2pval(x,v)
%FPVAL Chi-square distribution p-value function.
%   P = CHI2PVAL(X,V) returns the upper tail of the chi-square cumulative
%   distribution function with V degrees of freedom at the values in X.  If X
%   is the observed value of a chi-square test statistic, then P is its
%   p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also CHI2CDF, CHI2INV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.

%   Copyright 2009 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:29:08 $

if nargin < 2
    error(message('stats:chi2pval:TooFewInputs'));
end

[errorcode,x,v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:chi2pval:InputSizeMismatch'));
end

% Return NaN for out of range parameters.
v(v <= 0) = NaN;
x(x < 0) = 0;

p = gammainc(x/2,v/2,'upper');

end
