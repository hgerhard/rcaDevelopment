function [errorEllipse,amplBounds] = getErrorEllipse(xyData,withinSubj,ellipseType)
%[errorEllipse] = getErrorEllipse(xyData,[withinSubj],[ellipseType],[makePlot])
%   user provides xyData, an Nx2 matrix of 2D data of N samples
%   opt. input withinSubj is a logical specifying whether or not the
%       samples in xyData are within subject. If true, then the two
%       dimensions are pooled for a circular estimate of variance a la
%       Victor & Mast (19XX). Otherwise errors along the two dimensions are
%       considered separately.
%   opt. input ellipseType can be 'SEM' '95CI' or a string specifying
%       a different percentage CI formated following: '%.1fCI'. Default is
%       'SEM'.
%   opt. input makePlot is a logical specifying whether or not to
%       generate a plot of the data & ellipse & eigen vectors (draw at
%       length of 1 std dev, i.e. sqrt(corresponding eigenvalue))
%
%   *NOTE: For EEG measurements in the Fourier domain, xyData rows should
%   be: [real,imag].
%
%   The function uses the eigen value decomposition to find the two
%   perpandicular axes aligned with the data in xyData, where the
%   eigenvalues of each eigen vector correspond to the variances along each
%   direction.
% 
% ### add functionality for withinSubj error bars a la Victor & Mast

if nargin<2 || isempty(withinSubj)
    withinSubj = false;
end
if nargin<3 || isempty(ellipseType)
    ellipseType = 'SEM';
end
if withinSubj
    warning('Code does not yet support within Subject error bars that pool errors over real \& imag components.\n');
    withinSubj = false;
end

dims = size(xyData);
N = dims(1);
if dims(2) ~= 2
    error('input data must be a matrix of 2D row samples');
end
if N < 2
    error('input data must contain at least 2 samples');
end

srIx = 1;
siIx = 2;

[eigenvec, eigenval] = eig(cov([xyData(:,srIx),xyData(:,siIx)]));
[orderedVals,eigAscendIx] = sort(diag(eigenval));
smallest_eigenvec = eigenvec(:,eigAscendIx(1));
largest_eigenvec = eigenvec(:,eigAscendIx(2));
smallest_eigenval = orderedVals(1);
largest_eigenval = orderedVals(2);
phi = atan2(largest_eigenvec(2), largest_eigenvec(1));
% This angle is between -pi and pi, shift to 0 and 2pi:
if(phi < 0)
    phi = phi + 2*pi;
end

theta_grid = linspace(0,2*pi);

switch ellipseType
    case '1STD'
        a = sqrt(largest_eigenval); 
        b = sqrt(smallest_eigenval);
    case '2STD'
        a = 2*sqrt(largest_eigenval); 
        b = 2*sqrt(smallest_eigenval);
    case 'SEMarea'
        a = sqrt(largest_eigenval/sqrt(N)); % scales the ellipse's area by sqrt(N)
        b = sqrt(smallest_eigenval/sqrt(N));
    case {'SEMellipse' 'SEM'} % ### DEFAULT!
        a = sqrt(largest_eigenval)/sqrt(N); % contour at stdDev/sqrt(N)
        b = sqrt(smallest_eigenval)/sqrt(N);
    case '95CI'
        if ~withinSubj
            Tsq = ((N-1)*2)/(N-2) * finv( 0.95, 2, N - 2 );
            a = sqrt(largest_eigenval*Tsq);
            b = sqrt(smallest_eigenval*Tsq);
        else
            % ### Victor & Mast version..
        end
    otherwise
        if strcmp(ellipseType(end-1:end),'CI')
          if ~withinSubj
              critVal = str2double(ellipseType(1:end-2))./100;
              if critVal < 1 && critVal > 0
                  Tsq = ((N-1)*2)/(N-2) * finv( critVal, 2, N - 2 );
                  a = sqrt(largest_eigenval*Tsq);
                  b = sqrt(smallest_eigenval*Tsq);
              else
                  error('CI range must be on the interval (0, 100). Please see the help!')
              end
          else
              % ### Victor & Mast version..
          end
        else
            error('You entered an invalid error ellipse type. Please see the help!')
        end
end

% the ellipse in x and y coordinates
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
errorEllipse = [ellipse_x_r;ellipse_y_r]' * R;

%Shift to be centered on mean coordinate
meanXy = mean(xyData);
meanXyAmp = norm(meanXy);
errorEllipse = bsxfun(@plus,errorEllipse,meanXy);

% find vector lengths of each point on the ellipse
norms = nan(1,length(errorEllipse));
for pt = 1:length(errorEllipse)
    norms(pt) = norm(errorEllipse(pt,:));
end
[minNorm,minNormIx] = min(norms);
[maxNorm,maxNormIx] = max(norms);
ellipseExtremes = [minNorm,maxNorm];

if (sign(max(errorEllipse(:,1))) ~= sign(min(errorEllipse(:,1)))) && (sign(max(errorEllipse(:,2))) ~= sign(min(errorEllipse(:,2))))
    % the ellipse overlaps with the origin
    amplBounds = [0, maxNorm];
else
    amplBounds = ellipseExtremes;
end













