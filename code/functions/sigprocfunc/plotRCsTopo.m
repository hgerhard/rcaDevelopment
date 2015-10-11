function plotRCsTopo(A,nComp,symmetricColorbars,alignPolarityToRc1,showColorbar,forceOzPositive)
% plotRCsTopo(A,nComp,[symmetricColorbars],[alignPolarityToRc1],[showColorbar],[forceOzPositive])
% 
% make a set of topographic maps of the A, weight matrices, for each RC. In
% this plot, the color range used is the same for each topographic map so
% they can be compared easily. The color range is fixed so that zero values
% take on the midpoint of the color range.
%
% if desired, the polarity of RC2 .. RCN are aligned with RC1 such that the
% color used to plot the most extreme direction of RC1 (+ or -) is the same
% color used to plot the most extreme direction of each of the following
% RCs

if nargin<3 || isempty(symmetricColorbars)
    symmetricColorbars = true;
end
if nargin<4 || isempty(alignPolarityToRc1)
    alignPolarityToRc1 = true;
end
if nargin<5 || isempty(showColorbar)
    showColorbar = false;
end
if nargin<6 || isempty(forceOzPositive)
    forceOzPositive = false;
end

if symmetricColorbars
    % for a consistent colorbar across RCs:
    colorbarLimits = [min(A(:)),max(A(:))];
    newExtreme = max(abs(colorbarLimits));
    colorbarLimits = [-newExtreme,newExtreme];
else 
    colorbarLimits = [];
end

figure;

if alignPolarityToRc1
    extremeVals = [min(A); max(A)];
    for rc = 1:nComp
        [~,f(rc)]=max(abs(extremeVals(:,rc)));
    end
    s = ones(1,nComp);
    if f(1)==1 % largest extreme value of RC1 is negative
        s(1) = -1; % flip the sign of RC1 so its largest extreme value is positive
        f(1) = 2;
    end
    if forceOzPositive && (s(1)*A(75,1))<0
        fprintf('The extremes of A(:,1) are %1.3f and %1.3f. You are flipping the sign so Oz will be positive.\n',extremeVals(1,1),extremeVals(2,1));
        s(1) = -(s(1));
        if f(1) == 1
            f(1) = 2;
        else
            f(1) = 1;
        end
    end
        
    for rc = 2:nComp
        if f(rc)~=f(1)
            s(rc) = -1;
        end
    end
else
    s = ones(1,nComp);
end

if forceOzPositive && (s(1)*A(75,1))<0    
    extremeVals = [min(A); max(A)];
    fprintf('The extremes of A(:,1) are %1.3f and %1.3f. You are flipping the sign so Oz will be positive.\n',extremeVals(1,1),extremeVals(2,1));
    s(1) = -(s(1));
end

for c=1:nComp
    subplot(1,nComp,c);    
    plotOnEgi(s(c).*A(:,c),colorbarLimits,showColorbar);
    title(['RC' num2str(c)]);
    axis off;
end

% make extra plot of colorbar for figure making purposes later ###