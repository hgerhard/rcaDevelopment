function plotRCsTopo(A,nComp,symmetricColorbars,alignPolarityToRc1,showColorbar)
% plotRCsTopo(A,nComp,[symmetricColorbars],[alignPolarityToRc1],[showColorbar])
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

if nargin<3
    symmetricColorbars = true;
end
if nargin<4
    alignPolarityToRc1 = true;
end
if nargin<5
    showColorbar = false;
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
    s(1) = 1;
    for rc = 2:nComp
        if f(rc)~=f(1)
            s(rc) = -1;
        end
    end
else
    s = ones(1,nComp);
end

for c=1:nComp
    subplot(1,nComp,c);    
    plotOnEgi(s(c).*A(:,c),colorbarLimits,showColorbar);
    title(['RC' num2str(c)]);
    axis off;
end
