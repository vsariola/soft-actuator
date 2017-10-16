function retCurvature = analyzeVideo(filename,showplots,skipping)    
    % analyzeVideo  Finds the curvature of a soft fluidic actuator
    %   from a video. The curvature is estimated by using machine vision
    %   to find the pixels on the lower edge of the actuator. To calibrate
    %   the scale (pixels per cm), asks the user to click two points
    %   in the video and input the known distance (in cm) between these
    %   points. Also asks the user to click on a point near the beginning
    %   (left most edge) of the actuator. NOTE: when clicking on the first
    %   point on the actuator edge, click slightly outside the aluminum
    %   frame area. Detecting the edge of the actuator against the black
    %   background is much more reliable.
    % 
    %   curvature = analyzeVideo('/path/to/video.MOV')
    %     Analyzes the video file and returns the curvature.
    %
    %   curvature = analyzeVideo(...,true)
    %     Also plots the pixels found on the lower edge for every 20th
    %     frame analyzed. Useful for debugging purposes; make sure that
    %     the code detects the pixels on the lower edge correctly.
    %
    %   curvature = analyzeVideo(...,true,skipping)
    %     skipping is the number of frames between every analyzed frame.
    %     By default, skipping is 50 i.e. if the video is 50 frames per
    %     second, 1 frame per second is analyzed. Setting skipping to 1
    %     analyzes every frame.

    if (nargin < 2)
        showplots = 0;
    end
    
    if (nargin < 3)
        % The videos are 50 fps, so analyze 1 frame per second
        skipping = 50;
    end

    v = VideoReader(filename);
    i = 1;
    j = 1;
    N = ceil(v.Duration * v.FrameRate / skipping);
    
    retCurvature = zeros(N,1);
    centimetersPerPixel = [];
    startingPoint = [];
    
    tic;
    while hasFrame(v)
        while hasFrame(v) && j > 0
            img = readFrame(v);
            j = j - 1;
        end        
        
        if isempty(centimetersPerPixel)
            centimetersPerPixel = askCentimetersPerPixel(255 - img);
        end     
        
        if isempty(startingPoint)
            startingPoint = askStartingPoint(255 - img);
        end

        nthFrame = mod(i,20) == 1;
        
        curvatureInvMeters = findCurvature(...
            img,...
            startingPoint,...            
            centimetersPerPixel,...
            nthFrame && showplots);
        retCurvature(i) = curvatureInvMeters;
        if (nthFrame || i == N)
            fprintf('Analyzed frame %d / %d (%.2f%% done)\n',i,N,i*100/N);
        end
        j = skipping;
        i = i + 1;    
    end 
    toc
end

function retCurvature = findCurvature(img,x0,centimetersPerPixel,showPlot)
    % Length of the curve to be looked for. Actuator should be somewhat
    % shorter than this, or the edge finding will continue beyond the
    % end of the actuator
    lengthCentimeters = 6;    
    lengthPixels = lengthCentimeters/centimetersPerPixel;
    
    % Number of points to be found along the edge
    numSteps = 10;
    
    % The number of pixels to be checked for each 
    scanWidthPixels = 20;
    
    % Initial guess for curvature, in units of m^-1
    curvatureInvMeters = 5;

    grayImg = double(rgb2gray(img));
    p = findEdgePixels(grayImg,x0,lengthPixels,numSteps,scanWidthPixels);
    
    % https://en.wikipedia.org/wiki/Menger_curvature#Definition
    % and https://math.stackexchange.com/a/516224
    x = p(1:(end-2),:);
    y = p(2:(end-1),:);
    z = p(3:end,:);
    xy = x-y;
    zy = z-y;
    A = (xy(:,2).*zy(:,1)-xy(:,1).*zy(:,2))/2;
    xz = x-z;
    n = @(x) sqrt(sum(x.^2,2));
    c = 4*A ./ (n(xy) .* n(zy) .* n(xz));
    
    retCurvature = mean(c)/centimetersPerPixel*100;
    
    if (showPlot)
       imshow(img);
       hold on;       
       plot(p(:,1),p(:,2),'bx');
       hold off;
       drawnow;
    end   
end

function p = findEdgePixels(grayImg,x0,lengthPixels,numSteps,scanwidthPixels)        
    % alpha is a filtering parameter for the edge detection, [0,1]
    % Smaller values are less sensitive to misdetection but respond to 
    % the curving of the edge slower.
    alpha = 0.5; 
    % sigmaPosPixels is weighting applied to the detected edge. Higher 
    % weight is given to center i.e. an edge that follows straight line
    % from the previous edge    
    sigmaPosPixels = 10; 
    % sigmaDerPixels is the pixel constant of the low pass filter which is
    % applied to the image before taking the first derivative to find 
    % the edge
    sigmaDerPixels = 1;    
    x = x0; % Starting point
    u = [1 0]; % The first scan line is perfectly vertical    
        
    p = zeros(numSteps,2);
    delta = lengthPixels/(numSteps-1);
    scanSteps = scanwidthPixels*2-1; % The edge position is found with half pixel accuracy
    pixelsPerScanStep = scanwidthPixels / scanSteps;
        
    sigmaDerSteps = round(sigmaDerPixels / pixelsPerScanStep);
    % Using a filter longer than 6*sigma or 3*sigma on both sides is
    % unneccesary. Note that ('der',1) parameter performs gaussian low pass
    % filtering and derivative estimation in same filtering step 
    derFilter = ndgauss(round(6*sigmaDerSteps+1),sigmaDerSteps,'der',1); 
        
    sigmaPosSteps = sigmaPosPixels / pixelsPerScanStep;
    weights = ndgauss(scanSteps,sigmaPosSteps)';    
    
    for i = 1:numSteps
        % u is current direction parallel to the edge, v is perpendicular
        v = [-u(2) u(1)];
        p1 = x + v * (scanwidthPixels-1) / 2;
        p2 = x - v * (scanwidthPixels-1) / 2;
        p(i,:) = scanEdgeAlongLine(grayImg,p1,p2,scanSteps,derFilter,weights);
        
        if (i > 1)
            % When we have found at least 2 points, update the 
            % direction which is considered parallel to the edge
            nu = p(i,:) - p(i-1,:);
            nu = nu / norm(nu,2);
            % Filter the change in u, to avoid a single misdetection 
            % ruining everything
            u = (1-alpha) * u + alpha * nu;
            u = u / norm(u,2);        
        end
               
        x = p(i,:) + u * delta;
    end       
end

function p = scanEdgeAlongLine(grayImg,p1,p2,N,derFilter,weights)    
    xx = linspace(p1(1),p2(1),N);
    yy = linspace(p1(2),p2(2),N);
    arr = interp2(grayImg,xx,yy);
    filt = conv(arr,derFilter,'same');
    weighted = filt .* weights;
    [~,ind] = max(weighted);        
    p = [xx(ind) yy(ind)];
end

function ret = askCentimetersPerPixel(img)
    while(1)
        imshow(img);
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        title('Click on a first point on the ruler');
        p1 = ginput(1); 
        hold on;
        plot(p1(1,1),p1(1,2),'r*');
        title('Click on a second point on the ruler');
        p2 = ginput(1);
        plot(p2(1,1),p2(1,2),'r*');
        p = [p1;p2];
        plot(p(:,1),p(:,2),'b-');
        hold off;
        choice = questdlg('Is this scale bar OK?','OK?','OK','Pick again','Quit','OK');
        if (strcmp(choice,'OK'))
            centimeters = inputdlg('And how many centimeters is this?');
            centimeters = str2double(centimeters{1});
            pixels = norm(p(1,:) - p(2,:));
            ret = centimeters / pixels;
            close(gcf);
            return;
        end
        if (strcmp(choice,'Quit'))
            error('User did not give necessary information for computing pixels per centimeter, exiting');
        end
    end       
end

function ret = askStartingPoint(img)
    while(1)
        imshow(img);
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        title('Click on a point on the bottom edge of the actuator, on the side mounted to the aluminum frame. The point should be slightly outside the aluminum frame area.');
        p = ginput(1);        
        hold on;
        plot(p(:,1),p(:,2),'r*');
        hold off;
        choice = questdlg('Is this starting point OK?','OK?','OK','Pick again','Quit','OK');
        if (strcmp(choice,'OK'))
            ret = p;
            close(gcf);
            return;
        end
        if (strcmp(choice,'Quit'))
            error('User did not give necessary information for starting point, exiting');
        end
    end       
end

%% code from ndgauss w/ license

% Copyright (c) 2010, Avan Suinesiaputra
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function h = hermite(n,x)
    % HERMITE: compute the Hermite polynomials.
    % 
    %   h = hermite(n)
    %   h = hermite(n,x)
    % 
    % Inputs:
    %   - n is the order of the Hermite polynomial (n>=0).
    %   - x is (optional) values to be evaluated on the resulting Hermite
    %     polynomial function.
    % 
    % There are two possible outputs:
    % 1. If x is omitted then h is an array with (n+1) elements that contains
    %    coefficients of each Hermite polynomial term.
    %    E.g. calling h = hermite(3)
    %    will result h = [8 0 -12 0], i.e. 8x^3 - 12x
    % 
    % 2. If x is given, then h = Hn(x) and h is the same size of x.
    %    E.g., H2(x) = 4x^2 - 2
    %    calling h = hermite(2,[0 1 2])
    %    will result h = [-2 2 14]
    % 
    % More information:
    % - about the Hermite polynomial: http://mathworld.wolfram.com/HermitePolynomial.html
    % - some examples of this function:
    % http://suinotes.wordpress.com/2010/05/26/hermite-polynomials-with-matlab/
    % 
    % Authors: Avan Suinesiaputra (avan.sp@gmail.com)
    %          Fadillah Z Tala    (fadil.tala@gmail.com)

    % rev.
    % 26/05/2010 - first creation.
    %            - bug fixed: error when hermite(0,x) is called (x isn't empty)
    % 24/09/2010 - bug fixed: the size of x does match with y in line 50.
    %              (thanks to Shiguo Peng)

    % check n
    if( n<0 ), error('The order of Hermite polynomial must be greater than or equal to 0.'); end

    % again check n is an integer
    if( 0~=n-fix(n) ), error('The order of Hermite polynomial must be an integer.'); end

    % call the hermite recursive function.
    h = hermite_rec(n);

    % evaluate the hermite polynomial function, given x
    if( nargin==2 )
        y = h(end) * ones(size(x));
        p = 1;
        for i=length(h)-1:-1:1
            y = y + h(i) * x.^p;
            p = p+1;
        end

        % restore the shape of y, the same as x
        h = reshape(y,size(x));
    end
end

function h = hermite_rec(n)
    % This is the reccurence construction of a Hermite polynomial, i.e.:
    %   H0(x) = 1
    %   H1(x) = 2x
    %   H[n+1](x) = 2x Hn(x) - 2n H[n-1](x)

    if( 0==n ), h = 1;
    elseif( 1==n ), h = [2 0];
    else

        h1 = zeros(1,n+1);
        h1(1:n) = 2*hermite_rec(n-1);

        h2 = zeros(1,n+1);
        h2(3:end) = 2*(n-1)*hermite_rec(n-2);

        h = h1 - h2;

    end
end

function [g,varargout] = ndgauss(hsize,sigma,varargin)
    % NDGAUSS: create ND gaussian kernel in any derivative order.
    %
    %   g = ndgauss(hsize,sigma);
    %   [g,xi,yi,..] = ndgauss(hsize,sigma);
    %
    % Inputs:
    %   - hsize is N-length of kernel size.
    %   - sigma is N-length of standard deviations (gaussian widths).
    %
    % Outputs:
    %   - g is the kernel. The size of g is hsize(1) x hsize(2) x ... hsize(n).
    %   - [xi,yi,..] is the grid to create the kernel.
    %
    % Options:
    % 1. 'der', <gaussian_derivatives>. Default is 0.
    %    The option must be in N-length array. Each defines which derivative
    %    order.
    %    Example:
    %    a. to create d^2G(x,y) / dxdy, then 'der' = [1 1].
    %    b. to create d^2G(x,y) / dy^2, then 'der' = [0 2].
    %
    % 2. 'normalize', 1 | 0. Default is 0.
    %    If this option is set to 1, the kernel will be divided by its sum that
    %    makes the sum of kernel is 1.
    %    Warning: gaussian derivative kernel may not be sum to 1, e.g.,
    %    dG(x,y)/dxdy.
    %
    % Examples:
    % 1. To create 33x33 2D gaussian kernel with different sigma:
    %    zi = ndgauss([33 33],[2*sqrt(2) 1.2]);
    %    imagesc(zi); colormap(gray); axis image;
    %
    % 2. To create several 1D gaussian kernel in different derivative orders:
    %    leg = cell(1,5);
    %    for i=0:4
    %        yi(:,i+1) = ndgauss(33,2*sqrt(2),'der',i);
    %        leg{i+1} = sprintf('g%d(x)',i);
    %    end
    %    plot(yi); 
    %    legend(leg);
    %
    % 3. To create 15x15 Laplacian of Gaussian of width 2:
    %    Lg = ndgauss([15 15],[2 2],'der',[2 0]) + ...
    %         ndgauss([15 15],[2 2],'der',[0 2]);
    %    figure; colormap(gray);
    %    subplot(1,2,1); imagesc(Lg); axis image;
    %    subplot(1,2,2); surf(Lg); axis vis3d;
    %
    % Authors:
    %   Avan Suinesiaputra - avan dot sp at gmail dot com.
    %   Fadillah Z Tala - fadil dot tala at gmail dot com.

    % rev:
    % 27/06/2010 - first creation.

    % default options
    opt.der = [];
    opt.normalize = 0;

    % get options
    for i=1:2:length(varargin)
        if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
        else error('Unknown found.'); end
    end

    % check der & sigma
    if( isempty(opt.der) ), opt.der = zeros(size(hsize)); end

    % check sizes
    dim = length(hsize);
    if( dim ~= length(sigma) ) 
        error('Dimension mismatches.'); 
    end

    % check derivative order
    if( any(opt.der)<0 || any(opt.der-fix(opt.der))~=0 )
        error('Derivative is invalid.'); 
    end

    % check values of hsize & sigma
    if( any(hsize<=0) ), error('One of the kernel size is negative.'); end
    if( any(sigma<=0) ), error('One of the sigma is negative.'); end

    % half kernel size
    sz = (hsize - 1)/2;

    % create N-dimensional grid
    if( 1==dim )
        X = (-sz:sz)';
        varargout = {X};
    else
        sa = ''; sr = ''; T = {};
        for i=1:dim
            sa = sprintf('%s,%f:%f',sa,-sz(i),sz(i));
            sr = sprintf('%s,T{%d}',sr,i);
        end
        eval(sprintf('[%s] = ndgrid(%s);',sr(2:end),sa(2:end)));
        X = zeros(numel(T{1}),dim);
        for i=1:dim
            X(:,i) = reshape(T{i},[],1);
        end
        varargout = {T};
        clear sa sr T;
    end

    % normalized 1D gaussian function
    gfun = @(x,s) exp(-(x.*x)/(2*s.*s)) ./ (sqrt(2*pi)*s);

    % create kernel
    for i=1:dim
        c = sigma(i) * sqrt(2);
        gx = gfun(X(:,i),sigma(i));
        gx(gx<eps*max(X(:,i))) = 0;
        Hn = hermite(opt.der(i),X(:,i) ./ c);
        X(:,i) = Hn .* gx .* (-1/c).^opt.der(i);
    end
    g = prod(X,2);

    % normalize kernel, but derivative kernel may not sum to 1.
    if( opt.normalize )
        sumg = sum(g);
        if( 0~=sumg )
            g = g ./ sumg;
        end
    end

    % reshape kernel
    if( dim>1 )
        g = reshape(g,hsize);
    end
end