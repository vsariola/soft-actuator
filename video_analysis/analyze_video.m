function results = analyze_video(filename,length_pixels,show_results)
    if (nargin < 2)
        length_pixels = 200;
    end

    if (nargin < 3)
        show_results = 0;
    end

    v = VideoReader(filename);
    i = 1;
    N = v.Duration * v.FrameRate;
    results = zeros(N,1);
    millimetersPerPixel = [];
    startingPoint = [];
    
    while hasFrame(v)
        img = readFrame(v);
        
        if isempty(millimetersPerPixel)
            millimetersPerPixel = askMillimetersPerPixel(255 - img);
        end     
        
        if isempty(startingPoint)
            startingPoint = askStartingPoint(255 - img);
        end

        curvature_pixels = find_curvature_pixels(...
            img,...
            startingPoint,....
            mod(i,100) == 0 && showresults);
        results(i) =  * millimetersPerPixel;
        if (mod(i,100) == 0)
            fprintf('Analyzing frame %d / %d (%.2f%% done)\n',i,N,i*100/N);
        end
        i = i + 1;    
    end
end

function ret = find_curvature(img,x0,length_pixels,do_plot)
    theta = linspace(0.001,1.3*pi,1000);
    k = theta./length_pixels;
    x = sin(theta)./k+x0(1);
    y = (1-cos(theta))./k+x0(2);
    g = double(rgb2gray(abs(img2-img1)))/255; 
    c = interp2(g,x,y);
    i = firstPeak(c);
    if (do_plot)
        imshow(img2);
        hold on;
        plot(x,y);
        plot(x(i),y(i),'*');
        hold off;
        drawnow;
    end
    ret = k(i);
    
    function ret = firstPeak(arr)
        sigma = 16;
        kern = ndgauss(3*sigma,sigma,'der',1);
        kern = kern / sqrt(sum(kern.^2));
        filt = conv(arr,kern,'same');
        [py,px] = findpeaks(-filt);
        px = px(py > 0.5);
        ret = px(end);
    end 
end

function ret = askMillimetersPerPixel(img)
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
        title('Click on a point near the left edge of the actuator. The point should be outside the frame area.');
        p = ginput(1);        
        hold on;
        plot(p(:,1),p(:,2),'r*');
        hold off;
        choice = questdlg('Is this starting point OK?','OK?','OK','Pick again','Quit','OK');
        if (strcmp(choice,'OK'))
            ret = p;
            return;
        end
        if (strcmp(choice,'Quit'))
            error('User did not give necessary information for starting point, exiting');
        end
    end       
end

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