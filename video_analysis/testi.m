function ret = testi(img1,img2,x0,x1,doplot)
    l = norm(x1-x0);
    theta = linspace(0.001,1.3*pi,1000);
    k = theta./l;
    x = sin(theta)./k+x0(1);
    y = (1-cos(theta))./k+x0(2);
    g = double(rgb2gray(abs(img2-img1)))/255; 
    c = interp2(g,x,y);
    i = firstPeak(c);
    if (doplot)
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