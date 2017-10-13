pixelsPerM = (750 - 653)/0.03;

img1 = imread('bg.png');
v = VideoReader('temp/MVI_4346.MOV');
i = 1;
N = v.Duration * v.FrameRate;
k = zeros(N,1);

while hasFrame(v)
    img2 = readFrame(v);
    
    k(i) = testi(img1,img2,[765 302],[1068 346],mod(i,100) == 0) * pixelsPerM;
    if (mod(i,100) == 0)
        fprintf('%d / %d\n',i,N);
    end
    i = i +1;    
end