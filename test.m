
alpha=0.05;
iteration=500;
fwhm = 8;
sd = fwhm/sqrt(8*log(2));

radius = 200;

% 2D
fp = zeros(iteration,1);
for i=1:iteration
    rI = imgaussfilt(randn(radius,radius), sd, 'Padding', 'symmetric');
    rI = rI(fwhm:end-fwhm,fwhm:end-fwhm)/std(rI(:));
    
    prft = rft_xhi2(rI, fwhm, alpha);
    fp(i) = double(any(prft(:)));
end
fprintf('False positive rate for 2D : %.3f (expected %.3f)\n',...
    mean(fp), alpha);


% 3D 
fp = zeros(iteration,1);
for i=1:iteration
    rI = imgaussfilt(randn(radius,radius,radius), sd, 'Padding', 'symmetric');
    rI = rI(fwhm:end-fwhm,fwhm:end-fwhm,fwhm:end-fwhm)/std(rI(:));
    
    prft = rft_xhi2(rI, fwhm, alpha);
    fp(i) = double(any(prft(:)));
end
fprintf('False positive rate for 3D : %.3f (expected %.3f)\n',...
    mean(fp), alpha);