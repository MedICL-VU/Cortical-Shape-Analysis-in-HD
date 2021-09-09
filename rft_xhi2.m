function [prft] = rft_xhi2(xhi,fwhm,alpha)

Z = linspace(0,7,1000);

resels = prod(size(xhi)/fwhm);

if numel(size(xhi))==2
    expEC = expected_ec2d_z(Z, resels);
else
    expEC = expected_ec3d_z(Z, resels);    
end

tmp = (Z > 2) & (expEC<=alpha);
Z = Z(tmp);
alphaTH = Z(1);
prft = xhi>=alphaTH;

end

% 
% function ec = expected_ec2d_xhi(xhi, reselcount)
% %     ec = (resel * (4 * log(2)) * ((2*pi)**(-3./2)) * z) * exp((z ** 2)*(-0.5));
% end
%     
    
function ec = expected_ec2d_z(z, reselcount)
    ec = reselcount .* ((4 * log(2))/((2*pi)^(3./2.))) .* z .* exp(-(z.^2)/2);
end

function ec = expected_ec3d_z(z,reselcount)
    ec =  reselcount .* (((4*log(2))^(3/2)) / ((2*pi)^2)) .* exp(-(z.^2)/2) .* (z.^2 - 1);
end