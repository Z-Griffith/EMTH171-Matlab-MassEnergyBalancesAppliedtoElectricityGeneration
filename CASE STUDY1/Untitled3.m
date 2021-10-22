clear
clc
f = @(v) m.*g.*sin(ang).*v + f.*m.*g.*cos(ang).*v -...
    ((b.*(r_g.*v)/r_p)-((r_g.^2.*v.^2)/(r_p.^2)));   
d = @(v) m.*g.*sin(ang) + f.*m.*g.*cos(ang)-((b.*r_g)/r_p)-((2*r_g.^2.*v)/(r_p.^2));

b = f([2])