function f = ncyl(D,L,H,W);
   f = 2*pi/acot(-2*H/W + sqrt(D^2/W^2-1)) * 2*L/sqrt(3)/W;