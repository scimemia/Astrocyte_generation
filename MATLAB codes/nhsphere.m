function f = nhsphere(R,H,W);
   f = 2*pi/(3*acos((2*R^2+2*H^2-4*H*sqrt(R^2-W^2/3)-W^2)/(4*R^2+4*H^2-8*H*sqrt(R^2-W^2/3)-W^2))-pi);