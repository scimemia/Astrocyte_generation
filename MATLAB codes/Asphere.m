function f = Asphere(R,W);
   f = 2*pi*R^2-3*R^2*(2*atan(3*R/sqrt(3*R^2-W^2))-W/sqrt(3)/R*asin(sqrt(3)*W/sqrt(12*R^2-W^2)));