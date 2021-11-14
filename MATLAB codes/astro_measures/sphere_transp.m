number = zeros(10,1);
area = zeros(10,1);
r = 5;
d1 = 5;
w = 0.008;

for i=1:11;
h = 0.003 + (i-1)*0.00035;
number(i) = 4*pi/(3*acos((2*r^2+2.*h.^2-4*h*sqrt(r^2-w^2/3)-w^2)/(4*r^2+4*h.^2-8*h*sqrt(r^2-w^2/3)-w^2))-pi);
Adelta(i) = 2*pi*r^2-3*r^2*(2*atan(3*r/sqrt(3*r^2-w^2))-w/sqrt(3)/r*asin(sqrt(3)*w/sqrt(12*r^2-w^2)));
Ats(i) = number(i)*Adelta(i);
end;

Acap = 2*pi*r*(r-sqrt(r^2-d1^2/4));

Scaled_number = number * (4*pi*r^2 - 3*Acap)/(4*pi*r^2);
Scaled_area = Ats * (4*pi*r^2 - 3*Acap)/(4*pi*r^2);