function leaf(v1,v2,d);
%      a = (x2-x1)^2+(y2-y1)^2+(z2-z1)^2;
%      b = 2*(x1*(x2-x1)+y1*(y2-y1)+z1*(z2-z1));
%      c = x1^2+y1^2+z1^2-d^2;
%      t = (-b+sqrt(b^2-4*a*c))/2/a;
%      x = x1+t*(x2-x1);
%      y = y1+t*(y2-y1);
%      z = z1+t*(z2-z1);

     v = v1+d*(v2-v1)/norm(v2-v1);
     
     %direction of the original branch
     u1 = (v2-v1)/norm(v2-v1);
     %generate random perpendicular direction
     u2(1) = rand-0.5;
     u2(2) = rand-0.5;
     u2(3) = -(u1(1)*u2(1)+u1(2)*u2(2))/u1(3);
     w = 2*u2/norm(u2);
     %line([v(1),v(1)+w(1)],[v(2),v(2)+w(2)],[v(3),v(3)+w(3)],'Color','k','LineWidth',1); 