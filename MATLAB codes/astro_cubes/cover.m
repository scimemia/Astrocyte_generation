function CVbranch = cover(u,v,LFT,RGT,stp,PT)
    
%     LFT = -10;
%     RGT = 10;
%     stp = 0.1;
%     PT = floor((RGT-LFT)/stp);
    CVbranch = zeros(PT,PT,PT);
     
    x = LFT:stp:RGT;
    y = LFT:stp:RGT;
    z = LFT:stp:RGT;

     a = floor((u(1)-LFT)/stp);
     b = floor((v(1)-LFT)/stp);
     lftx = max(0,min(a,b));
     rgtx = min(PT,max(a,b));
     for i=lftx+1:rgtx+1;
         if (i==lftx+1); y1 = u(2); else; y1 = u(2) + (x(i)-u(1))/(v(1)-u(1))*(v(2)-u(2));
         end
         if(i==rgtx+1); y2 = v(2); else; y2 = u(2) + (x(i+1)-u(1))/(v(1)-u(1))*(v(2)-u(2));
         end
         a = floor((y1-LFT)/stp);
         b = floor((y2-LFT)/stp);
         lfty = max(0,min(a,b));
         rgty = min(PT,max(a,b));
         for j = lfty+1:rgty+1;
             if (i==lftx+1); z1 = u(3); else;
                 if (j==lfty+1); z1 = u(3) + (x(i)-u(1))/(v(1)-u(1))*(v(3)-u(3));
                 else; z1 = u(3) + (y(j)-u(2))/(v(2)-u(2))*(v(3)-u(3)); 
                 end
             end
             if (i==rgtx+1); z2 = v(3); else;
                 if (j==rgty+1); z2 = u(3) + (x(i+1)-u(1))/(v(1)-u(1))*(v(3)-u(3));
                 else; z2 = u(3) + (y(j+1)-u(2))/(v(2)-u(2))*(v(3)-u(3));
                 end
             end

             a = floor((z1-LFT)/stp);
             b = floor((z2-LFT)/stp);
             lftz = max(0,min(a,b));
             rgtz = min(PT,max(a,b));
             for k = lftz+1:rgtz+1;
               CVbranch(i,j,k) = 1;
             end
         end;
     end;
