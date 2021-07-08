function primary_branch(U)
global W r0 R c c0 d D N M Csholl Nsholl step h0 h1 Hcount NT SA ST VOL LGTH brch Ntip Stip SAtip
global G_sec L_sec G_prim L_prim xi lambda klength count f col td;

%basic segment length for this specific primary branch;
c0 = c*rand;

%initialize all measures to zero for primary branch; 
%these are to be added to the cell totals in the main code;

%number of branches per level;
G_prim = zeros(N+1,1);

%total length per level;
L_prim = zeros(N+1,1);

%surface area;
SA = zeros(M+N,1);
SAtip = 0;

%volume;
VOL = 0;

%total branch length;
LGTH = 0;

%number of branch points;
brch = 0;

%number of transporters (depends on cytoplasmic insertion depth);
NT = zeros(M+N,Hcount+1);
Ntip = zeros(1,Hcount+1);

%transporter area (depends on cytoplasmic insertion depth);
ST = zeros(M+N,Hcount+1);
Stip = zeros(1,Hcount+1);

%choose the tip diameter minimum (threshold) value for the primary branch,
%as a random number in the range tm to tM;
tM = 0.15;
tm = 0.02;
td = rand*(tM-tm)+tm;

%initiate current length of primary branch (for calculating coordinates);
s = 0;

xi = zeros(M,1);
lambda = zeros(M,1);
count = 0;

%for each added segment;
for j = 1:M;
%only continue adding segments as long as the smaller child branch 
%passes the tip diameter min threshold;
    if (D(j)/((1+f^(3/2))^(2/3))>td)

        %segments have length between c0 and c0+R;
        csi = c0+R*rand;

        %define the start and the end of each segment;
        u = (r0+s)*U; 
        v = (r0+s+csi)*U;

        %update the total length with the primary branch (for computing coordinates);
        s = s + csi;

        %plot the segment of the primary branch;
        line([v(1),u(1)],[v(2),u(2)],[v(3),u(3)],'Color',col(1,:),'LineWidth',1.5);

        %add the branch to the sholl count;
              for sect = 1:1:Nsholl-1;
                  Rsholl = r0 + sect*step;
                  if (u(1)^2+u(2)^2+u(3)^2 < Rsholl^2) && (v(1)^2+v(2)^2+v(3)^2 > Rsholl^2);
                      Csholl(sect) = Csholl(sect)+1;
                  end;
              end;   
              
        if (j<M);      
            %add the number of transporters in the segment to the total
            for hi = 1:Hcount+1
                h = h0 + (h1-h0)*(hi-1)/Hcount;
                NT(1,hi) = NT(1,hi) + ncyl(D(j),csi,h,W);
                %add the transporter area of the segment to the prim br
                ST(1,hi) = ST(1,hi) + ncyl(D(j),csi,h,W)*Acyl1(D(j),W);
            end

            %add the surface area of the segment to the total
            SA(1) = SA(1) + D(j)*pi*csi;
            %add the volume of the segment to the total       
            VOL = VOL + pi*csi*D(j)^2/4;
        
        else
            %add the number of transporters in the segment to the total
            for hi = 1:Hcount+1
                h = h0 + (h1-h0)*(hi-1)/Hcount;
                NT(1,hi) = NT(1,hi) + ncyl(D(j),csi-D(j)/2,h,W);
                Ntip(hi) = Ntip(hi) + nhsphere(D(j)/2,h,W);
                %add the transporter area of the segment to the prim br
                ST(1,hi) = ST(1,hi) + ncyl(D(j),csi-D(j)/2,h,W)*Acyl1(D(j),W);
                Stip(hi) = Stip(hi) + nhsphere(D(j)/2,h,W)*Asphere(D(j)/2,W)
            end

            %add the surface area of the segment to the total
            SA(1) = SA(1) + D(j)*pi*(csi-D(j)/2);
            SAtip = SAtip + 2*pi*D(j)^2/4
            %add the volume of the segment to the total       
            VOL = VOL + pi*(csi-D(j)/2)*D(j)^2/4 + 2*pi*D(j)^3/8/3;
        end;
        
        %add the length of the segment to the total
        LGTH = LGTH + csi;
        %add the length of the segment to the prim br
        L_prim(1) = L_prim(1) + csi;

       %call the construction of secondary branch that starts at level j;
       secondary_branch(u,v,j);  
       
       
       %when finished, add the number of branches off the secondary branch and their length to the
       %respective level vector (level 1 on sec branch is actually level 2);
       G_prim(2:N+1) = G_prim(2:N+1) + G_sec;
       L_prim(2:N+1) = L_prim(2:N+1) + L_sec;
       count = count + 1;
       xi(j) = csi;
       lambda(j) = klength;
    else
        return
    end
end
