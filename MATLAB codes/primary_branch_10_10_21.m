function primary_branch_10_10_21(U)
global W r0 R c c0 D N M Csholl Nsholl step h0 h1 Hcount tipdiam
global NT SA SA_th ST VOL LGTH brch Ntip Stip SAtip VOLtip 
global G_sec L_sec G_prim L_prim 
global D xi lambda klength count col td ct thold
global freq_leaves nr_leaves lh fh cth Diamhr Vh Vhtip SAh SAhtip NTh NTtiph STh STtiph

%basic segment length for this specific primary branch;
c0 = c*rand;

%fraction bias in the 3/2 rule (random for each cell between 1 and 1.5);
f = 1+0.5*rand;

%initialize all measures to zero for primary branch; 
%these are to be added to the cell totals in the main code;

%number of branches per level;
G_prim = zeros(N+1,1);
G_prim_th = zeros(N+1,1);

%total length per level;
L_prim = zeros(N+1,1);
L_prim_th = zeros(N+1,1);

%surface area;
SA = zeros(M+N,1);
SAtip = 0;
SA_th = zeros(M+N,1);
SAtip_th = 0;

%volume;
VOL = 0;
VOLtip = 0;
VOL_th = 0;
VOLtip_th = 0;

%total branch length;
LGTH = 0;
LGTH_th = 0;

%number of branch points;
brch = 0;
brch_th = 0;

%number of transporters (depends on cytoplasmic insertion depth);
NT = zeros(M+N,Hcount+1);
Ntip = zeros(1,Hcount+1);
NTh = zeros(1,Hcount+1);
NTtiph = zeros(1,Hcount+1);

%transporter area (depends on cytoplasmic insertion depth);
ST = zeros(M+N,Hcount+1);
Stip = zeros(1,Hcount+1);
STh = zeros(1,Hcount+1);
STtiph = zeros(1,Hcount+1);

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

j = 1;
G_prim(1) = 0;
G_prim_th(1) = 0;
%for each added segment;
while j < M+1
        %segments have length between c0 and c0+R;
        csi = c0+R*rand;

        %define the start and the end of each segment;
        u = (r0+s)*U; 
        v = (r0+s+csi)*U;

        %update the total length with the primary branch (for computing coordinates);
        s = s + csi;

        %plot the segment of the primary branch;
        line([v(1),u(1)],[v(2),u(2)],[v(3),u(3)],'Color',col(1,:),'LineWidth',1.5);
        
        %add the branch point to the total number;
        brch = brch + 1;
        if (D(j)>thold) && (csi>thold)
        brch_th = brch_th + 1;
        %add the branch to the sholl count;
              for sect = 1:1:Nsholl-1
                  Rsholl = r0 + sect*step;
                  if (u(1)^2+u(2)^2+u(3)^2 < Rsholl^2) && (v(1)^2+v(2)^2+v(3)^2 > Rsholl^2)
                      Csholl(sect) = Csholl(sect)+1;
                  end
              end 
        end
        
          %HERE, COMPUTE NUMBER OF LEAVES BASED ON THE LENGTH OF THE SEGMENT 
          %AND ON THE DISTANCE BETWEEN THEM
          %DRAW THEM ON THE SEGMENT (JUST FOR ILLUSTRATION)                  
            nr_leaves = floor(csi/freq_leaves);
            Dh = D(j);
            fh = 0.5 + 20*(1+rand)*Dh;
            for hr = 1:nr_leaves
                leaf(u,v,hr*freq_leaves);
                diamhr = Dh/((1+fh^(3/2))^(2/3));
                Diamhr(cth) = diamhr;
                cth = cth+1;
                Dh = fh*Dh/((1+fh^(3/2))^(2/3)); 
                    SAh = SAh + pi*diamhr*lh
                    SAhtip = SAhtip + 2*pi*diamhr^2/4;
                    Vh = Vh + pi*(lh-diamhr/2)*diamhr^2/4;
                    Vhtip = Vhtip + 2*pi*diamhr^3/8/3;
                
                if diamhr > 3*h1
                for hi=1:Hcount+1
                            h = h0 + (h1-h0)*(hi-1)/Hcount;
                            hout = h1 - h;
                            NTh(hi) = NTh(hi) + ncyl(diamhr,lh-hout,h,W);
                            NTtiph(hi) = NTtiph(hi) + nhsphere(diamhr/2,h,W);
                            STh(hi) = STh(hi) + ncyl(diamhr,lh-hout,h,W)*Acyl1(diamhr,W);
                            STtiph(hi) = STtiph(hi) + nhsphere(diamhr/2,h,W)*Asphere(diamhr/2,W);
                end
                STh(1)
                end
            end
                
        if (j<M)
                %add the surface area of the segment to the total
                SA(1) = SA(1) + D(j)*pi*csi;
                %add the volume of the segment to the total       
                VOL = VOL + pi*csi*D(j)^2/4;        
                %add the length of the segment to the total
                LGTH = LGTH + csi;
                %add the length of the segment to the prim br
                L_prim(1) = L_prim(1) + csi;
                   
            %add the number of transporters in the segment to the total
            for hi = 1:Hcount+1
                h = h0 + (h1-h0)*(hi-1)/Hcount;
                NT(1,hi) = NT(1,hi) + ncyl(D(j),csi,h,W);
                %add the transporter area of the segment to the prim br
                ST(1,hi) = ST(1,hi) + ncyl(D(j),csi,h,W)*Acyl1(D(j),W);
            end
        
        %for terminal branches, gotta add the tips
        else
                ct = ct + 1;
                tipdiam(ct) = D(j) + 1;
                %add the surface area of the segment to the total
                SA(1) = SA(1) + D(j)*pi*(csi-D(j)/2);
                SAtip = SAtip + 2*pi*D(j)^2/4;
                %add the volume of the segment to the total       
                VOL = VOL + pi*(csi-D(j)/2)*D(j)^2/4;
                VOLtip = VOLtip + 2*pi*D(j)^3/8/3;
                %add the length of the segment to the total
                LGTH = LGTH + csi;

                %add the number of transporters in the segment to the total
                for hi = 1:Hcount+1
                    h = h0 + (h1-h0)*(hi-1)/Hcount;
                    NT(1,hi) = NT(1,hi) + ncyl(D(j),csi-D(j)/2,h,W);
                    Ntip(hi) = Ntip(hi) + nhsphere(D(j)/2,h,W);
                    %add the transporter area of the segment to the prim br
                    ST(1,hi) = ST(1,hi) + ncyl(D(j),csi-D(j)/2,h,W)*Acyl1(D(j),W);
                    Stip(hi) = Stip(hi) + nhsphere(D(j)/2,h,W)*Asphere(D(j)/2,W);
                end
        end

       %call the construction of secondary branch that starts at level j;
       secondary_branch_10_10_21(u,v,j);  
       
       
       %when finished, add the number of branches off the secondary branch and their length to the
       %respective level vector (level 1 on sec branch is actually level 2);
       G_prim(2:N+1) = G_prim(2:N+1) + G_sec;
       L_prim(2:N+1) = L_prim(2:N+1) + L_sec;
       count = count + 1;
       xi(j) = csi;
       lambda(j) = klength;  
    D(j) = Dh; 
    
    %only continue adding segments as long as the smaller child branch 
    %passes the tip diameter min threshold;
    if (D(j)/((1+f^(3/2))^(2/3))>td)
        j = j+1;
    else; j = M+1;
    end
end