function primary_branch_cube_cover(U)
global r0 R c c0 D N M tipdiam;
global SA VOL LGTH brch SAtip VOLtip ;
global G_sec L_sec G_prim L_prim;
global xi lambda klength count col td ct;
global freq_leaves nr_leaves lh fh cth Diamhr Vh Vhtip SAh SAhtip;
global LFT RGT stp PT;
global CV1 CV2 Volcore;

%basic segment length for this specific primary branch;
c0 = c*rand;

%fraction bias in the 3/2 rule (random for each cell between 1 and 1.5);
f = 1+0.5*rand;

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
VOL = zeros(M+N,1);
VOLtip = 0;

%total branch length;
LGTH = 0;

%number of branch points;
brch = 0;

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

        %add branch to cover;
        if (j<3)
          CV1 = CV1 + cover(u,v,LFT,RGT,stp,PT);
          Volcore = Volcore + pi*(D(j))^2/4*csi;
        else
          CV2 = CV2 + cover(u,v,LFT,RGT,stp,PT);
        end
        

          %HERE, COMPUTE NUMBER OF LEAVES BASED ON THE LENGTH OF THE SEGMENT 
          %AND ON THE DISTANCE BETWEEN THEM
          %DRAW THEM ON THE SEGMENT (JUST FOR ILLUSTRATION)                  
            nr_leaves = floor(csi/freq_leaves);
            Dh = D(j);
            fh = 0.5 + 3*(1+rand)*Dh;
            for hr = 1:nr_leaves
                leaf(u,v,hr*freq_leaves);
                diamhr = Dh/((1+fh^(3/2))^(2/3));
                Diamhr(cth) = diamhr;
                cth = cth+1;
                Dh = fh*Dh/((1+fh^(3/2))^(2/3)); 
                    SAh = SAh + pi*diamhr*(lh-diamhr/2);
                    SAhtip = SAhtip + 2*pi*diamhr^2/4;
                    Vh = Vh + pi*(lh-diamhr/2)*diamhr^2/4;
                    Vhtip = Vhtip + 2*pi*diamhr^3/8/3;
            end
                
        if (j<M)
                %add the surface area of the segment to the total
                SA(1) = SA(1) + D(j)*pi*csi;
                %add the volume of the segment to the total       
                VOL(1) = VOL(1) + pi*csi*D(j)^2/4;        
                %add the length of the segment to the total
                LGTH = LGTH + csi;
                %add the length of the segment to the prim br
                L_prim(1) = L_prim(1) + csi;
                   
        %for terminal branches, gotta add the tips
        else
                ct = ct + 1;
                tipdiam(ct) = D(j) + 1;
                %add the surface area of the segment to the total
                SA(1) = SA(1) + D(j)*pi*(csi-D(j)/2);
                SAtip = SAtip + 2*pi*D(j)^2/4;
                %add the volume of the segment to the total       
                VOL(1) = VOL(1) + pi*(csi-D(j)/2)*D(j)^2/4;
                VOLtip = VOLtip + 2*pi*D(j)^3/8/3;
                %add the length of the segment to the total
                LGTH = LGTH + csi;
        end

       %call the construction of secondary branch that starts at level j;
       secondary_branch_cube_cover(u,v,j);  
       
       
       %when finished, add the number of branches off the dary branch and their length to the
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