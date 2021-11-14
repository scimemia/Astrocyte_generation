function secondary_branch_cube_cover(u,v,level)
global R c0 d N G_sec L_sec tipdiam
global SA VOL LGTH brch SAtip VOLtip;
global D diam col ct td klength;
global freq_leaves nr_leaves lh fh cth Diamhr Vh Vhtip SAh SAhtip;
global LFT RGT stp PT;
global CV1 CV2 Volcore;

%initialize vectors with number of branches and branch length / generation;
G_sec = zeros(N,1);
L_sec = zeros(N,1);

%initialize all diameters to zero, for all possible child branches;
diam = zeros(2^N,1);
%initialize the diameter of the starting secondary branch to what was
%determined in the main code;
diam(1) = d(level);

%construct secondary segment;
%this segment is level+1 away from the soma;
        %all secondary segments make an angle theta1 with primary brance;
        theta1 = pi/6;
        %pick a random 3D direction for the segment;
        V = rand(3,1)-0.5;
        V = V'/norm(V);
        w = (v-u)/norm(v-u);
        dir = cos((-1)^level*theta1)*w + sin((-1)^level*theta1)*V; 
        %length of the segment is from c0 to c0+R;
        csi = c0+R*rand;
        %add the segment length to the total length of the level;
        klength = csi;
        %define the end point of the segment;
        x = v + csi*dir/norm(dir);
        %plot the secondary segment;
        line([v(1),x(1)],[v(2),x(2)],[v(3),x(3)],'Color',col(2,:),'LineWidth',1.5);
        
        %add the branch point to the total number;
        brch = brch + 1;
        %add branch to the volume cover;
        if (level==1)
          CV1 = CV1 + cover(v,x,LFT,RGT,stp,PT);
          Volcore = Volcore + pi*(diam(1))^2/4*csi;
        else
          CV2 = CV2 + cover(v,x,LFT,RGT,stp,PT);
        end
          
            
          %surface area and length (for the total cell measures);
          SA(2) = SA(2) + diam(1)*pi*csi;
          LGTH = LGTH + csi; 
          %volume (for the total cell measures);
          VOL(2) = VOL(2) + pi*csi*diam(1)^2/4;
          
          %adjust the volume to include the branch point transition;
          hbp = diam(1)*tan(theta1);
          Delta = (D(level)^1.5+diam(1)^1.5)^(2/3);
          VOL(1) = VOL(1) - pi*hbp*diam(1)^2/4 + 1/3*pi*hbp*(D(level)^2+Delta^2+D(level)*Delta);

%%%%%%%%%%%%%%%%%%%%%
%CONTINUED BRANCHING
%%%%%%%%%%%%%%%%%%%%%
%prepare vectors for the N subsequent branch levels (max of 2^N segments);
%the start and the end of each segment are recorded in the vectors of coordinates s and t;
    s = zeros(2^N-1,3); t = zeros(2^N-1,3); 
    %length of each segment;
    leng = zeros(2^N,1); 
    %generation of each segment;
    gen = zeros(2^N,1);
    
    %the start and end of the secondary segment are already defined;
    s(1,:) = u; t(1,:) = v; 
    %the length of the first segment is csi;
    leng(1) = csi; 
    %the generation of the secondary segment is 1
    %(will get shifted by one when added to the primary branch profile)
    gen(1) = 1;
    
%i indexes all binary branches sprouting from this secondary branch
for i = 1:2^N-1
%fraction bias in the 3/2 rule (random for each cell between 1 and 1.5);
f = 1+0.5*rand;
            %Compute the number of leaves based on the segment length and
            %inter-leaf distance, and call "leaf" to draw 
            csi = norm(t(i,:)-s(i,:));
            if (s(i,1)~=1i)
            nr_leaves = floor(csi/freq_leaves);
            Dh = diam(i);
            fh = 0.5 + 3*(1+rand)*Dh;
            for hr = 1:nr_leaves
                leaf(s(i,:),t(i,:),hr*freq_leaves);
                diamhr = Dh/((1+fh^(3/2))^(2/3));
                Diamhr(cth) = diamhr;
                cth = cth+1;
                Dh = fh*Dh/((1+fh^(3/2))^(2/3)); 
                    SAh = SAh + pi*diamhr*(lh-diamhr/2);
                    SAhtip = SAhtip + 2*pi*diamhr^2/4;
                    Vh = Vh + pi*(lh-diamhr/2)*diamhr^2/4;
                    Vhtip = Vhtip + 2*pi*diamhr^3/8/3;

            end
            end
        diam(i) = Dh;
        
        %the angle theta increases with the generation
        theta = pi/6 + gen(i)*pi/80;
        %the termination condition: tip diameter of the smaller child > td;
        if ((diam(i)/((1+f^(3/2))^(2/3))> td) && (s(i,1)~=1i))
        %pick up the start and end points of the mother branch in the
        %vectors u and v;
        u = s(i,:); v = t(i,:);
        
        %add the branch point to the total count;
        brch = brch + 1;

         
        %BRANCH SEGMENT 
        %in a random direction, angle theta with the mother branch direction, and length csi;
            V = rand(3,1)-0.5;
            V = V'/norm(V);
            w = (v-u)/norm(v-u);
            dir = cos((-1)^i*theta)*w + sin((-1)^i*theta)*V;
            csi = c0+R*rand;
            xa = v + csi*dir/norm(dir);
            %compute the diameter of new segment;
            diam(2*i+1) = diam(i)/((1+f^(3/2))^(2/3));
            %assign the generation to the new segment;
            gen(2*i+1) = gen(i)+1;

            %plot the branch segment;
            line([v(1),xa(1)],[v(2),xa(2)],[v(3),xa(3)],'Color',col(gen(2*i+1),:),'LineWidth',1.5); 
            
            %add branch to cover;
            CV2 = CV2 + cover(v,xa,LFT,RGT,stp,PT);

            %if this is not a terminal segment
            if (diam(2*i+1)/((1+f^(3/2))^(2/3))>td)
                SA(gen(2*i+1)+1) = SA(gen(2*i+1)+1) + diam(2*i+1)*pi*csi;
                VOL(gen(2*i+1)+1) = VOL(gen(2*i+1)+1) + pi*csi*diam(2*i+1)^2/4;
                LGTH = LGTH + csi;
                
                %adjust the volume to include the branch point transition;
                hbp = diam(2*i+1)*tan(theta);
                Delta = (diam(i)^1.5+diam(2*i+1)^1.5)^(2/3);
                VOL(gen(i)) = VOL(gen(i)) - pi*hbp*diam(2*i+1)^2/4 + 1/3*pi*hbp*(diam(i)^2+Delta^2+diam(i)*Delta);
                                  
            %if it is a terminal segment    
            else
            ct = ct + 1;
            tipdiam(ct) = diam(2*i+1);
            %add the segment measures to the total counts;
                SA(gen(2*i+1)+1) = SA(gen(2*i+1)+1) + diam(2*i+1)*pi*(csi-diam(2*i+1)/2);
                SAtip = SAtip + 2*pi*diam(2*i+1)^2/4;
                VOL(gen(2*i+1)+1) = VOL(gen(2*i+1)+1) + pi*(csi-diam(2*i+1)/2)*diam(2*i+1)^2/4;
                VOLtip = VOLtip + 2*pi*diam(2*i+1)^3/8/3;
                LGTH = LGTH + csi;
                
                %adjust the volume to include the branch point transition;
                VOL(gen(i)) = VOL(gen(i)) - pi*hbp*diam(2*i+1)^2/4 + 1/3*pi*hbp*(diam(i)^2+Delta^2+diam(i)*Delta);
            end
                 
                    
        %BRANCH SEGMENT
        %in the same direction as the mother branch;     
            xb = v + csi*w;
            %calculate diameter of new segment, with bias f;    
            diam(2*i) = f*diam(i)/((1+f^(3/2))^(2/3));
            %assign the level (same as the parent);
            gen(2*i) = gen(i);
            
            %plot the segment;
            line([v(1),xb(1)],[v(2),xb(2)],[v(3),xb(3)],'Color',col(gen(2*i),:),'LineWidth',1.5);
            
            %add branch to cover;
            CV2 = CV2 + cover(v,xb,LFT,RGT,stp,PT);
                    
            %if this is not a terminal branch    
            if (diam(2*i)/((1+f^(3/2))^(2/3))>td)
                        %add the segment measures to the total counts;
                        SA(gen(2*i)+1) = SA(gen(2*i)+1) + diam(2*i)*pi*csi;
                        VOL(gen(2*i)+1) = VOL(gen(2*i)+1) + pi*csi*diam(2*i)^2/4;
                        LGTH = LGTH + 1;
                        

             %if this is a terminal segment
            else
            ct = ct + 1;
            tipdiam(ct) = diam(2*i);
                        %add the segment measures to the total counts;
                        SA(gen(2*i)+1) = SA(gen(2*i)+1) + diam(2*i)*pi*(csi-diam(2*i)/2);
                        SAtip = SAtip + 2*pi*diam(2*i)^2/4;
                        VOL(gen(2*i)+1) = VOL(gen(2*i)+1) + pi*(csi-diam(2*i)/2)*diam(2*i)^2/4;
                        VOLtip = VOLtip + 2*pi*diam(2*i)^3/8/3;
                        LGTH = LGTH + 1;
            end

            %record the starts and ends of the two new branches in the
            %vectors s and t;
            s(2*i,:) = v; s(2*i+1,:) = v; 
            t(2*i,:) = xb; t(2*i+1,:) = xa;
            %set the mother generation to zero (since we no longer need that)
            gen(i) = 0;
            %update the lengths of the branches;
            leng(2*i) = leng(i)+csi; leng(2*i+1) = csi;
                
        else 
        %if branching does not occur, give the points complex values, to also prevent future branching;
        %the tip ciam vaue is added to the vectore tipdiam;
         s(2*i,1) = 1i; s(2*i+1,1) = 1i;
        end     
end

%for each level, search for how many branches G_sec(k) are in each level k 
%sum up their lengths, to get the total length for each level L_sec(k);
%recall that levels are shifted: level 1 here is actually level 2.
for k=1:N
    S = size(find(gen==k));
    G_sec(k) = S(1);
        if S(1)>0
        L_sec(k) = sum(leng(find(gen==k)));
        end
end

        