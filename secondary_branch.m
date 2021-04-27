function secondary_branch(u,v,level)
global theta W r0 R c0 N Csholl Nsholl step G_sec L_sec gen Ntip Stip SAtip;
global Hcount h0 h1 NT SA ST VOL LGTH brch klength d diam f col ct tipdiam td;

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
        %define the end point of the segment;
        x = v + csi*dir;
        %plot the secondary segment;
        line([v(1),x(1)],[v(2),x(2)],[v(3),x(3)],'Color',col(2,:),'LineWidth',1.5);
        %add the branch point to the total number;
        brch = brch+1;
        %add the segment length to the total length of the level;
        klength = csi;
        %add the segment to the sholl profile;
          for sect = 1:1:Nsholl-1
              Rsholl = r0 + sect*step;
              if (v(1)^2+v(2)^2+v(3)^2 < Rsholl^2) && (x(1)^2+x(2)^2+x(3)^2 > Rsholl^2)
                  Csholl(sect) = Csholl(sect)+1;
              end
          end
          %calculate measures for the new segment;
                  %transporter number and transporter area (depend on insertion depth)
                  for hi=1:Hcount+1
                      h = h0 + (h1-h0)*(hi-1)/Hcount;
                      NT(2,hi) = NT(2,hi) + ncyl(diam(1),csi,h,W);
                      ST(2,hi) = ST(2,hi) + ncyl(diam(1),csi,h,W)*Acyl1(diam(1),W);
                  end
                  %surface area, volume and length (for the total cell measures);
                  SA(2) = SA(2) + diam(1)*pi*csi;
                  VOL = VOL + pi*csi*diam(1)^2/4;
                  LGTH = LGTH + csi;
      
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
        %the angle theta increases with the generation
        theta = pi/6 + gen(i)*pi/80;
        %the termination condition: tip diameter of the smaller child > td;
        if ((diam(i)/((1+f^(3/2))^(2/3))> td) && (s(i,1)~=1i))
        %if branching occurs
            %add the branch point to the total count;
            brch = brch+1;
            %pick up the start and end points of the mother branch in the
            %vectors u and v;
            u = s(i,:); v = t(i,:);
            %construct a branch segment in a random direction, angle theta with the
            %mother branch direction, and length csi;
            V = rand(3,1)-0.5;
            V = V'/norm(V);
            w = (v-u)/norm(v-u);
            dir = cos((-1)^i*theta)*w + sin((-1)^i*theta)*V;    
            csi = c0+R*rand;
            xa = v + csi*dir;
            %compute the diameter of new segment;
            diam(2*i+1) = diam(i)/((1+f^(3/2))^(2/3));
            %assign the generation to the new segment;
            gen(2*i+1) = gen(i)+1;
            %plot the branch segment;
            line([v(1),xa(1)],[v(2),xa(2)],[v(3),xa(3)],'Color',col(gen(2*i+1),:),'LineWidth',1.5); 
            %add the branch segment to the sholl profile;
                    for sect = 1:1:Nsholl-1
                        Rsholl = r0 + sect*step;
                        if (v(1)^2+v(2)^2+v(3)^2 < Rsholl^2) && (xa(1)^2+xa(2)^2+xa(3)^2 > Rsholl^2)
                            Csholl(sect) = Csholl(sect)+1;
                        end
                    end
                    
                    if (diam(2*i+1)/((1+f^(3/2))^(2/3))>td)
                        %add the segment measures to the total counts;
                        for hi=1:Hcount+1
                            h = h0 + (h1-h0)*(hi-1)/Hcount;
                            NT(gen(2*i+1)+1,hi) = NT(gen(2*i+1)+1,hi) + ncyl(diam(2*i+1),csi,h,W);
                            ST(gen(2*i+1)+1,hi) = ST(gen(2*i+1)+1,hi) + ncyl(diam(2*i+1),csi,h,W)*Acyl1(diam(2*i+1),W);
                        end
                        SA(gen(2*i+1)+1) = SA(gen(2*i+1)+1) + diam(2*i+1)*pi*csi;
                        VOL = VOL + pi*csi*diam(2*i+1)^2/4;
                    else
                        %add the segment measures to the total counts;
                        for hi=1:Hcount+1
                            h = h0 + (h1-h0)*(hi-1)/Hcount;
                            NT(gen(2*i+1)+1,hi) = NT(gen(2*i+1)+1,hi) + ncyl(diam(2*i+1),csi-diam(2*i+1)/2,h,W);
                            Ntip(hi) = Ntip(hi) + nhsphere(diam(2*i+1)/2,h,W);
                            ST(gen(2*i+1)+1,hi) = ST(gen(2*i+1)+1,hi) + ncyl(diam(2*i+1),csi-diam(2*i+1)/2,h,W)*Acyl1(diam(2*i+1),W);
                            Stip(hi) = Stip(hi) + nhsphere(diam(2*i+1)/2,h,W)*Asphere(diam(2*i+1)/2,W);
                        end
                        SA(gen(2*i+1)+1) = SA(gen(2*i+1)+1) + diam(2*i+1)*pi*(csi-diam(2*i+1)/2);
                        SAtip = SAtip + 2*pi*diam(2*i+1)^2/4;
                        VOL = VOL + pi*(csi-diam(2*i+1)/2)*diam(2*i+1)^2/4 + 2*pi*diam(2*i+1)^3/8/3;
                    end;
                    
                    LGTH = LGTH + csi;
                    

            %construct a branch segment in the same direction as the mother branch;     
            xb = v + csi*w;
            %calculate diameter of new segment, with bias f;
            diam(2*i) = f*diam(i)/((1+f^(3/2))^(2/3));
            %assign the level (same as the parent);
            gen(2*i) = gen(i);
            %plot the segment;
            line([v(1),xb(1)],[v(2),xb(2)],[v(3),xb(3)],'Color',col(gen(2*i),:),'LineWidth',1.5);
            %add the segment to the sholl profile;
                    for sect = 1:1:Nsholl-1
                        Rsholl = r0 + sect*step;
                        if (v(1)^2+v(2)^2+v(3)^2 < Rsholl^2) && (xb(1)^2+xb(2)^2+xb(3)^2 > Rsholl^2)
                            Csholl(sect) = Csholl(sect)+1;
                        end
                    end
                    
                    if (diam(2*i)/((1+f^(3/2))^(2/3))>td)
                        %add the segment measures to the total counts;
                        for hi=1:Hcount+1
                            h = h0 + (h1-h0)*(hi-1)/Hcount;
                            NT(gen(2*i)+1,hi) = NT(gen(2*i)+1,hi) + ncyl(diam(2*i),csi,h,W);
                            ST(gen(2*i)+1,hi) = ST(gen(2*i)+1,hi) + ncyl(diam(2*i),csi,h,W)*Acyl1(diam(2*i),W);
                        end
                        SA(gen(2*i)+1) = SA(gen(2*i)+1) + diam(2*i)*pi*csi;
                        VOL = VOL + pi*csi*diam(2*i)^2/4;
                    else
                        %add the segment measures to the total counts;
                        for hi=1:Hcount+1
                            h = h0 + (h1-h0)*(hi-1)/Hcount;
                            NT(gen(2*i)+1,hi) = NT(gen(2*i)+1,hi) + ncyl(diam(2*i),csi-diam(2*i)/2,h,W);
                            Ntip(hi) = Ntip(hi) + nhsphere(diam(2*i)/2,h,W);
                            ST(gen(2*i)+1,hi) = ST(gen(2*i)+1,hi) + ncyl(diam(2*i),csi-diam(2*i)/2,h,W)*Acyl1(diam(2*i),W);
                            Stip(hi) = Stip(hi) + nhsphere(diam(2*i)/2,h,W)*Asphere(diam(2*i)/2,W);
                        end
                        SA(gen(2*i)+1) = SA(gen(2*i)+1) + diam(2*i)*pi*(csi-diam(2*i)/2);
                        SAtip = SAtip + 2*pi*diam(2*i)^2/4;
                        VOL = VOL + pi*(csi-diam(2*i)/2)*diam(2*i)^2/4 + 2*pi*diam(2*i)^3/8/3;
                    end;
                    
                    LGTH = LGTH + csi;

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
             if diam(i)>0
             tipdiam(ct) = diam(i);
             ct = ct+1;
             end
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

        