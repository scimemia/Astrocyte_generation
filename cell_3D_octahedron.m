global W r0 R c a d D N M Csholl Nsholl step Hcount h0 h1 NT SA ST VOL LGTH brch Ntip Stip SAtip;
global G_prim L_prim xi lambda count d f col ct tipdiam;


%%%%%%%%%%%%%%%%%%%%%%%%%
%Cell and branch geometry
%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of primary branches P;
P = 3;

%Radius of the soma r0;
r0 = 5;

%Maximum basic length of a segment along the primary branch;
c = 15*rand;

%Variability added to the basic length for each segment along promary branch;
R = 5*rand;

%exponent in the diameter rule (in case we want to use anything besides the 3/2 rule);
a = 3/2;

%max allowable number of branch points along a primary branch;
M = 15;

%max allowable number of subsequent branch levels on secondary branches;
N = 12;

%color coding for branch level, jet color map (blue to red)
col = [0.00	0.00	0.59
0.00	0.00	0.94
0.00	0.36	1.00
0.00	0.75	1.00
0.13	1.00	0.91
0.56	1.00	0.44
0.88	1.00	0.12
1.00	0.69	0.00
1.00	0.25	0.00
0.87	0.00	0.00
0.50	0.00	0.00
0.50	0.00	0.00];

%fraction bias in the 3/2 rule (random for each cell between 1 and 1.5);
f = 1+0.5*rand;

%create diameter vectors;
%D are the diameters of the M+1 segments along the primary branch;
%d are the diameters of the M secondary branches;
D = zeros(M+1,1);
d = zeros(M,1);

%initialize the diamter of the primary branch at the exit form the soma;
D(1) = 5;

%calculate diameters using the a=3/2 rule, and the bias fraction f;
for k = 1:M  
    D(k+1) = f*D(k)/((1+f^a)^(1/a));
    d(k) = D(k)/((1+f^a)^(1/a));
end

%%%%%%%%%%%%%%%%%%%%%
%Transporter geometry
%%%%%%%%%%%%%%%%%%%%%
%cytoplasmic insertion height (Hcount samples from h0=3nm to h1=6.5nm);
Hcount = 10;
h0 = 0.003;
h1 = 0.0065;
%transporter width 8nm);
W = 0.008;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing cell measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize sholl profile (with 4000 shells, step 0.05);
Nsholl = 4000;
step = 0.05;
Csholl = zeros(Nsholl+1,1);

%vector for number of branches per generation;
gn = zeros(N+1,1); 

%vector for total length per generation;
lg = zeros(N+1,1);

%reserve a large enough vector to record tip diameters (only the nonzero entries
%will be included in the statistical evauation)
tipdiam = zeros(7*2^N,1);
%initialize index for tipdiam (will increase each time a tip value is aadded;
ct = 1;

%total surface area of the branching tree;
SFA = zeros(M+N,1);
SFAtip = 0;

%total cell volume of the branching tree;
VOLM = 0;

%total branch length;
LENGTH = 0;
 
%number of branch points;
BR = 0; 

%total number of cell transporters in the branching tree (depends on cytoplasmic insertion);
NTR = zeros(M+N,Hcount+1);
NTRtip = zeros(1,Hcount+1);

%total area occupied by transporters in the branching tree (depends on cytoplasmic insertion);
SFT = zeros(M+N,Hcount+1);
SFTtip = zeros(1,Hcount+1);

XI = zeros(M,P);
LAMBDA = zeros(M,P);
COUNT = 0;

%directions of the primary branches (chosen randomly);
U = 0.5*[1 1 1];
V = rand(3,1)-U';
U1 = V/norm(V);
V = rand(3,1)-U';
U2 = V/norm(V);
V = rand(3,1)-U';
U3 = V/norm(V);
V = rand(3,1)-U';
U4 = V/norm(V);
V = rand(3,1)-U';
U5 = V/norm(V);
V = rand(3,1)-U';
U6 = V/norm(V);
V = rand(3,1)-U';
U7 = V/norm(V);

%fix frame size and camera angle;
figure;
az = pi*(rand-0.5);
el = pi*(rand-0.5);
axis([-100 100 -100 100 -100 100]);
hold on;
axis square;
view(az,el);
axis off;

%plot the soma;
[x,y,z] = sphere(50); surf(r0*x,r0*y,r0*z,'FaceColor',[0 0 0.59],'EdgeColor',[0 0 0.59]);
hold on;

%start with P branches in the first generation;
gn(1) = P;

%call the subroutine that builds each primary branch;
%all branch measures are computed, values are brought in as global variables, 
%and used to update the cell-wide measures;

primary_branch(U1);
for hi = 1:Hcount+1
    %transporter number (depends on cytoplasmic insertion depth);
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    %transporter area (depends on cytoplasmic insertion depth);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
%number of branches per level;
gn = gn+G_prim;
%total length per level;
lg = lg+L_prim;
%surface area;
SFA(:) = SFA(:) + SA(:);
SFAtip = SFAtip + SAtip;
%volume;
VOLM = VOLM +VOL;
%total branch length;
LENGTH = LENGTH + LGTH;
%number of branch points;
BR = BR + brch;
%lengths of secondary segments;
XI(:,1) = xi;
%lengths of primary segments;
LAMBDA(:,1) = lambda;
%total number of secondary branches (for the chi stats);
COUNT = COUNT + count;

primary_branch(U2);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,2) = xi;
LAMBDA(:,2) = lambda;
COUNT = COUNT + count;

primary_branch(U3);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,3) = xi;
LAMBDA(:,3) = lambda;
COUNT = COUNT + count;

primary_branch(U4);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,4) = xi;
LAMBDA(:,4) = lambda;
COUNT = COUNT + count;

primary_branch(U5);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,5) = xi;
LAMBDA(:,5) = lambda;
COUNT = COUNT + count;

primary_branch(U6);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,6) = xi;
LAMBDA(:,6) = lambda;
COUNT = COUNT + count;

primary_branch(U7);
for hi = 1:Hcount+1
    NTR(:,hi) = NTR(:,hi) + NT(:,hi);
    NTRtip(hi) = NTRtip(hi) + Ntip(hi);
    SFT(:,hi) = SFT(:,hi) + ST(:,hi);
    SFTtip(hi) = SFTtip(hi) + Stip(hi);
end
gn = gn+G_prim;
lg = lg+L_prim;
SFA = SFA + SA;
SFAtip = SFAtip + SAtip;
VOLM = VOLM +VOL;
LENGTH = LENGTH + LGTH;
BR = BR + brch;
XI(:,7) = xi;
LAMBDA(:,7) = lambda;
COUNT = COUNT + count;

%plot frame with sizing vectors;
px=[70 70 70];
vx=[20 0 0];
drawVector(px, vx, 'cyan', 1);
axis equal;

pz=[70 70 70];
vz=[0 0 20];
drawVector(pz, vz, 'black', 1);
axis equal;

py=[70 70 70];
vy=[0 20 0];
drawVector(py, vy, 'red', 1);
axis equal;



