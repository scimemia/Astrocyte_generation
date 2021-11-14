global r0 R c a d D N M f;
global SA VOL LGTH brch; 
global SAtip VOLtip;
global Diamhr Vh Vhtip SAh SAhtip cth freq_leaves;
global G_prim L_prim lh
global xi lambda count col ct tipdiam; 
global LFT RGT stp PT;
global CV1 CV2 BK Volcore;

Volcore = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Cell and branch geometry
%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of primary branches P;
P = 1;

%Radius of the soma r0;
r0 = 5;

%Maximum basic length of a segment along the primary branch;
c = 25*rand;

%Variability added to the basic length for each segment along promary branch;
R = 10*rand;

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


%create diameter vectors;
%D are the diameters of the M+1 segments along the primary branch;
%d are the diameters of the M secondary branches;
D = zeros(M+1,1);
d = zeros(M,1);

%initialize the diamter of the primary branch at the exit form the soma;
D(1) = 5;

%fraction bias in the 3/2 rule (random for each cell between 1 and 1.5);
f = 1+0.5*rand;

%calculate diameters using the a=3/2 rule, and the bias fraction f;
for k = 1:M  
    D(k+1) = f*D(k)/((1+f^a)^(1/a));
    d(k) = D(k)/((1+f^a)^(1/a));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Leaf distribution and measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance between leaves
freq_leaves = 5+3*rand;
%leaf length
lh = 0.7;

%initialize Surface Area (SA) and Volume (V) of leaves and their tips;
SAh = 0;
SAhtip = 0;
Vh = 0;
Vhtip = 0;

%initialize leaf counter and vector for leaf diameters;
cth = 1;
Diamhr = zeros(50000,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing cell measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vector for number of branches per generation;
gn = zeros(N+1,1); 

%vector for total length per generation;
lg = zeros(N+1,1);

%reserve a large enough vector to record branch tip diameters (only the nonzero entries
%will be included in the statistical evauation)
tipdiam = zeros(7*2^N,1);
%initialize count for tipdiam (will increase each time a branch tip is reached;
ct = 1;

%total surface area of the branching tree and tips;
SFA = zeros(M+N,1);
SFAtip = 0;

%total cell volume of the branching tree and tips;
VOLM = zeros(M+N,1);
VOLMtip = 0;

%total branch length;
LENGTH = 0;
 
%number of branch points;
BR = 0; 

XI = zeros(M,P);
LAMBDA = zeros(M,P);
COUNT = 0;

%cover parameters; 
LFT = -200; RGT = 200;
stp = 2;
PT = floor((RGT-LFT)/stp)+1;
CV1 = zeros(PT,PT,PT);
CV2 = zeros(PT,PT,PT);
BK = 1;


%fix frame size and camera angle;
figure;
az = pi*(rand-0.5);
el = pi*(rand-0.5);
%axis([-100 100 -100 100 -100 100]);
hold on;
axis square;
view(az,el);
axis off;

%plot the soma;
[x,y,z] = sphere(50); surf(r0*x,r0*y,r0*z,'FaceColor',[0 0 0.59],'EdgeColor',[0 0 0.59]);
hold on;

%start with P branches in the first generation;
gn(1) = P;

for primbr = 1:P;
        %directions of the primary branches (chosen randomly);
        U = 0.5*[1 1 1];
        V = rand(3,1)-U';
        U1 = V/norm(V);

        %call the subroutine that builds each primary branch;
        %all branch measures are computed, values are brought in as global variables, 
        %and used to update the cell-wide measures;

        primary_branch_cube_cover(U1);
        %number of branches per level;
        gn = gn + G_prim;
        %total length per level;
        lg = lg + L_prim;
        
        %surface area;
        SFA(:) = SFA(:) + SA(:);
        SFAtip = SFAtip + SAtip;
        %volume;
        VOLM(:) = VOLM(:) + VOL(:);
        VOLMtip = VOLMtip + VOLtip;
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

end

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

% figure;
% az = pi*(rand-0.5);
% el = pi*(rand-0.5);
% axis([-100 100 -100 100 -100 100]);
% hold on;
% axis square;
% view(az,el);
% axis off;
    x = LFT:stp:RGT;
    y = LFT:stp:RGT;
    z = LFT:stp:RGT;
hold on;
     for i = 1:PT
         for j = 1:PT
             for k = 1:PT
                if (CV2(i,j,k)>0)
                    plotcube([stp stp stp],[x(i) y(j) z(k)],0.1,'r');
                end
             end
         end
     end
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRINT COVER GENERATED MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%Bindocci paper estimate of core volume fraction
(4*pi*5^3/3+Volcore)/(4*pi*5^3/3+sum(VOLM+VOLMtip))

%Medvedev estimate of GV fraction
(sum(VOLM)-Volcore)/length(find(CV2(:)))/stp^3

