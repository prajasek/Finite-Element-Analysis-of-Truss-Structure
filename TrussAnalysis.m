%% A simple code for linear structural analysis
%

%% Begin Input (units used: forces=kN; lengths=mm)

syms dt1 dt2 dt3 dt4 dt5 dt6
syms Ft
dt2 = 0 ;
dt4 = 0 ;
dt5 = 0 ; 
dt6 = 0 ;

dt1 =  45;
dt3 =  -20;

dT = [dt1 dt2 dt3 dt4 dt5 dt6];



alpha = 1.2*10^-5 ;
% joint locations
coord = [0 0; 6 0; 0 6; 6 6]*1000;

% support conditions
bcs = [1 1; 0 1; 0 0; 0 0];

% member connnectivity
connect = [1 2; 1 3; 3 4; 2 4; 3 2; 1 4];

% properties
area =  [5000 5000 5000 5000 3000 3000];   
E = 200;  







% loads
c = E*alpha
Ptemp = c*[ 0 ,0;((dt5*area(5)/sqrt(2)) + dt1*area(1)) ,0;
(- dt3*area(3)) - ((dt5*area(5) /sqrt(2))),...
(dt2*area(2))+ ((dt5*area(5)/sqrt(2)));...
 (dt3*area(3))+ (dt6*area(6)/sqrt(2)) ,...
 dt4*area(4) + (dt6*area(6)/(sqrt(2)))];      %thermal loads

P = [0 0; 0 0; 10 0; 0 0];                    %mechanical loads

Ft    =  (E*alpha)* [dT(1)*area(1); dT(2)*area(2); dT(3)*area(3);...
                     dT(4)*area(4); dT(5)*area(5); dT(6)*area(6)];


% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculations
% Code from this point onward does not depend on number of DOF per joint,
% number of joints per member, type of member etc. specific to input data,
% except a call to a function to compute the element stiffness matrix

numJts = size(coord,1); % Total number of joints in the model
% The following is an example of an integrity check on the input data,
% however, no more such checks are done in the following. It is assumed
% that the input data is consistent.
if (numJts ~= size(bcs,1))
    fprintf('Error: number of coord must equal number of boundary cond\n')
    return
end

numMem = size(connect,1); % Total number of members in the model

numDofPerJt = size(bcs,2); % Number of DOF per node.

numJtsPerMem = size(connect,2); % Number of nodes per element

numDofPerMem = numJtsPerMem * numDofPerJt; % Number of DOF associated with
                                           % an element

% count and number the free DOF
numFreeDof = 0;
jtDofNum = zeros(size(bcs)); % an array to hold the joint DOF numbers
for m = 1:numJts
    for n = 1:numDofPerJt
        if (bcs(m,n) == 0) % i.e. a free DOF
            numFreeDof = numFreeDof + 1;
            jtDofNum(m,n) = numFreeDof;
        end
    end
end

memDofNum = zeros(numMem, numDofPerMem); % an array to hold DOF numbers 
                                         % associated with elements
for m = 1:numMem
    for n = 1:numJtsPerMem
        memDofNum(m,(n-1)*numDofPerJt+1:n*numDofPerJt) = ...
                                               jtDofNum(connect(m,n),:);
    end
end

% build stiffness matrix
K = zeros(numFreeDof); 
for m = 1:numMem
    % Compute B'*E*B for the member
    Kmem = get2DtrussStiffMat(coord(connect(m,:),:), E, area(m));
    
    % plug it into the right place in the stiffness matrix K
    maskActiveDof = memDofNum(m,:) > 0;
    indexActiveDof = memDofNum(m,maskActiveDof);
    
    K(indexActiveDof,indexActiveDof) = K(indexActiveDof,indexActiveDof)+...
        Kmem(maskActiveDof,maskActiveDof);
end

% load vector at the free DOF (is there a better way ??)
Pfree = zeros(numFreeDof,1);
for m = 1:numJts
    for n = 1:numDofPerJt
        if (jtDofNum(m,n) > 0)
            Pfree(jtDofNum(m,n)) = P(m,n);
        end
    end
end

Ptempfree = zeros(numFreeDof,1);
for m = 1:numJts
    for n = 1:numDofPerJt
        if (jtDofNum(m,n) > 0)
            Ptempfree(jtDofNum(m,n)) = Ptemp(m,n);
        end
    end
end
Ptotal = Pfree + Ptempfree;

% Solve
ufree = K\Ptotal;

% all displacements
u = zeros(size(bcs));
for m = 1:numJts
    for n = 1:numDofPerJt
        if (jtDofNum(m,n) > 0)
            u(m,n) = ufree(jtDofNum(m,n));
        end
    end
end

Fmemberthermal    =  (E*alpha)* [dt1*area(1); dt2*area(2); dt3*area(3);...
                                 dt4*area(4); dt5*area(5); dt6*area(6)];

% compute member forces
F = zeros(numMem,1);
for m = 1:numMem
    maskActiveDof = memDofNum(m,:) > 0;
    indexActiveDof = memDofNum(m,maskActiveDof);
    
    umem = zeros(numDofPerMem,1);
    umem(maskActiveDof) = ufree(indexActiveDof);
    
    [F(m)] = get2DtrussMemberForce(coord(connect(m,:),:),E, area(m), umem);
end
%% Results
%Displacements
u
