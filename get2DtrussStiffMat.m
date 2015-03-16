function Kmem = get2DtrussStiffMat(coord, E, area)

% Compute stiffness matrix of 2D truss member

x1 = coord(1,1);
y1 = coord(1,2);
x2 = coord(2,1);
y2 = coord(2,2);

L = sqrt( (x2-x1)^2 + (y2-y1)^2 );

theta = atan2(y2-y1,x2-x1);

B = [-cos(theta) -sin(theta) cos(theta) sin(theta)];

Kmem = B'*(E*area/L)*B;



function [Fmem] = get2DtrussMemberForce(coord, E, area, u )

% Compute stiffness matrix of 2D truss member
alpha = 1.2* 10^-5;

x1 = coord(1,1);
y1 = coord(1,2);
x2 = coord(2,1);
y2 = coord(2,2);

L = sqrt( (x2-x1)^2 + (y2-y1)^2 );

theta = atan2(y2-y1,x2-x1);

B = [-cos(theta) -sin(theta) cos(theta) sin(theta)];

Fmem = (E*area/L)*B*u; 



