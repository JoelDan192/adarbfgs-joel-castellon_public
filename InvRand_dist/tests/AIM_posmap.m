% Is R -ZRZ semi-positive definie?

n = 10;
R = rand(n,n);
R = R'*R;

% Z = S(S^TAS)S^T
p =4;
A = rand(n,n);
A = A'*A; 
S = rand(n,p);  AS = A*S;
Z = S*inv(S'*A*S)*S';
AZ = A*Z;
Rn = R - (AZ')*R*AZ;
% eig(R)
% eig(Rn)
X = rand(n,n);
X = X'*X;
Xn = X + (AZ')*(inv(A) -X)*AZ;
eig(X)
eig(Xn)
eig(inv(A) -X)
%% Norm test
norm(R)
norm(Rn)

norm(R,'fro')
norm(Rn,'fro')

%% For A pos def and W = inv(A)
n = 10;
R = rand(n,n);
R = R'*R;

% Z = S(S^TAS)S^T
p =4;
A = rand(n,n);
A = A'*A; 
S = rand(n,p);  AS = A*S;
Z = S*inv(S'*A*S)*S';
AZ = A*Z;
X = rand(n,n);
X = X'*X;
Xn = X +Z - (AZ')*X*AZ;
eig(X)
eig(Xn)
%eig(inv(A) -X)


%% Curvature tests

AS'*Rn*AS
eig(S'*A*S)
%
V = null(S');
V'*Rn*V
V'*R*V
eig(V'*R*V)

%% Visualize

f = @(x,y)(Rn(1,1)*x^2 +Rn(2,2)*y^2+2*Rn(1,2)*x*y ) 
ezsurf(f)