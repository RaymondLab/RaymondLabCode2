function [A, lamda] = EllipseFitByTaubin(XY)
%
%   Ellipse fit by Taubin's Method published in
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With 
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: A = [a b c d e f]' is the vector of algebraic 
%             parameters of the fitting ellipse:
%             ax^2 + bxy + cy^2 +dx + ey + f = 0
%             the vector A is normed, so that ||A||=1
%
%     Among fast non-iterative ellipse fitting methods, 
%     this is perhaps the most accurate and robust
%
%     Note: this method fits a quadratic curve (conic) to a set of points;
%     if points are better approximated by a hyperbola, this fit will 
%     return a hyperbola. To fit ellipses only, use "Direct Ellipse Fit".
%
% Conversion to standard form:
% a = A(1); b = A(2)/2; c = A(3); d = A(4)/2; f = A(5)/2; g = A(6);
% 
% center(1) = (c*d - b*f)/(b^2-a*c);
% center(2) = (a*f - b*d)/(b^2-a*c);
% 
% sem(1) = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
% sem(2) = sqrt( 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g) / ((b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c))));
% 
% if b == 0 && a < c
% phi = 0;
% elseif b == 0 && a > c
% phi = 0.5*pi;
% elseif b ~= 0 && a < c
% phi = 0.5* acot((a-c)/(2*b));
% else
% phi = 0.5*pi + 0.5* acot((a-c)/(2*b));
% end

centroid = mean(XY);   % the centroid of the data set

Z = [(XY(:,1)-centroid(1)).^2, (XY(:,1)-centroid(1)).*(XY(:,2)-centroid(2)),...
     (XY(:,2)-centroid(2)).^2, XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)];
M = Z'*Z/size(XY,1);

P = [M(1,1)-M(1,6)^2, M(1,2)-M(1,6)*M(2,6), M(1,3)-M(1,6)*M(3,6), M(1,4), M(1,5);
     M(1,2)-M(1,6)*M(2,6), M(2,2)-M(2,6)^2, M(2,3)-M(2,6)*M(3,6), M(2,4), M(2,5);
     M(1,3)-M(1,6)*M(3,6), M(2,3)-M(2,6)*M(3,6), M(3,3)-M(3,6)^2, M(3,4), M(3,5);
     M(1,4), M(2,4), M(3,4), M(4,4), M(4,5);
     M(1,5), M(2,5), M(3,5), M(4,5), M(5,5)];

Q = [4*M(1,6), 2*M(2,6), 0, 0, 0;
     2*M(2,6), M(1,6)+M(3,6), 2*M(2,6), 0, 0;
     0, 2*M(2,6), 4*M(3,6), 0, 0;
     0, 0, 0, 1, 0;
     0, 0, 0, 0, 1];

[V,D] = eig(P,Q);

[Dsort,ID] = sort(diag(D));
lamda = Dsort(1);

A = V(:,ID(1));
A = [A; -A(1:3)'*M(1:3,6)];
A4 = A(4)-2*A(1)*centroid(1)-A(2)*centroid(2);
A5 = A(5)-2*A(3)*centroid(2)-A(2)*centroid(1);
A6 = A(6)+A(1)*centroid(1)^2+A(3)*centroid(2)^2+...
     A(2)*centroid(1)*centroid(2)-A(4)*centroid(1)-A(5)*centroid(2);
A(4) = A4;  A(5) = A5;  A(6) = A6;
A = A/norm(A);

end  %  Taubin

