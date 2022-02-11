% source: https://au.mathworks.com/help/symbolic/find-tangent-plane-and-normal-line.html

clear; close all; clc
syms r [1 3] 
f = r*r.'
feqn = (f == 14); %radius is 14^(0.5)

fimplicit3(feqn)
axis equal
axis([-6 6 -6 6 -6 6])

fgrad = gradient(f,r)
size(fgrad)

% Define the equation for the tangent plane. Use the subs function to
% evaluate the gradient at the point r0
r0 = [-2,1,3];
r_v = [-1.5,0.8,2.5];

fplane = (r-r0)*subs(fgrad,r,r0);
%  plot
hold on
% scatter3(r0(1),r0(2),r0(3),'ro')
plot3([r0(1);r_v(1)],[r0(2);r_v(2)],[r0(3);r_v(3)],'r', 'LineWidth',2)


%equation for normal line
syms t
n = r0 - t*subs(fgrad,r,r0).'; % normal line

fplot3(n(1),n(2),n(3),[0 0.1],'b-','LineWidth',4)

syms n(t)
n(t) = r0 - t*subs(fgrad,r,r0).'; % normal line


V = [r0;r_v];
new = [n(0); n(0.1)];
%rotate about normal line
final_temp = V-dot(V,new).*new;
final = [n(0);final_temp(1,:)];
% final=final_temp;

A = r0;
B = r_v; %starting line
C = final(2,:); %normal line

S1 = B - A;
S2 = C - A;
Theta = atan2(norm(cross(S1, S2)), dot(S1, S2))

test = S2*sind(Theta);
%find distance between point and normal line
% plot3([r0(1); final(2,1)],[r0(2); final(2,2)],[r0(3); final(2,3)],'g', 'LineWidth',3)
plot3([r0(1); r0(1)+test(1)],[r0(2); r0(2)+test(2)],[r0(3); r0(3)+test(3)],'g', 'LineWidth',3)


%% method of finding angle between rotated line and original points
A = r0;
B = r_v; %starting line
C = final(2,:); %normal line

S1 = B - A;
S2 = C - A;
Theta = atan2(norm(cross(S1, S2)), dot(S1, S2))

