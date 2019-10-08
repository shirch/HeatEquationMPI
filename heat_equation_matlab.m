% Numerical scheme used is an implicit second order central
%difference in space(5-point difference)
%Specifying parameters
nx=502; %Number of steps in space(x)
ny=302; %Number of steps in space(y)
dx=500/(nx-1);; %Width of space step(x)
dy=300/(ny-1); %Width of space step(y)
x=0:dx:500; %Range of x(0,2) and
%specifying the grid points
y=0:dy:300; %Range of y(0,2) and
%Plotting the solution
u=importdata('C:\Users\shirc\Downloads\outputMatrix.txt');
surf(x,y,u');
shading interp
colorbar
title('2-D Heat Equation')
xlabel(' x \rightarrow')
ylabel('{\leftarrow} y')
zlabel('Solution profile (P) \rightarrow')