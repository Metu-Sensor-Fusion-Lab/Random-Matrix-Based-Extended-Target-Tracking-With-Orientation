function [extent] = plot_ellipse(par,line_style, color, line_width)
% Plotting Ellipse
% Input:
%        par:    the parameters of the ellipse [x_pos y_pos x_len y_len tet]
% Output:
%        extent: the handle of the plot
%

L = [cos(par(5)) -sin(par(5)); sin(par(5)) cos(par(5))]; 
alpha = 0:pi/100:2*pi;
x = par(3)*cos(alpha);
y = par(4)*sin(alpha);

rotated = L * [x; y];
xcoords = rotated(1,:)+ par(1);
ycoords = rotated(2,:)+ par(2);

extent = plot(xcoords,ycoords,'LineStyle',line_style,'color',color,'LineWidth',line_width);

end