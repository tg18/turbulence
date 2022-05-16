%% AYUSH BISEN 21105025
%% TURBULENCE QUESTION3
%%
clc;
clear all;
% File which includes three columns, y U U''
data= load ('velocity1.dat');  % load your velocity profile 1
%data= load ('velocity2.dat');  % load your velocity profile 2, Note: Kindly uncomment the data file for the velocity profile for which you want to plot the graphs.
%% ----------------------------------------

% Number of Chebyshev modes
n= 99; N1 = n+1; 
j = (0:n)'; X = cos(pi*j/n);
%-------------------------------------

% Scaling for mapping from -1 to 1
L = data(end,1);
scale=L/2; % Physical domain => [0, L]
Y=X*scale + scale; % Mapping from Physical to Chebyshev domain
%----------------------------------

% Velocity and its second derivative are being evaluated at scaled Y
U = interp1(data(:,1), data(:,2), Y);
U2 = interp1(data(:,1), data(:,3), Y);
fac = 1/scale;
%-------------------------------------

%% For singel value of alpha and Re, for example, Note : Kindly comment this section when plotting stability curve of figure(3)
% alpha = 0.25;
% length(alpha)
% Rey = 600;
%% -------------------------------

%% For a range of alpha and Reynolds number, for example, Note : Kindly uncomment this section when plotting stability curve of figure(3)
 alpha = 0.1:0.1:1.5;
 length(alpha)
 Rey = 10:5:3000;
%% -----------------------------------------

% Main program; try to understand before you run it
ci = zeros(length(alpha),length(Rey));
cr = zeros(length(alpha),length(Rey));
for ii = 1:length(alpha)
    ii
    for jj = 1:length(Rey)
        al = alpha(ii);
        R = Rey(jj);
        zi = sqrt(-1); a2 = al^2; a4 = a2^2; er = -200*zi;
        [D0,D1,D2,D3,D4] = Dmat(n); % Read about it in the book.
        D1 = fac*D1;
        D2 = (fac^2)*D2;
        D4 = (fac^4)*D4;
        B = (D2-a2*D0);
        A = (U*ones(1,N1)).*B-(U2*ones(1,N1)).*D0-(D4-2*a2*D2+a4*D0)/(zi*al*R);
        A = [er*D0(1,:); er*D1(1,:); A(3:n-1,:); er*D1(N1,:); er*D0(N1,:)];
        B = [D0(1,:); D1(1,:); B(3:n-1,:); D1(N1,:); D0(N1,:)];
        d = (inv(B)*A);
        [vv, c] = eig(d);                 % eigenvalues are being evaluated using eig function
         [mxci,I] = max(max(imag(c)));  % useful for contour plots
         ci(ii,jj)= mxci;
         cr(ii,jj) =real(c(I,I)); 

    end;
end

%x=sprintf('the most unstable eigenvalue for alpha = 0.25 and Re = 600 is %d',ci);
%disp(x);
%% plotting the absolute value of the eigenfunctions, |phi(y)|, corresponding to the wall-normal disturbance velocity

% figure(1)
% phi_y = vv(:,I);
% plot(abs(phi_y),Y)
% ylabel('Y-axis');
% xlabel('v'); 
% grid on
% 
% ylim([18 20])
% title('Y v/s v','color','k');
% 
% %% plotting the absolute value of the eigenfunctions, |phi(y)|, corresponding to the stream-wise disturbance velocity
% figure(2)
% phi_for_u = ((zi/al).*(D1*U))./U;
% plot(abs(phi_for_u ),Y)
% ylabel('Y-axis');
% xlabel('u'); 
% grid on
% ylim([18 20])
% title('Y v/s u ','Color','k');

%% Plotting the neutral stability curves in the plane of alpha and Re

figure(3)
[X Y] = meshgrid(Rey,alpha);
contour(X,Y,ci,[0.001 0.001],'LabelSpacing',350,'LineWidth',3)
colorbar 
xlabel('Re')
ylabel('alpha')
title('neutral stability curve','Color','k')
ylim([0 0.5])
grid on