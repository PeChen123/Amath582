% Clean workspace
 clear all; close all; clc
 % Preamble from CP1_sample.m.  The autograder will take care of any files
 % needed for this step.  Please only submit solution file.
 load('Kraken.mat'); %load data matrix
 L = 10; % spatial domain
 n = 64; % Fourier modes
 x2 = linspace(-L,L,n+1); 
 x = x2(1:n); 
 y =x; 
 z = x; %create 3D axis arrays with 64 points
 k = (2*pi/(2*L))*[0:(n/2 - 1), -n/2:-1]; %create frequency array and rescale them to be 2pi periodic 
 ks = fftshift(k); %shift values to order them correctly
 %Create 3D grids for both spatial domain and frequency domain 
 [X,Y,Z] = meshgrid(x,y,z);
 [Kx,Ky,Kz] = meshgrid(ks,ks,ks);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % sum all realizations in frequency space
 % after your for loop ends, save the sum as variable A1
 Unt_sum = zeros(n,n,n);
 for j = 1:49
     Un(:, :, :) = reshape(Kraken(:, j), n, n, n);
     Unt_sum = Unt_sum + fftn(Un);     
 end

 A1 = Unt_sum;
 % Average the sum over the 49 realizations (i.e., A1/49) and save as A2
 A2 = A1/49;
 % find the peak frequencies in x, y, and z directions; i.e., find the
 % max in each direction of the normalized sum A2.
 % save these variables as A3, A4, and A5
Max = max(A2,[],'all');
LI = find(A2==Max);
ind = [LI];
[kx0,ky0,kz0] = ind2sub([n,n,n],ind);
 A3 = k(kx0);
 A4 = k(ky0);
 A5 = k(kz0);
%create an appropriate Gaussian filter and save it as A6
ki = Kx(kx0,ky0,kz0);
kj = Ky(kx0,ky0,kz0);
kk = Kz(kx0,ky0,kz0);
tau = 1/L;
gFilter  = exp(-tau*((Kx - ki).^2 + (Ky - kj).^2 +(Kz - kk).^2));
A6 = gFilter;
% Using the peak frequencies for the filtered signal, estimate the x, y, and z coordinates of the Kraken over time and save as A7, A8, A9
Kraken1 = zeros(3,49);
for j = 1:49
    Un = reshape(Kraken(:, j), n, n, n);
    Untf = fftn(Un).*gFilter;
    unf = ifftn(Untf);
    Max1 = max(unf,[],'all');
    LI1 = find(unf==Max1);
    ind = [LI1];
    [Kraken1_x,Kraken1_y,Kraken1_z] = ind2sub([n,n,n],ind);
    Kraken1(1,j) = X(Kraken1_x,Kraken1_y,Kraken1_z);
    Kraken1(2,j) = Y(Kraken1_x,Kraken1_y,Kraken1_z);
    Kraken1(3,j) = Z(Kraken1_x,Kraken1_y,Kraken1_z);   
end
A7 = Kraken1(2,:); % x coordinates 
A8 = Kraken1(1,:); % y coordinates
A9 = Kraken1(3,:); % z coordinates

% Plot the location in x-y-z space over time for your report (not for the autograder)
figure(1)
plot3(Kraken1(2,:), Kraken1(1,:),Kraken1(3,:) ,'b-o'), grid on;
title('location of Kraken '), 
xlabel('x'), ylabel('y'), zlabel('z')

figure(2)
threshold = 0.6;
isosurface(X,Y,Z,abs(unf)/max(abs(unf(:))),threshold)
xlabel('x'), ylabel('y'), zlabel('z')
axis([min(x) max(x) min(y) max(y) min(z) max(z)]),
grid on 

% Plot the projection onto the x-y plane for your reprot (not for the autograder)
figure(3)
plot(Kraken1(2,:), Kraken1(1,:),'b-o');
title('projection onto the x-y plane'), 
xlabel('x'), ylabel('y'),

% Save a table of the x-y-z coordinates for your report (not for the autograder)
T = array2table(Kraken1);
sympref('FloatingPointOutput',1)
latex_table = latex(sym(Kraken1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include all helper functions below since the autograder can only accept
