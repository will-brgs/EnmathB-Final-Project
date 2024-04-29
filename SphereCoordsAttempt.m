%% Enmath B Final Project:
%% Introduction

% * Author:                   Will Burgess
% * Class:                    ESE 319
% * Date:                     Created 4/28/2024, Last Edited 
%% Housekeeping
close all
clc
clear
filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. EnMath B\Final Project\Figures";

%% Objective 1: Find the Solution to the BVP
% All relevant caclulations are done by hand and carried out in the
% objective 2 calculation
%% Objective 2: Evaluate the First Six Nonzero Terms of u(r,theta)
% Results are shown in hand calculation - Matlab is used to assist in the
% integration and constant coefficient.

obj2_outputs = zeros(6,1);
% Index 0
for n = 0:5
index = n + 1;
constant = ((2.*n) + 1)/(2.^(n+1));
integralIn = @(x) (1/6) .* besselk(0,sqrt((1+x)/2)) .* legendreP(n, x);
integralOut = integral(integralIn,-1,1);
obj2_outputs(index) = integralOut * constant;
end
% Direct inspection of the obj2_outputs variable is applied to hand
% calculations

%% Objective 3: Genereate a Figure Displaying Concentration inside Sphere

%independent Vars
numvals = 100;
theta = linspace(0,2*pi,numvals); % Covers entire sphere, use radians
r = linspace(0,3,numvals); % diamater is 6cm
phi = linspace(-1*pi, pi, numvals);
outputMerged = zeros(numvals,numvals);

% fig1 = figure();
% sgtitle('Contribution to Concentration for each Nonzero Index of N')
for n = 0:5
% Output Calculation for u(r,theta)
figindex = n + 1;
constant = ((2.*n) + 1)/(2.^(n+1));
integralIn = @(x) (1/6) .* besselk(0,sqrt((1+x)/2)) .* legendreP(n, x);
integralOut = integral(integralIn,-1,1);
outputTerm = (constant .* integralOut .* r.^n);

%perform matrix multiplication to obtain unique range for each value of r & theta accross numvals
outputTerm = (outputTerm.' * legendreP(n,cos(theta))); 
outputMerged = outputMerged + outputTerm;

% [meshR, meshTheta, meshPhi] = meshgrid(r,theta, phi);
[x,y,z] = sph2cart(theta, phi, r);
[meshX, meshY, meshZ] = meshgrid(x,y,z);
% x = reshape(x, [], numvals);
% y = reshape(y, [], numvals);
% z = reshape(z, [], numvals);

% Plotting Mesh Contributions
% subplot(3,2,figindex)
% Z = squeeze(outputTerm);
% surf(meshR .* cos(meshTheta), meshR .* sin(meshTheta), Z, 'edgecolor','none');
% xlabel('X Position (cm)'),ylabel('Y Position (cm)'),title(['n = ' num2str(n)])
% colormap('jet');
% bar = colorbar;
% caxis([-0.1, 0.14])
% ylabel(bar, 'Quantity Concentration');
% view(0,90);
% grid off
end
%% Generate Merged Mesh Figure
fig2 = figure();
%Z = squeeze(outputMerged);
surf(meshX, meshY, meshZ, outputMerged, 'edgecolor','none');
xlabel('X Position (cm)'),ylabel('Y Position (cm)'),title({'Quantity Concentration of a Uniform Slice of the Sphere' 'Merged Results of First 6 Nonzero Coefficients'})
colormap('jet');
bar = colorbar;
ylabel(bar, 'Quantity Concentration');
grid off
%view(0,90);


%% Export Figures
% exportgraphics(fig1, fullfile(filepath, 'Independent Output.jpg'), 'resolution', 300);
% exportgraphics(fig2, fullfile(filepath, 'Merged Output.jpg'), 'resolution', 300);
