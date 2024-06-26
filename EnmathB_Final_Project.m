%% Enmath B Final Project:
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 319
% * Date:                     Created 4/28/2024, Last Edited 5/6/2024
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
constant = ((2.* n) + 1)/(2.^(n+1));
integralIn = @(x) (1/6) .* besselk(0,sqrt((1+x)/2)) .* legendreP(n, x);
integralOut = integral(integralIn,-1,1);
obj2_outputs(index) = integralOut .* constant;
end
% Direct inspection of the obj2_outputs variable is applied to hand
% calculations
%% Objective 3: Genereate a Figure Displaying Concentration inside Sphere
%independent Vars
numvals = 1000;
theta = linspace(0,2*pi,numvals); % Covers entire sphere
r = linspace(0,2,numvals); % diamater is 4cm
outputMerged = zeros(numvals,numvals);

[meshR, meshTheta] = meshgrid(r, theta);

fig1 = figure();
sgtitle('Contribution to Concentration for each Nonzero Index of N')
for n = 0:5
% Output Calculation for u(r,theta)
figindex = n + 1;
constant = (2.*n + 1)/(2.^(n+1));
integralIn = @(x) (1/6) .* besselk(0,sqrt((1+x)/2)) .* legendreP(n, x);
integralOut = integral(integralIn,-1,1);
outputTerm = constant .* integralOut .* r.^n;
outputTerm = (legendreP(n,cos(theta))).' * outputTerm;
outputMerged = outputMerged + outputTerm;

% Plotting Mesh Contributions

subplot(3,2,figindex)
surf(meshR .* cos(meshTheta), meshR .* sin(meshTheta), outputTerm, 'edgecolor','none');
xlabel('X Position (cm)'),ylabel('Y Position (cm)'),title(['n = ' num2str(n)])
colormap('jet');
bar = colorbar;
caxis([-0.1, 0.1])
ylabel(bar, 'Quantity Concentration'); 
view(0,90);
end

%% Merged Mesh Figure
fig2 = figure();
surf(meshR .* cos(meshTheta), meshR .* sin(meshTheta), outputMerged, 'edgecolor','none');
xlabel('X Position (cm)'),ylabel('Y Position (cm)'),title({'Quantity Concentration of a Uniform Slice of the Sphere' 'Merged Results of First 6 Nonzero Coefficients'})
colormap('jet');
bar = colorbar;
ylabel(bar, 'Quantity Concentration');
view(0,90);
%% Export Figures
exportgraphics(fig1, fullfile(filepath, 'Independent Output.jpg'), 'resolution', 600);
exportgraphics(fig2, fullfile(filepath, 'Merged Output.jpg'), 'resolution', 600);