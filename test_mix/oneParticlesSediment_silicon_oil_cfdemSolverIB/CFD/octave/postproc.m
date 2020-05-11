% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
particle_vel_1 = load('../../DEM/post/velocity_particle_1.txt');
linienstaerke = 1;
MarkerGroesse = 4;

figure(2)
h = plot(particle_vel_1(:, 1), particle_vel_1(:, 4), '*');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 1.4 -0.14 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Comparison of the settling velocity of two particles', 'FontSize', 12)
legend('following particle', 'leading particle')
set(gca, 'FontSize', 12)
print('vel_y_two_particles.png')
