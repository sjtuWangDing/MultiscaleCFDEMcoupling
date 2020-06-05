% postproc.m: loads the data necessary and saves the graphs
% (compared to the reults by glowinski)

clear;
clc;
close all;

% read data from simulation
vel = load('../../DEM/post/velocity_particle_1.txt');
linienstaerke = 1;
MarkerGroesse = 4;

figure(2)
h = plot(vel(:, 1), vel(:, 4), '*');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 1.4 -0.14 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Settling velocity using cfdemSolverMix', 'FontSize', 12)
legend('cfdemSolverMix')
set(gca, 'FontSize', 12)
print('settling_direction_velocity.png')
