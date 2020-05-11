% postproc.m: loads the data necessary and saves the graphs

clear;
clc;
close all;

% read data from simulation
vel_IB = load('../oneParticlesSediment_silicon_oil_cfdemSolverIB/DEM/post/velocity_particle_1.txt');
vel_mix = load('../oneParticlesSediment_silicon_oil_cfdemSolverMix/DEM/post/velocity_particle_1.txt');
linienstaerke = 1;
MarkerGroesse = 4;

figure(2)
h = plot(vel_IB(:, 1), vel_IB(:, 4), '*',
         vel_mix(:, 1), vel_mix(:, 4), '+');
set(h,'LineWidth', linienstaerke, 'MarkerSize', MarkerGroesse);
set(gca, 'FontSize', 12)
axis([0.0 1.4 -0.14 0.0])
xlabel('time (s)')
ylabel('z-veloctiy (m/s)')
title('Comparison of the settling velocity of two solvers', 'FontSize', 12)
legend('cfdemSolverIB', 'cfdemSolverMix')
set(gca, 'FontSize', 12)
print('settling_direction_velocity.png')
