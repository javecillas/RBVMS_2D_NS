%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Residual-based Variational Multi-scale Formulation of
%%%%%% 2D Incompressible Navier-Stokes Equation
%%%%%% Javier A. Avecillas
%%%%%% Department of Civil and Environmental Engineering, UIUC
%%%%%% Version 0.1 - Last update: 12/14/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

% Read mesh
pts = readmatrix('../mesh/p.txt');
con = readmatrix('../mesh/t.txt');
bcs = readmatrix('../mesh/bcs.txt');

%% Create material
ro      = 1.0;
mu      = 0.02;
kappa   = 1.0;
tf      = 0.01;
nstp    = 2;
isSUPG  = true;
mat     = material(1,ro,mu,kappa,tf,nstp,isSUPG);

%% Velocity
% Create mesh for velocity
meshV = mesh();
meshV.setNodes(pts);
% Create tri3 linear triangular elements
intRule = 3;
meshV.setTri3(mat,con,intRule);
% Set velocity and pressure global DoF numbering
meshV.setDoFVP();
% Read boundary conditions for velocity
meshV.readBC(bcs);
meshV.setBC();

%% Analysis
% Create analysis
myAnal = nonlintrans(meshV,mat);

%% Print mesh to VKT
%
if ~exist('../paraview','dir')
    mkdir('../paraview');
end
%
myAnal.print2VTK('../paraview/sol_step_0.vtk',myAnal.Uh);

%% Solve linear system
myAnal.solve();
