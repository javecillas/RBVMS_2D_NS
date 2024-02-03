%% Javier A. Avecillas-Leon
%  CEE 598
%  Fall 2023

%%
classdef material < handle
    %MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % ID
        ID = -1;
        % Density
        ro = -1;
        % Convection constant
        mu = -1;
        % Diffusion constant
        kappa = -1;
        % Constant
        CL = 1.0;

        % Simulation time
        tf = 1;
        % Number of steps
        nstp = 1;
        % Time increment
        dt;

        % SUPG
        isSUPG = true;
    end
    
    methods
        % Constructor
        function obj = material(iD,ro,mu,kappa,tf,nstp,SUPG)
            obj.ID = iD;
            obj.ro = ro;
            obj.mu = mu;
            obj.kappa = kappa;
            
            % Time parameters
            obj.tf   = tf;
            obj.nstp = nstp;
            obj.dt = tf/nstp;

            % SUPG
            obj.isSUPG = SUPG;
        end
        %
    end
end

