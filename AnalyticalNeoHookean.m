% A MATLAB code to compute Stress tensor for neo-Hookean material given Deformation gradient tensor
% by Sina Taghizadeh
% November 2022
author = "Sina Taghizadeh";
license = "GNU GPL Version 3.0";

% Deformation gradient tensor
%F = [1,0,0;0,1,0.4;0,0,1]; %for simple shear
F = [1,0,0;0,1,0;0,0,1.2]; %for simple extension

% Material parameters
shear_modulus = 1;
bulk_modulus = 0.78;

% Calculate Stress
C1=shear_modulus / 2;
D1= 2 / bulk_modulus;

j=det(F); %Determinant of the deformation gradient tensor
FBar=F*j^(-1/3); %the distortion gradient
BBar=FBar*transpose(FBar); %the deviatoric left Cauchy-Green deformation tensor
Sigma_SimpleExtension=2/j*C1*(BBar-(1/3)*trace(BBar)*eye(3))+2/D1*(j-1)*eye(3)
%Sigma_SimpleShear=2/j*C1*(BBar-(1/3)*trace(BBar)*eye(3))+2/D1*(j-1)*eye(3)

