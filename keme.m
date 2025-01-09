function [ke,me] = keme(Le,rho,c)
%basic_FEM_acoustic
%ke me is the function for obatianing the ke and me matrix for wave propogation solution Summary of this function goes here
%   Detailed explanation goes here
ke = 1/rho*[1/Le -1/Le
            -1/Le 1/Le];
me = 1/(rho*c^2)*[Le/3 Le/6
                   Le/6 Le/3];
end