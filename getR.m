function [R, material, K_] = getR(K)
% Function returns the five strings that are closest in stiffness to the
% desired normalized stiffness
T = readtable('string_database.xlsx');
string_material = T.Material;
string_radius = T.Radius_mm_;
string_norm_stiffness = T.NormalizedStiffness_N_;

a = abs(K-string_norm_stiffness);
[~,I] = mink(a,5);

R = string_radius(I);
K_ = string_norm_stiffness(I);
material = string_material{I};
