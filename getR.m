function [R, material, K_] = getR(K)
% Function returns the five strings that are closest in stiffness to the
% desired normalized stiffness
T = readtable('string_database.xlsx');
string_material = T.Material;
string_radius = T.Radius_mm_;
string_norm_stiffness = T.NormalizedStiffness_N_;

a = abs(K-string_norm_stiffness);
[~,I] = mink(a,3);

R = string_radius(I);
K_ = string_norm_stiffness(I);
material = cell(3,1);
for i = 1:length(I)
    material{i} = string_material{I(i)};
end
end