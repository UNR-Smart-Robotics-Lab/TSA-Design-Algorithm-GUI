function [R, material, K_] = getR(K, min_max)
% Function returns the five strings that are closest in stiffness to the
% desired normalized stiffness
T = readtable('string_database.xlsx');
string_material = T.Material;
string_radius = T.Radius_mm_;
string_norm_stiffness = T.NormalizedStiffness_N_;
%
a = (K-string_norm_stiffness);
if min_max == 0 % The MINIMUM STIFFNESS
    a = a(a>0);
elseif min_max == 1 % The MAXIMUM STIFFNESS
    a = a(a<0);
end
%
material = cell(3,1);
if ~isempty(a)
    [~,I] = mink(a,3);
    R = string_radius(I);
    K_ = string_norm_stiffness(I);
    
    for i = 1:length(I)
        material{i} = string_material{I(i)};
    end
else
    R = [1,1,1]; % default
    K_ = K*[1,1,1];
    for i = 1:3
        material{i} = 'NA';
    end
end
end