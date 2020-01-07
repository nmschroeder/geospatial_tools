%% MATLAB function to compute area ratio of region of sphere bounded by given 
% latitude and longitude

% Inputs ---------------------------------------------------------------%

% lat1: the lower latitude boundary of quadrangle

% lon1: the lower longitude boundary of quadrangle

% lat2: the upper latitude boundary of quadrangle

% lon2: the upper longitude boundary of quadrangle

% Outputs -------------------------------------------------------------%

% area_ratio: the ratio of the entire surface area of the Earth represented
% by the quadrangle

function area_ratio = area_sphere(lat1, lon1, lat2, lon2)
        whole_sphere = 4*pi; % assumes unit radius of 1
        % Correct for different direction of measurement for phi (from
        % z-axis) and convert to radians
        phi2 = 90 - lat1; phi2 = phi2*(2*pi)/360;
        phi1 = 90 - lat2; phi1 = phi1*(2*pi)/360;
        theta1 = lon1*(2*pi)/360;
        theta2 = lon2*(2*pi)/360;
        f=@(phi,theta)(sin(phi));
        area_ratio = 1/whole_sphere*integral2(f,phi1,phi2,theta1,theta2);
end
