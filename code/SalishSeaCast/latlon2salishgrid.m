
function [x,y] = latlon2salishgrid(lat,lon)

% Convert latitude and longitude to Salish Sea grid coordinates.
% Lon and Lat may be arrays of the same size. returned x and y will be the
% same size

home = fileparts(which(mfilename));
load([home '/Bathymetry.mat'])

lat_grid = ubcSSnBathymetryV21x2d08.latitude;
lon_grid = ubcSSnBathymetryV21x2d08.longitude;

[m,n] = size(lat_grid);

% Initialize arrays

x = nan .* ones(size(lon));
y = nan .* ones(size(lat));

% Calculate the distance between the center of each grid cell and the given
% lat x lon coordinates

for ii = 1:numel(x)

    for i = 1:m
        for j = 1:n
            d(i,j) = pdist2([lon_grid(i,j), lat_grid(i,j)], [lon(ii),lat(ii)],'euclidean');
        end
    end

    % Find the grid cell with the shortest distance between its center and the
    % coordinates

    % Note: Salish Sea Cast grid y refers to rows and x refers to columns
    [y(ii),x(ii)] = find(abs(d) == min(abs(d),[],'all'));

end
end