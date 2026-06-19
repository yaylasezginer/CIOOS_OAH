function savedir = compile_SSC_sensor_network(siteID, downloadsDir)
% 
% Compile Salish Sea Cast virtual mooring files downloaded from Salish Sea
% Cast ERDDAP. Compile into a single .nc file
% 
% INPUTS
% siteID - virtual mooring siteID - matches naming convention for
% downloaded files. - string
% downlaodsDir - file directory for ERDDAP downloaded datasets
%
% OUTPUTS 
% single compiled .nc file saved to child folder 'compiled' within
% downloads parent directory

df = dir(downloadsDir);
files = {df.name};
site_fn = files(contains(files, siteID));

var = [{'total_alkalinity'}, {'dissolved_inorganic_carbon'}, {'dissolved_oxygen'},{'temperature'},{'salinity'}];

for v = 1:numel(var)
    
    var_keep = contains(site_fn, var{v});
    fn = site_fn(var_keep);
    savedir = [downloadsDir '/compiled/' siteID '2015_2025h_' var{v} '.nc'];
    
    cat_mat = [];
    time = [];
    tic
    for i = 1:numel(fn)
        filename = [downloadsDir '/' fn{i}];
        data = ncread(filename, var{v});
        time = [time; ncread(filename, 'time')];
        cat_mat = cat(4, cat_mat, data); % concatenate along 4th dim
    end
    depth = ncread(filename, 'depth');
    [x,y,z,t] = size(cat_mat);

    %Write concatenated data to netcdf file
    nccreate(savedir, var{v}, 'dimensions',{'x',x,'y',y,'z',z,'t',t},'FillValue','disable')
    nccreate(savedir, 'time', 'dimensions',{'t',t},'FillValue', 'disable'); 
    ncwrite(savedir, 'time', time);
    nccreate(savedir, 'depth', 'dimensions',{'z',z}, 'FillValue','disable'); 
    ncwrite(savedir, 'depth', depth);
    ncwrite(savedir, var{v}, cat_mat);
    fprintf('-> %s (%.1f secs)\n',savedir,toc);
end