
function [data_export, data_export_nc, metadata, stats] = ExportFilesCreate(savedir, siteID, siteLat, siteLon, siteDepth, sensorDir, loDir, sscDir, params_sensor, params, units, DataURL, loURL, sscURL)

% read data from downloaded files from sensor and model (SalishSeaCast and LiveOcean) data providers 
% Bring datasets to the same temporal resolution and reformat into a single
% matrix. Return matrix as a csv and nc filetype. 
% Note: This fxn is designed for single depth fixed points
% Note: update L66 with local directory for nctoolbox
%
% INPUTS:
% savedir: folder directory for export files (string)
% siteID: Name of site (string)
% siteLat: Lat in degree decimal (scalar)
% siteLon: Lon in degree decimal (scalar)
% sensorDir: filepath and filename for sensor data (string)
% loDir: filepath and filename for liveocean data (string)
% sscDir: filepath and filename for SalishSeaCast data (string)
% params_sensor: parameter names given by data provider for variables measured by sensor (cell array)
% params: standardized names for parameters measured by sensor. May vary from parameter name given
%   by data provider depending on controlled vocabulary usage (cell array). 
%   possible elements of params: time, temp, sal, dO2, DIC, TA, pCO2, pH
% units: sensor measurement units (cell array equal in size to params)
%   time: possible units: s_1970, s_1950, string format (copy string format
%   exactly), dn (matlab dn)
%   temp: possible units: degC, K
%   sal: possible units: psu, g/kg, ppt
%   dO2: possible units: mmol/m3, umol/kg, umol/L, mL/L, mg/L, %
%   DIC: possible units: mmol/m3, umol/kg, umol/L
%   TA: possible units: mmol/m3, umol/kg, umol/L
%   pCO2: possible units: uatm 
% DataURL: URL for in-situ data source citation page (string)
% loURL: URL for LiveOcean data access/citation page (string)
% sscURL: URL for ssc data access page (string)
%
% OUTPUTS:
% data_export: reformatted data matrix identical to data contents of csv export file. Interpolated values are not included (table)
% data_export_csv: reformatted data matrix identical nc export file
% contents. Includes interpolated values. For CIOOS Pacific internal use
% only. (table)
% metadata: metadata associated with export file data. 
% 
% Data will always be returned in standardized units:
% time: time in UTC
% T: degC
% S: g/kg
% DIC: mmol/m3
% TA:mmol/m3
% dO2: uM
% pH: ()
% pCO2: uatm
% saturation state: ()
% 
% Notes: prior to running, salishseacast data should be compiled into a
% single .nc file


%% ------------------------------------------------------------------------
% Initialize outputs
% -------------------------------------------------------------------------

data_export = table();
data_export_nc = struct();

%% ------------------------------------------------------------------------
% DETERMINE FILE TYPES AND READ DATA
% -------------------------------------------------------------------------

% Prepare to read .nc datafiles
addpath('/Users/yaylasezginer/Documents/MATLAB/nctoolbox') 
setup_nctoolbox

% 1. Read sensor data
[~,~,sensorFileType] = fileparts(sensorDir);

switch sensorFileType
    case '.csv'
        sensorData_x = struct(readtable(sensorDir)); 
    case '.mat'
        sensorData_x = load(sensorDir); 
    case '.nc'
        sensorData_nc = ncdataset(sensorDir); 
        % reformat to a struct
        varSensor = sensorData_nc.variables;
        sensorData_x = struct();
        for v = 1:numel(varSensor)
            sensorData_x.varSensor{v} = squeeze(sensorData_nc.data(varSensor{v}));
        end
end
% Convert sensor data parameter names from given names to standardized
% vocabulary
for v = 1:numel(params_sensor)
    sensorData.(params{v}) = sensorData_x(params_sensor{v});
end

% 2. Read LiveOcean data (always in .nc format). 3D data to be reformatted to
% vector from matrix. If site doesn't overlap with LO model domain, return
% empty vectors

varLo = [{'TIC'},{'oxygen'},{'alkalinity'},{'salt'},{'temp'}]; % native LO variable names
stdNames = [{'DIC'},{'dO2'},{'TA'},{'sal'},{'temp'}]; % controlled vocabulary variable names
loData = struct();
if ~isempty(loDir)
    loData_nc = ncdataset(loDir);
    loData.time = loData_nc.data.('ocean_time'); % s since 1970
    loData.dn = datenum(1970,1,1,0,0,loData.time);
    [~,loDepth_ind] = min(abs(siteDepth - loData_nc.data.('z_rho')));
    for v = 1:numel(varLo)
        loData.(stdNames{v}) = squeeze(loData_nc.data(varLo{v}));
        loData.(stdNames{v}) = loData.(stdNames{v})(loDepth_ind,:)';
    end
else 
    loData.dn = [];
    for v = 1:numel(varLo)
        loData.(stdNames{v}) = [];
    end
end

% 3. Read SSC data (always in .nc format). If site doesn't overlap with LO model domain, return
% empty vectors

varSsc = [{'dissolved_inorganic_carbon'},{'dissolved_oxygen'},{'total_alkalinity'},{'salinity'},{'temperature'}];
sscData = struct();
if ~isempty(sscDir)
    sscData_nc = ncdataset(sscDir);
    sscData.time = sscData_nc.data.('time'); % s since 1970
    sscData.dn = datenum(1970,1,1,0,0,sscData.time);
    [~,sscDepth_ind] = min(abs(siteDepth - sscData_nc.data.('depth')));
    for v = 1:numel(varSsc)
        sscData.(stdNames{v}) = squeeze(sscData_nc.data(varSsc{v}));
        sscData.(stdNames{v}) = sscData.(stdNames{v})(sscDepth_ind,:)';
    end
else
    sscData.dn = [];
    for v = 1:numel(varSsc)
        sscData.(stdNames{v}) = [];
    end
end

%% ------------------------------------------------------------------------
% Ensure sensor data is in standardized units
% -------------------------------------------------------------------------

addpath /Users/yaylasezginer/Documents/MATLAB/GSW-Matlab-master/Toolbox
sensorData.pres = sw_pres(siteDepth, siteLat); % units: db
sensorData.dens = sw_dens(sensorData.sal, sensorData.temp, sensorData.pres); % kg/m3

units_out = cell(size(units));

% Time

time_ind = strcmp(params, 'time');
if ~isempty(time_ind)
    time_unit = units{time_ind};
    switch time_unit
        case 's_1970'
            sensorData.dn = datenum(1970, 1,1,0,0,sensorData.time);
        case 'dn'
            sensorData.dn = sensorData.time;
        case 's_1950'
            sensorData.dn = datenum(1950, 1,1,0,0,sensorData.time);
        case 'string'
            sensorData.dn = datenum(sensorData.time, 'string'); % replace string with formatIn
    end
    units_out{time_ind} = 'UTC';
end


% Temperature
temp_ind = strcmp(params,'temp');
if ~isempty(temp_ind)
    temp_unit = units{temp_ind};
    switch temp_unit
        case 'degC'
            sensorData_std.temp = sensorData.temp;
        case 'K'
            sensorData_std.temp = sensorData.temp - 273.15;
    end
    units_out{temp_ind} = 'degC';
end

% Salinity
sal_ind = strcmp(params,'sal');
if ~isempty(sal_ind)
    sal_unit = units{sal_ind};
    switch sal_unit
        case 'psu'
            sensorData_std.sal = gsw_SP_from_SA(sensorData.sal,sensorData.pres,siteLon,siteLat);
        case 'g/kg'
            sensorData_std.sal = sensorData.sal;
        case 'ppt'
            sensorData_std.sal = sensorData.sal;
    end
    units_out{sal_ind} = 'g/kg';
end

% Dissolved O2
dO2_ind = strcmp(params,'dO2');
if ~isempty(dO2_ind)
    dO2_unit = units{dO2_ind};
    switch dO2_unit
        case 'mmol/m3'
            sensorData_std.dO2 = sensorData.dO2;
        case 'umol/L'
            sensorData_std.dO2 = sensorData.dO2;
        case 'umol/kg'
            sensorData_std.dO2 = sensorData.dO2 .* 1000 ./ sensorData.dens;
        case 'mL/L'
            sensorData_std.dO2 = 44.6596 .* sensorData.dO2;
        case 'mg/L'
            sensorData_std.dO2 = 31.25 .* sensorData.dO2;
        case 'mg/kg'
            sensorData_std.dO2 = (sensorData.dO2 .* sensorData.dens)./31.998;
        case 'mL/kg'
            sensorData_std.dO2 = (sensorData.dO2 .* sensorData.dens)./22.392;
    end
    units_out{dO2_ind} = 'mmol/m3';
end

% DIC
DIC_ind = strcmp(params, 'DIC');
if ~isempty(DIC_ind)
    DIC_unit = units{DIC_ind};
    switch DIC_unit
        case 'mmol/m3'
            sensorData_std.DIC = sensorData.DIC;
        case 'umol/kg'
            sensorData_std.DIC = sensorData.DIC .* 1000 ./ sensorData.dens;
        case 'umol/L'
            sensorData_std.DIC = sensorData.DIC;
    end
    units_out{DIC_ind} = 'mmol/m3';
end

% TA
TA_ind = strcmp(params, 'TA');
if ~isempty(TA_ind)
    TA_unit = units{TA_ind};
    switch TA_unit
        case 'mmol/m3'
            sensorData_std.TA = sensorData.TA;
        case 'umol/kg'
            sensorData_std.TA = sensorData.TA .* 1000 ./ sensorData.dens;
        case 'umol/L'
            sensorData_std.TA = sensorData.TA;
    end
    units_out{TA_ind} = 'mmol/m3';
end


%% ------------------------------------------------------------------------
% Calculate carbonate system parameters for LiveOcean and SalishSeaCast
% -------------------------------------------------------------------------

% Sal input to CO2SYS in psu
% TA and DIC input to CO2SYS in umol/kgSw

% LiveOcean carbonate system 
if ~isempty(loData.dn)
    loData.pres = sw_pres(loData.z_rho,siteLat); % units: db
    loData.dens = sw_dens(loData.salt, loData.temp, loData.pres); % kg/m3
    loData.sal_psu = gsw_SA_from_SP(loData.sal, loData.pres, siteLon, siteLat);
    [DATA,~,~]=CO2SYS(loData.TA*1000./loData.dens,loData.DIC.*1000./loData.dens, ...
        1,2,loData.salt_psu,loData.temp,loData.temp,loData.pres,loData.pres,0,0,0,0,1,10,1,2,2);
    loData.pCO2 = DATA(:,4);
    loData.pH = DATA(:,21);
    loData.omega_arag = DATA(:,18);
    loData.omega_calc = DATA(:,17);
else
    loData.pCO2 = [];
    loData.pH = [];
    loData.omega_arag = [];
    loData.omega_calc = [];
end

% SSC carbonate system 

if ~isempty(sscData.dn)
    sscData.pres = sw_pres(sscData.depth, siteLat); % units: db
    sscData.dens = sw_dens(sscData.S, sscData.T, sscData.pres); % kg/m3
    sscData.sal_psu = gsw_SA_from_SP(sscData.salt, sscData.pres, siteLon, siteLat);
    [DATA,~,~]=CO2SYS(sscData.TA*1000./sscData.dens,sscData.DIC.*1000./sscData.dens, ...
        1,2,sscData.sal_psu,sscData.temp,sscData.temp,sscData.pres,sscData.pres,0,0,0,0,1,10,1,2,2);
    sscData.pCO2 = DATA(:,4);
    sscData.pH = DATA(:,21);
    sscData.omega_arag = DATA(:,18);
    sscData.omega_calc = DATA(:,17);
else
    sscData.pCO2 = [];
    sscData.pH = [];
    sscData.omega_arag = [];
    sscData.omega_calc = [];
end

%% ------------------------------------------------------------------------
% Bring data to consistent temporal resolution avoiding gaps
% -------------------------------------------------------------------------

compiled_dn = unique([sscData.dn; loData.dn; sensorData_std.dn]);
data_export.time = datestr(compiled_dn);
data_export_nc.time = (compiled_dn - datenum(1970,1,1))./datenum(0,0,0,0,0,1); % convert back to s since 1970

[~,iSSC, iCompiledxSSC] = intersect(sscData.dn, compiled_dn);
[~,iLO, iCompiledxLO] = intersect(loData.dn, compiled_dn);
[~,iSensor, iCompiledxSensor] = intersect(sensorData_std.dn, compiled_dn);

% Adjustable interpolation thresholds 

maxGapDays = 24/24; %  1 hr = 1/24 day (adjustable)
maxMissing = 0.3; % (1-maxMissing) = maximum allowable fraction of data to be interpolated

for i = 1:numel(params)
    if strcmp(params{i}, 'time')
        continue
    end
    % create uniform sized arrays 
    data_export.([params{i} '_sensor']) = nan .* ones(size(data_export.time_UTC));
    data_export.([params{i} '_SSC']) = nan .* ones(size(data_export.time_UTC));
    data_export.([params{i} '_LiveOcean']) = nan .* ones(size(data_export.time_UTC));
    
    % Fill in arrays with data corresponding with correct datetime stamp
    % for each datastream 
    data_export.([params{i} '_LiveOcean'])(iCompiledxLO) = loData.(params{i})(iLO);
    data_export.([params{i} '_SSC'])(iCompiledxSSC) = loData.(params{i})(iSSC);
    data_export.([params{i} '_sensor'])(iCompiledxSensor) = sensorData_std.(params{i})(iSensor);

    data_export_nc.([params{i} '_LiveOcean']) = data_export.([params{i} '_LiveOcean']);
    data_export_nc.([params{i} '_SSC']) = data_export.([params{i} '_SSC']);

    % Interpolate high resolution sensor data to correspond with exact LO
    % and SSC time stamps - but avoid interpolating across large data gaps
    % (for example due to sensor servicing etc...). Only .nc export file
    % updated with modified, interpolated data (for internal CIOOS use
    % only)

    % Identify valid (non-NaN) indices
    valid = find(~isnan(data_export.([params{i} '_sensor'])));
    if any(valid)

        gapInd = find(diff(compiled_dn(valid)) >= maxGapDays);
        gapInd = [gapInd; numel(valid)];
        % dataset indices from first valid (non-nan) datapoint to last valid data pt before gap
        ind = valid(1):valid(gapInd(1)-1);

        for chunk = 1:numel(gapInd)-1
            % Ensure there are enough non-nan points to linear interpolate
            if numel(ind) > 3 & sum(~isnan(data_export.([params{i} '_sensor'])(ind))) > maxMissing*numel(ind)
                % fill nans within chunk between gaps
                data_export_nc.([params{i} '_sensor'])(ind) = fillmissing(data_export.([params{i} '_sensor'])(ind),'linear');
                % update indices to next chunk
                ind = valid(gapInd(chunk)+1):valid(gapInd(chunk+1)-1);
            else
                % update indices to next chunk, don't attempt to fill nans 
                ind = valid(gapInd(chunk)+1):valid(gapInd(chunk+1)-1);
            end
        end

    else
        error('No valid sensor data detected for interpolation')
    end

%% ------------------------------------------------------------------------
% Calculating Model Evaluation Stats & Plotting/Saving figures
% -------------------------------------------------------------------------
    obs = data_export_nc.([params{i} '_sensor']);
    stats.LiveOcean.Residuals.(params{i}) = obs(iCompiledxLO) - data_export_nc.([params{i} '_LiveOcean'])(iCompiledxLO);
    stats.SSC.Residuals.(params{i}) = obs(iCompiledxSSC) - data_export_nc.([params{i} '_SSC'])(iCompiledxSSC);
    stats.LiveOcean.RMSE.(params{i}) = sqrt(sum((stats.LiveOcean.Residuals.(params{i})).^2)./numel(iCompiledxLO));
    stats.SSC.RMSE.(params{i}) = sqrt(sum((stats.SSC.Residuals.(params{i})).^2)./numel(iCompiledxSSC));

    sp = figure;
    sp(1) = subplot(3,2,1);
    plot(compiled_dn, data_export.([params{i} '_sensor']),'k.'); hold on
    plot(compiled_dn, data_export.([params{i} '_LiveOcean']),'r.')
    ylabel([params{i} ' ' units{i}])
    legend('Sensor', 'LiveOcean')
    datetick('x')
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)

    sp(2) = subplot(3,2,2);
    plot(compiled_dn, data_export.([params{i} '_sensor']),'k.'); hold on
    plot(compiled_dn, data_export.([params{i} '_SSC']),'r.')
    ylabel([params{i} ' ' units{i}])
    legend('Sensor','SSC')
    datetick('x')
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)
   
    sp(3) = subplot(3,2,3);
    plot(compiled_dn, zeros(numel(compiled_dn),1),'k--'); hold on
    plot(compiled_dn, stats.LiveOcean.Residuals.(params{i}),'r.')
    ylabel(['Residual ' units{i}])
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)
    datetick('x')
  
    sp(4) = subplot(3,2,4);
    plot(compiled_dn, zeros(numel(compiled_dn),1),'k--'); hold on
    plot(compiled_dn, stats.SSC.Residuals.(params{i}),'r.')
    ylabel(['Residual ' units{i}])
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)
    datetick('x') 

    sp(5) = subplot(3,2,5);
    plot([min(obs) max(obs)], [min(obs) max(obs)],'k--'); hold on
    binnedscatter(obs(iCompiledxLO), data_export_nc.([params{i} '_LiveOcean'])(iCompiledxLO))
    ylabel(['LiveOcean ' units{i}])
    xlabel(['Sensor ' units{i}])
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)
    
    sp(6) = subplot(3,2,6);
    plot([min(obs) max(obs)], [min(obs) max(obs)],'k--'); hold on
    binnedscatter(obs(iCompiledxSSC), data_export_nc.([params{i} '_SSC'])(iCompiledxSSC))
    ylabel(['SSC ' units{i}])
    xlabel(['Sensor ' units{i}])
    set(gca,'xgrid','on','ygrid',on,'FontSize',15)  
    sgtitle([siteID ' Lat: ' siteLat, ' Lon: ' siteLon, 'Depth: ' siteDepth])
end

saveas(sp,[savedir '/figures/' siteID '/' params{i} '.png'])

%% ------------------------------------------------------------------------
% Format metadata & save export files
% -------------------------------------------------------------------------

fileSavePath = [savedir '/' siteID];

metadata = cell(10,2+numel(params));
metadata(1:end-1,1) = [{'Site_ID'}; {'Latitude'};{'Longitude'};{'Depth'};{'InSitu_DataSource'},...
    {'LiveOceanModel'},{'SalishSeaCastModel'},{'Variables'};{'Units'}];
metadata{1,2} = siteID; metadata{2,2} = siteLat; metadata{3,2} = siteLon; metadata{4,2} = siteDepth;
metadata{5,2} = DataURL;
metadata{6,2} = loURL;
metadata{7,2} = sscURL;
metadata(8,2:end) = params;
metadata(9,2:end) = units_out;

% Export CSV file for external download
writetable(cell2table(metadata),[fileSavePath '.csv'],'WriteVariableNames',false)
writetable(data_export, [fileSavePath '.csv'], 'WriteMode',"append",'WriteVariableNames',true)

T = struct2table(data_export_nc);
nonan = sum(~isnan(table2array(T(:,2:end))),2) > numel(var_liveocean); % keep rows if there is more data than just the sensor data

% Formatting for internal .nc files Global attributes

ncwriteatt([fileSavePath '.nc'],'/','Site_ID',siteID);
ncwriteatt([fileSavePath '.nc'],'/','Latitude',siteLat);
ncwriteatt([fileSavePath '.nc'],'/','Longitude',siteLon);
ncwriteatt([fileSavePath '.nc'],'/','Depth',siteDepth);

nccreate([fileSavePath '.nc'], 'time','dimensions',{'t',numel(data_export_nc.time)})
ncwrite([fileSavePath '.nc'], 'time',data_export_nc.time)
ncwriteatt([fileSavePath '.nc'],'time','units','s since 1970,01,01 00:00 (UTC)')

for i = 1:numel(params)
    if strcmp(params{i}, 'time')
        continue
    end

    nccreate([fileSavePath '.nc'],[params{i} '_sensor'],'dimensions',{'x',numel(data_export_nc.time)},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_sensor'],data_export_nc.([params{i} '_sensor']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_sensor'], 'source', DataURL)

    nccreate([fileSavePath '.nc'],[params{i} '_SSC'],'dimensions',{'x',numel(data_export_nc.time)},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_SSC'],data_export_nc.([params{i} '_SSC']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_SSC'], 'source', sscURL)    

    nccreate([fileSavePath '.nc'],[params{i} '_LiveOcean'],'dimensions',{'x',numel(data_export_nc.time)},'FillValue','disable')
    ncwrite([fileSavePath '.nc'],[params{i} '_LiveOcean'],data_export_nc.([params{i} '_LiveOcean']))
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'units', units_out{i})
    ncwriteatt([fileSavePath '.nc'],[params{i} '_LiveOcean'], 'source', loURL)    

end

end