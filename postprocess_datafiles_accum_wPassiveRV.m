%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_refPTO_accum_wPassiveRV")
        load(files(j).name)

        q_permMean_array(iVar) = q_permMean;
        PP_WEC_array(iVar) = PP_WEC;
        PP_wp_array(iVar) = PP_wp;
        PP_rv_array(iVar) = PP_rv;
        PP_hPRV_array(iVar) = PP_hPRV;
        PP_roPRV_array(iVar) = PP_roPRV;
        dpdt_max_array(iVar) = dpdt_max;
        dpdt_97_array(iVar) = dpdt_97;

    end

end

%% Find indices for missing data files
if 1

files = dir;
nfiles = size(files,1);
notDone = 1:675;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_refPTO_accum_wPassiveRV")
        load(files(j).name)
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        
    end

end
jobArrayStr = num2str(notDone(1));
for j = 2:length(notDone)
    jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
end

if 1
for j = 1:length(notDone)
    iVar = notDone(j);
    q_permMean_array(iVar) = nan;
    PP_WEC_array(iVar) = nan;
    PP_wp_array(iVar) = nan;
    PP_rv_array(iVar) = nan;
    PP_hPRV_array(iVar) = nan;
    PP_roPRV_array(iVar) = nan;
    dpdt_max_array(iVar) = nan;
    dpdt_97_array(iVar) = nan;

    Vtotal_missing(j) = Vtotal_mesh(iVar);
    X_missing(j) = X_mesh(iVar);
    kv_missing(j) = kv_mesh(iVar);
    
end
end

end

% Display study variable values with missing data
tbl = table(notDone',1e3*Vtotal_missing',X_missing',kv_missing','VariableNames', {'iVar', 'Total Volume (L)','portion at RO inlet', 'kv'});
disp(tbl);

%% Transform data to 3D variable mesh
I = length(kv);
J = length(Vtotal);
K = length(X);


test = 1;
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*K*(i-1) + K*(j-1) + k;

            q_permMean_3D(i,j,k) = q_permMean;
            PP_WEC_3D(i,j,k) = PP_WEC_array(m);
            PP_wp_3D(i,j,k) = PP_wp_array(m);
            PP_rv_3D(i,j,k) = PP_rv_array(m);
            PP_hPRV_3D(i,j,k) = PP_hPRV_array(m);
            PP_roPRV_3D(i,j,k) = PP_roPRV_array(m);
            dpdt_max_3D(i,j,k) = dpdt_max_array(m);
            dpdt_97_3D(i,j,k) = dpdt_97_array(m);

            Vtotal_3D(i,j,k) = Vtotal_mesh(m);
            X_3D(i,j,k) = X_mesh(m);
            kv_3D(i,j,k) = kv_mesh(m);
            test = (Vtotal_mesh(m) == Vtotal(j)) ...
                && (X_mesh(m) == X(k)) ...
                && (kv_mesh(m) == kv(i)) ...
                && test;
        end
    end
end

if ~test; error('indexing incorrect'); end
clearvars test


%% Plot average power loss as a function of total accumulator volume for distribution (color) and valve coefficient (line type)
 % select indices to plot
    % distribution
    iiK = 1:4:K;
    nK = length(iiK);
    % valve coeff.
    iiI = 1:2:I;
    nI = length(iiI);
   
  % selct variable to plot
  Y = PP_rv_3D;

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

bottomEdge = 1;
leftEdge = 3;
width = 3.5625; % one column: 3+9/16, two column: 7.5
height = 3.25;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

% dummy plots for legend
iLeg = 0;

 % color - total volume
for k = 1:nK
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['fraction at RO inlet = ',num2str(X(iiK(k)))]);
end

for i = 1:nI
    plot(-99*[1, 0.5],-99*[1, 0.5],'k','LineStyle', linestyles{i});
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['k_v = ',num2str(kv(iiI(i))*sqrt(1000)),'(L/s/kPa^{1/2})']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vtotal*1e3,1e-3*Y(iiI(i),:,iiK(k)), ...
            'LineStyle', linestyles{i}, ...
            'Color',color(k,:), ...
            'LineWidth',1);
        p(k,i).HandleVisibility='off';
    end
end


xlabel('volume (L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Mean Power Loss With Passive Ripple Control'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xlim([0 1e3*max(Vtotal)])
ylim([0 1e-3*max(Y(:))])

%% Plot rate of change as a function of total accumulator volume for distribution (color) and valve coefficient (line type)
 % select indices to plot
    % distribution
    iiK = 1:4:K;
    nK = length(iiK);
    % valve coeff.
    iiI = 1:2:I;
    nI = length(iiI);
   
  % selct variable to plot
  Y = dpdt_max_3D;
  % Y = dpdt_97_3D;

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

bottomEdge = 1;
leftEdge = 3;
width = 3.5625; % one column: 3+9/16, two column: 7.5
height = 3.25;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

% dummy plots for legend
iLeg = 0;

 % color - total volume
for k = 1:nK
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['fraction at RO inlet = ',num2str(X(iiK(k)))]);
end

for i = 1:nI
    plot(-99*[1, 0.5],-99*[1, 0.5],'k','LineStyle', linestyles{i});
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['k_v = ',num2str(kv(iiI(i))*sqrt(1000)),'(L/s/kPa^{1/2})']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vtotal*1e3,1e-3*Y(iiI(i),:,iiK(k)), ...
            'LineStyle', linestyles{i}, ...
            'Color',color(k,:), ...
            'LineWidth',1);
        p(k,i).HandleVisibility='off';
    end
end

% plot target limit
iLeg = iLeg+1;;
plot(Vtotal([1 end])*1e3,1e-3*par.control.dpdt_ROmax*[1 1])
legLabels(iLeg) = convertCharsToStrings( ...
        ['target limit']);

xlabel('volume (L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('rate of change in pressure (kPa/s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Maximum Rate of Change in Pressure With Passive Ripple Control'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xlim([0 1e3*max(Vtotal)])
% ylim([0 1e-3*max(Y(:))])