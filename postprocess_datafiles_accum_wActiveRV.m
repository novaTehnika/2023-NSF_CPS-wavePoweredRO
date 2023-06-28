%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_refPTO_accum_wActiveRV")
        load(files(j).name)

        q_permMean_array(iVar) = q_permMean;
        PP_WEC_array(iVar) = PP_WEC;
        PP_wp_array(iVar) = PP_wp;
        PP_rv_array(iVar) = PP_rv;
        PP_hPRV_array(iVar) = PP_hPRV;
        PP_roPRV_array(iVar) = PP_roPRV;
        PP_aPRV_array(iVar) = PP_aPRV;
        PP_bPRV_array(iVar) = PP_bPRV;
        dpdt_max_array(iVar) = dpdt_max;
        dpdt_97_array(iVar) = dpdt_97;

    end

end

%% Find indices for missing data files
if 0

files = dir;
nfiles = size(files,1);
notDone = 1:4560;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_refPTO_accum_wActiveRV")
        load(files(j).name)
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        
    end

end
jobArrayStr = num2str(notDone(1));
for j = 2:length(notDone)
    jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);

end

end
%% Transform data to 3D variable mesh
I = length(Vtotal);
J = length(X);
K = length(kv);

test = 1;
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*I*(k-1) + I*(j-1) + i;

            q_permMean_3D(i,j,k) = q_permMean;
            PP_WEC_3D(i,j,k) = PP_WEC_array(m);
            PP_wp_3D(i,j,k) = PP_wp_array(m);
            PP_rv_3D(i,j,k) = PP_rv_array(m);
            PP_hPRV_3D(i,j,k) = PP_hPRV_array(m);
            PP_roPRV_3D(i,j,k) = PP_roPRV_array(m);
            PP_aPRV_3D(i,j,k) = PP_aPRV_array(m);
            PP_bPRV_3D(i,j,k) = PP_bPRV_array(m);
            dpdt_max_3D(i,j,k) = dpdt_max_array(m);
            dpdt_97_3D(i,j,k) = dpdt_97_array(m);

            Vtotal_3D(i,j,k) = Vtotal_mesh(m);
            X_3D(i,j,k) = X_mesh(m);
            kv_3D(i,j,k) = kv_mesh(m);
            test = (Vtotal_mesh(m) == Vtotal(i)) ...
                && (X_mesh(m) == X(j)) ...
                && (kv_mesh(m) == kv(k)) ...
                && test;
        end
    end
end

if ~test; error('indexing incorrect'); end
clearvars test


%% Plot average power results for one turndown fraction but multiple turndown ratio
black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

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

for k = 1:K
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['min. load frac. = ',num2str(lbFrac(k))]);
end

scatter(-99*[1 1],-99*[1 1],50, ...
        'filled','x','LineWidth',2,'MarkerEdgeColor','k');
iLeg = iLeg+1;
legLabels(iLeg) = convertCharsToStrings('MPC algorithm');

plot(-99*[1 1],-99*[1 1],'-k','LineWidth',1)
iLeg = iLeg+1;
legLabels(iLeg) = convertCharsToStrings('Optimal Coulomb damping');

plot(-99*[1 1],-99*[1 1],'--k','LineWidth',1)
iLeg = iLeg+1;
legLabels(iLeg) = "Fixed Coulomb damping";

% plot real data

for k = 1:K
    s(k) = scatter(1e-6*T_max,1e-3*PP_yearlyAve(:,k),50, ...
        'filled','x','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    s(k).HandleVisibility='off';
    
    p(k) = plot(1e-6*Tmax_Coulomb,1e-3*PP_yearlyAveCoulomb(k,:),'-','Color',color(k,:),'LineWidth',1);
    p(k).HandleVisibility='off';
end
p(K+1) = plot(1e-6*Tmax_Coulomb,1e-3*PP_yearlyAveCoulomb(end,:),'--k','LineWidth',1);
p(K+1).HandleVisibility='off';

xlabel('torque, max (MNm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Mean Power Capture, Yearly Average'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels)
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xlim([0 1e-6*max(T_max)])
