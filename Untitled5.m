close all

binhm = 10;

figure
x = rskinj(:, 1);
y = kskinj(:, 1);
hist3([x,y], [binhm binhm])
xlabel('{\bf r_{skin}} (m)')
ylabel('{\bf k_{skin}} (mD)')
zlabel('{\bf Frequency}')
% title('LAYER 1');
hold on
N = hist3([x,y], [binhm binhm]);
title('LAYER 1')
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
colormap('hot') % Change color scheme 
colorbar % Display colorbar
h.ZData = -2*max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];

% figure
% x = rskinj(:, 3);
% y = kskinj(:, 3);
% hist3([x,y], [binhm binhm])
% xlabel('{\bf r_{skin}} (m)')
% ylabel('{\bf k_{skin}} (mD)')
% zlabel('{\bf Frequency}')
% % title('LAYER 1');
% hold on
% N = hist3([x,y], [binhm binhm]);
% title('LAYER 3')
% N_pcolor = N';
% N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
% xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
% yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
% h = pcolor(xl,yl,N_pcolor);
% colormap('hot') % Change color scheme 
% colorbar % Display colorbar
% h.ZData = -2*max(N_pcolor(:))*ones(size(N_pcolor));
% ax = gca;
% ax.ZTick(ax.ZTick < 0) = [];

figure
x = rskinj(:, 2);
y = kskinj(:, 2);
hist3([x,y], [binhm binhm])
xlabel('{\bf r_{skin}} (m)')
ylabel('{\bf k_{skin}} (mD)')
zlabel('{\bf Frequency}')
% title('LAYER 2');
hold on
N = hist3([x,y], [binhm binhm]);
title('LAYER 2')
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
colormap('hot') % Change color scheme 
colorbar % Display colorbar
h.ZData = -4*max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];

rskinj(rskinj < rw) = rw;
kskinj(kskinj < 0) = 0;
SJ_guardar = zeros(Ne, nlayers);
for i = 1:Ne
    SJ_guardar(i, :) = (kj./kskinj(i, :)-1).*log(rskinj(i, :)./rw);
end

% kj_true = [600 600];
% kskinj_true = [240 600];
% rskinj_true = [0.25 0];
% Sj_true = (kj_true./kskinj_true-1).*log(rskinj_true./rw);
% Sj_true = [3.14 3.14];

binh = 100;
figure
grid on
hold on
for j = 1:nlayers
    histogram(SJ_guardar(:, j), binh)
    histogram(Sj_true(j)*ones(1, Ne), 2)
end
histogram(seq*ones(1, Ne), 2)
xlim([0 max(max(Sj_true)+1, 10)]);
legend({'S_1 (estimated)','S_1 (true)', 'S_2 (estimated)', 'S_2 (true)', 'S_3 (estimated)', 'S_3 (true)', 'S_{eq}'});
xlabel('{\bf Skin factor}', 'FontSize', 15);
ylabel('{\bf Frequency}', 'FontSize', 15);

% figure
% % x = rskinj(:, 3);
% % y = kskinj(:, 3);
% x = SJ_guardar(:, 1);
% y = SJ_guardar(:, 2);
% hist3([x,y], [bin bin])
% % xlabel('{\bf r_{skin}} (m)')
% % ylabel('{\bf k_{skin}} (mD)')
% xlabel('{\bf S_{1}}')
% ylabel('{\bf S_{2}}')
% zlabel('{\bf Frequency}')
% % title('LAYER SKIN HISTOGRAM');
% hold on
% N = hist3([x,y], [bin bin]);
% N_pcolor = N';
% N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
% xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
% yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
% h = pcolor(xl,yl,N_pcolor);
% colormap('hot') % Change color scheme 
% colorbar % Display colorbar
% h.ZData = -2*max(N_pcolor(:))*ones(size(N_pcolor));
% ax = gca;
% ax.ZTick(ax.ZTick < 0) = [];

% b = max(max(max(d)), max(max(d_initial)));
% a = min(min(min(d(11:end, :))), min(min(d_initial(11:end, :))));

b = ceil(log10(max(max(d(tr:end,:))))); b = 10^b;
a = floor(log10(min(min(d(tr:end,:))))); a = 10^a;

figure
loglog(t, d_initial(:, 1), 'Color', [.6 .6 .6])
hold on
loglog(t, d(:, 1), 'y')
loglog(t, pwf, '*-r')
loglog(t, d_initial, 'Color', [.5 .5 .5])
plot(t, d, 'y')
loglog(t, mean(d, 2), 'b', 'LineWidth', 6) 
loglog(t, pwf, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1 0 0])
grid on
xlabel('{\bf time (h)}', 'FontSize', 15)
ylabel('{\bf P (kgf/cm^2)}', 'FontSize', 14)
legend({'Initial Ensemble', 'Final Ensemble', 'Observed data'}, 'FontSize', 12, 'Location', 'southeast');
axis([t(11) t(end) a b]);