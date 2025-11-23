% run_moving_proton_AMR_fixed.m
% ========================================================================
% Moving proton: uniform vs AMR vs analytical comparison (fixed & improved)
% - Single-file script with local helper functions at the end
% - Animates proton moving horizontally, performs AMR per timestep from uniform base,
% visualizes uniform / AMR / analytical surfaces, cross-section (following
% proton) with node markers, logs errors and node counts.
% - Uses log color scale for surfaces to avoid saturation.
% - Uses correct rho = - \nabla^2 phi_ana for consistency (positive at center).
% - Resets AMR mesh each timestep to adapt to current position.
% - Added AMR mesh visualization with proton circle and annotations.
% - Normalized units for numerical stability (coef=1, L=1, R0=0.01).
% - Changed to top-view surfaces for "2D color surface".
% - Increased speed for noticeable movement without bouncing.
% - Removed unnecessary bouncing since horizontal movement stays within bounds.
% ========================================================================
clear; close all; clc;
%% ------------------------- physical / domain params --------------------
L = 1; % domain size (normalized)
R0 = 0.01; % proton "radius" regularization (normalized)
coef = 1; % Coulomb prefactor (normalized, set 1/(4*pi*eps0)=1)
% Analytical potential (regularized Coulomb)
phi_ana = @(x,y,xc,yc) coef ./ sqrt((x-xc).^2 + (y-yc).^2 + R0^2);
%% ------------------------- simulation params ---------------------------
Nx = 100; Ny = 100; % uniform mesh resolution
MAX_AMR_ITER = 6; % AMR iterations per timestep
THETA = 0.30; % Dörfler marking fraction
tol_eta_stop = 0.12; % stop AMR early when normalized estimator L2 < tol
Nt = 40; % number of time steps (frames)
plot_every = 1; % plot frequency (1 => every step)
% Proton initial position and velocity (classical, horizontal)
proton_pos = [0.2*L, 0.5*L];
speed = 0.6 * L; % total displacement over simulation ~0.6 L
theta_move = 0; % direction angle (radians) (0 -> increasing x)
vx = speed * cos(theta_move);
vy = 0; % horizontal movement
% Tolerances / solver regularization
reg_rel = 1e-12; % relative regularization added to diagonal if needed
%% ----------------------- meshes initialization -------------------------
[nodes_uni, elem_uni] = generate_initial_mesh(L, Nx, Ny);
N_uni = size(nodes_uni,1);
%% ----------------------- bookkeeping arrays ---------------------------
nodecount_amr_history = zeros(Nt,1);
L2_history_uni = nan(Nt,1);
L2_history_amr = nan(Nt,1);
Linf_history_uni = nan(Nt,1);
Linf_history_amr = nan(Nt,1);
proton_traj = zeros(Nt,2);
proton_vel = zeros(Nt,2);
%% ----------------------- cross-section probe (follows proton) ----------
probe_N = 500;
%% ----------------------- figure setup (single figure, tiled) ----------
hFig = figure('Name','Moving proton: Uniform vs AMR vs Analytical','Units','normalized','Position',[0.05 0.06 0.9 0.8]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
% Precreate axes handles for performance; we'll update data each timestep
axUniform3D = nexttile(1); title(axUniform3D,'Uniform (numeric)'); axis(axUniform3D,'equal'); view(axUniform3D,0,90);
axAMR3D = nexttile(2); title(axAMR3D,'AMR (numeric)'); axis(axAMR3D,'equal'); view(axAMR3D,0,90);
axAnalytic3D= nexttile(3); title(axAnalytic3D,'Analytical (regularized)'); axis(axAnalytic3D,'equal'); view(axAnalytic3D,0,90);
axCrossPot = nexttile(4); title(axCrossPot,'Cross-section: potentials'); xlabel(axCrossPot,'x'); ylabel(axCrossPot,'\phi');
axCrossErr = nexttile(5); title(axCrossErr,'Cross-section: pointwise abs error (log)'); xlabel(axCrossErr,'x'); ylabel(axCrossErr,'|error|');
axMesh = nexttile(6); title(axMesh,'AMR Mesh'); axis(axMesh,'equal');
% Color limits for log scaling (no clipping, use log)
min_log_phi = -6; % floor for log scale
max_phi = coef / R0; % max at center
%% ----------------------- main time-stepping loop ----------------------
for it = 1:Nt
    % update proton classical motion (no bounce needed, stays in [0.2L, 0.8L])
    if it == 1
        % initial proton_pos already set above
    else
        proton_pos = proton_pos + [vx, vy] * (1/(Nt-1)); % scale step-length so total traverse = speed
    end
    fprintf('\n=== Time step %d / %d proton = (%.3e, %.3e) ===\n', it, Nt, proton_pos(1), proton_pos(2));
    proton_traj(it,:) = proton_pos;
    proton_vel(it,:) = [vx, vy];
    % ---------------- 1) UNIFORM SOLVE --------------------------------
    Nuni = size(nodes_uni,1);
    [A_u, b_u, ~, ~, ~] = assemble_system(nodes_uni, elem_uni, Nuni, size(elem_uni,1), proton_pos, @(x,y) point_rho(x,y,proton_pos,R0,coef));
    onB = (abs(nodes_uni(:,1)) < 1e-16) | (abs(nodes_uni(:,1)-L) < 1e-16) | ...
          (abs(nodes_uni(:,2)) < 1e-16) | (abs(nodes_uni(:,2)-L) < 1e-16);
    free_u = find(~onB);
    A_ff = A_u(free_u, free_u);
    b_f = b_u(free_u);
    % diagonal regularization if needed
    dmax = max(abs(diag(A_ff)));
    if dmax == 0, dmax = 1; end
    alpha = reg_rel * dmax;
    A_ff_reg = A_ff + alpha * speye(size(A_ff));
    phi_uni = zeros(Nuni,1);
    try
        phi_uni(free_u) = A_ff_reg \ b_f;
    catch
        phi_uni(free_u) = lsqr(A_ff_reg, b_f, 1e-10, 2000);
    end
    % ---------------- 2) AMR solve (reset mesh, a few refinement iterations) ----------
    nodes_amr = nodes_uni;
    elem_amr = elem_uni;
    edge_mid_map = containers.Map('KeyType','char','ValueType','int32'); % reset map each time
    N_amr_initial = size(nodes_amr,1);
    for amr_it = 1:MAX_AMR_ITER
        N_amr_now = size(nodes_amr,1);
        Ne_amr_now = size(elem_amr,1);
        [A_a, b_a, dN_dx_list, dN_dy_list, Area_list] = assemble_system(nodes_amr, elem_amr, N_amr_now, Ne_amr_now, proton_pos, @(x,y) point_rho(x,y,proton_pos,R0,coef));
        onB_a = (abs(nodes_amr(:,1)) < 1e-16) | (abs(nodes_amr(:,1)-L) < 1e-16) | ...
               (abs(nodes_amr(:,2)) < 1e-16) | (abs(nodes_amr(:,2)-L) < 1e-16);
        free_a = find(~onB_a);
        if isempty(free_a)
            warning('AMR: no free nodes -> break');
            break;
        end
        A_ff_a = A_a(free_a, free_a);
        b_f_a = b_a(free_a);
        dmax_a = max(abs(diag(A_ff_a)));
        if dmax_a == 0, dmax_a = 1; end
        alpha_a = reg_rel * dmax_a;
        A_ff_a = A_ff_a + alpha_a * speye(size(A_ff_a));
        phi_amr_full = zeros(N_amr_now,1);
        try
            phi_amr_full(free_a) = A_ff_a \ b_f_a;
        catch
            phi_amr_full(free_a) = lsqr(A_ff_a, b_f_a, 1e-10, 2000);
        end
        % error estimator (normalized)
        eta = error_estimation(nodes_amr, elem_amr, phi_amr_full, dN_dx_list, dN_dy_list, Area_list);
        % if eta empty or below tol, accept
        if isempty(eta) || norm(eta) < tol_eta_stop
            phi_amr = phi_amr_full;
            break;
        end
        marked = mark_elements(eta, THETA);
        if isempty(marked)
            phi_amr = phi_amr_full;
            break;
        end
        marked_closed = closure_neighbours(marked, elem_amr);
        [nodes_amr, elem_amr, edge_mid_map] = refine_mesh(nodes_amr, elem_amr, marked_closed, edge_mid_map);
        % continue with refined mesh
    end
    % ensure phi_amr exists & consistent (if not, recompute once more)
    if ~exist('phi_amr','var') || length(phi_amr) ~= size(nodes_amr,1)
        N_amr_now = size(nodes_amr,1);
        [A_a, b_a, ~, ~, ~] = assemble_system(nodes_amr, elem_amr, N_amr_now, size(elem_amr,1), proton_pos, @(x,y) point_rho(x,y,proton_pos,R0,coef));
        onB_a = (abs(nodes_amr(:,1)) < 1e-16) | (abs(nodes_amr(:,1)-L) < 1e-16) | ...
               (abs(nodes_amr(:,2)) < 1e-16) | (abs(nodes_amr(:,2)-L) < 1e-16);
        free_a = find(~onB_a);
        A_ff_a = A_a(free_a, free_a);
        b_f_a = b_a(free_a);
        dmax_a = max(abs(diag(A_ff_a))); if dmax_a==0, dmax_a=1; end
        A_ff_a = A_ff_a + reg_rel*dmax_a*speye(size(A_ff_a));
        phi_amr = zeros(N_amr_now,1);
        phi_amr(free_a) = A_ff_a \ b_f_a;
    end
    new_nodes = size(nodes_amr,1) - N_amr_initial;
    nodecount_amr_history(it) = size(nodes_amr,1);
    fprintf('New nodes created by AMR: %d\n', new_nodes);
    % ---------------- Analytical potentials at nodes --------------------
    phi_ana_nodes_uni = phi_ana(nodes_uni(:,1), nodes_uni(:,2), proton_pos(1), proton_pos(2));
    phi_ana_nodes_amr = phi_ana(nodes_amr(:,1), nodes_amr(:,2), proton_pos(1), proton_pos(2));
    % ---------------- Error metrics (exclude near singular) -------------
    valid_region = @(nodes, r) (nodes(:,1) > 0.05*L) & (nodes(:,1) < 0.95*L) & ...
                               (nodes(:,2) > 0.05*L) & (nodes(:,2) < 0.95*L) & ...
                               (r > 5*R0);
    r_uni = sqrt((nodes_uni(:,1)-proton_pos(1)).^2 + (nodes_uni(:,2)-proton_pos(2)).^2);
    r_amr = sqrt((nodes_amr(:,1)-proton_pos(1)).^2 + (nodes_amr(:,2)-proton_pos(2)).^2);
    mask_u = valid_region(nodes_uni, r_uni);
    mask_a = valid_region(nodes_amr, r_amr);
    if any(mask_u)
        diff_u = phi_uni(mask_u) - phi_ana_nodes_uni(mask_u);
        L2_history_uni(it) = sqrt(mean(diff_u.^2));
        Linf_history_uni(it) = max(abs(diff_u));
    else
        L2_history_uni(it) = NaN; Linf_history_uni(it) = NaN;
    end
    if any(mask_a)
        diff_a = phi_amr(mask_a) - phi_ana_nodes_amr(mask_a);
        L2_history_amr(it) = sqrt(mean(diff_a.^2));
        Linf_history_amr(it) = max(abs(diff_a));
    else
        L2_history_amr(it) = NaN; Linf_history_amr(it) = NaN;
    end
    fprintf('Nodes: uniform=%d, AMR=%d | L2_err_uni=%.3e L2_err_amr=%.3e\n', size(nodes_uni,1), size(nodes_amr,1), L2_history_uni(it), L2_history_amr(it));
    % ---------------- Cross-section line centered at proton_y -----------
    % choose line spanning [proton_x - halfspan, proton_x + halfspan] but inside [0 L]
    halfspan = max(0.15*L, 0.3*(L/Nx)); % either a chunk or grid spacing-based
    x_left = max(0, proton_pos(1)-halfspan);
    x_right = min(L, proton_pos(1)+halfspan);
    line_x = linspace(x_left, x_right, probe_N)';
    line_y = proton_pos(2) * ones(size(line_x));
    % Interpolate along the probe: use triangulation -> barycentric (triinterp)
    TR_uni = triangulation(elem_uni, nodes_uni);
    TR_amr = triangulation(elem_amr, nodes_amr);
    phi_uni_line = triinterp(TR_uni, phi_uni, [line_x, line_y]);
    phi_amr_line = triinterp(TR_amr, phi_amr, [line_x, line_y]);
    phi_ana_line = phi_ana(line_x, line_y, proton_pos(1), proton_pos(2));
    % Fill NaNs by local interpolation so plots don't have holes
    phi_uni_line_plot = fill_line_nans(line_x, phi_uni_line);
    phi_amr_line_plot = fill_line_nans(line_x, phi_amr_line);
    phi_ana_line_plot = phi_ana_line;
    % No clipping for cross-section, but avoid negative for log if needed
    % Pointwise error along line (absolute)
    err_line_uni = abs(phi_uni_line_plot - phi_ana_line_plot);
    err_line_amr = abs(phi_amr_line_plot - phi_ana_line_plot);
    L2_line_uni = sqrt(mean(err_line_uni(~isnan(err_line_uni)).^2));
    Linf_line_uni= max(err_line_uni(~isnan(err_line_uni)));
    L2_line_amr = sqrt(mean(err_line_amr(~isnan(err_line_amr)).^2));
    Linf_line_amr= max(err_line_amr(~isnan(err_line_amr)));
    fprintf(' Cross-section L2 (uni/amr) = %.3e / %.3e\n', L2_line_uni, L2_line_amr);
    % Node markers close to line
    node_band_factor = 0.6;
    dx_uniform = L / Nx;
    tol_line_uni = node_band_factor * dx_uniform;
    idx_nodes_line_uni = find(abs(nodes_uni(:,2)-line_y(1)) <= tol_line_uni);
    x_nodes_line_uni = nodes_uni(idx_nodes_line_uni,1);
    phi_nodes_line_uni = phi_uni(idx_nodes_line_uni);
    ys_amr_unique = unique(round(nodes_amr(:,2), 12));
    if numel(ys_amr_unique) > 1
        dy_amr = median(diff(sort(ys_amr_unique)));
    else
        dy_amr = dx_uniform;
    end
    tol_line_amr = node_band_factor * max(dy_amr, dx_uniform*0.25);
    idx_nodes_line_amr = find(abs(nodes_amr(:,2)-line_y(1)) <= tol_line_amr);
    x_nodes_line_amr = nodes_amr(idx_nodes_line_amr,1);
    phi_nodes_line_amr = phi_amr(idx_nodes_line_amr);
    % ------------------ Update plots (only every plot_every steps) -----
    if mod(it, plot_every) == 0
        % --- Uniform surface (log scale)
        axes(axUniform3D); cla(axUniform3D);
        phi_nodes_uni_vis = max(phi_uni, 10^min_log_phi);
        try
            hU = trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), zeros(size(nodes_uni,1),1), phi_nodes_uni_vis, 'EdgeColor','none', 'Parent', axUniform3D);
            shading(axUniform3D,'interp');
        catch
            % fallback: patch
            trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), zeros(size(nodes_uni,1),1), phi_nodes_uni_vis,'EdgeColor','none'); shading interp;
        end
        colormap(axUniform3D, jet);
        set(axUniform3D, 'ColorScale', 'log');
        title(axUniform3D, sprintf('Uniform numeric (N=%d)', size(nodes_uni,1)));
        view(axUniform3D,0,90); axis(axUniform3D,'equal'); xlim(axUniform3D,[0 L]); ylim(axUniform3D,[0 L]);
        caxis(axUniform3D,[10^min_log_phi max_phi]);
        colorbar(axUniform3D);
        % --- AMR surface
        axes(axAMR3D); cla(axAMR3D);
        phi_nodes_amr_vis = max(phi_amr, 10^min_log_phi);
        trisurf(elem_amr, nodes_amr(:,1), nodes_amr(:,2), zeros(size(nodes_amr,1),1), phi_nodes_amr_vis, 'EdgeColor','none', 'Parent', axAMR3D);
        shading(axAMR3D,'interp'); colormap(axAMR3D, jet);
        set(axAMR3D, 'ColorScale', 'log');
        title(axAMR3D, sprintf('AMR numeric (N=%d)', size(nodes_amr,1)));
        view(axAMR3D,0,90); axis(axAMR3D,'equal'); xlim(axAMR3D,[0 L]); ylim(axAMR3D,[0 L]);
        caxis(axAMR3D, [10^min_log_phi max_phi]);
        colorbar(axAMR3D);
        % --- Analytical surface (sample coarse grid for speed)
        axes(axAnalytic3D); cla(axAnalytic3D);
        Ng = 120;
        gx = linspace(0,L,Ng); gy = linspace(0,L,Ng);
        [GX,GY] = meshgrid(gx,gy);
        GA = phi_ana(GX, GY, proton_pos(1), proton_pos(2));
        GA_vis = max(GA, 10^min_log_phi);
        surf(axAnalytic3D, GX, GY, zeros(size(GX)), GA_vis, 'EdgeColor','none'); shading(axAnalytic3D,'interp');
        colormap(axAnalytic3D, jet);
        set(axAnalytic3D, 'ColorScale', 'log');
        title(axAnalytic3D,'Analytical (regularized)');
        view(axAnalytic3D,0,90); axis(axAnalytic3D,'equal'); xlim(axAnalytic3D,[0 L]); ylim(axAnalytic3D,[0 L]);
        caxis(axAnalytic3D, [10^min_log_phi max_phi]);
        colorbar(axAnalytic3D);
        % --- Cross-section potentials & node markers
        axes(axCrossPot); cla(axCrossPot);
        plot(axCrossPot, line_x, phi_ana_line_plot, 'r:', 'LineWidth', 2); hold(axCrossPot,'on');
        plot(axCrossPot, line_x, phi_uni_line_plot, 'b-', 'LineWidth', 1.4);
        plot(axCrossPot, line_x, phi_amr_line_plot, 'g--', 'LineWidth', 1.4);
        % node markers
        if ~isempty(idx_nodes_line_uni)
            scatter(axCrossPot, x_nodes_line_uni, phi_nodes_line_uni, 18, 'o', 'MarkerEdgeColor','k','MarkerFaceColor','none');
        end
        if ~isempty(idx_nodes_line_amr)
            scatter(axCrossPot, x_nodes_line_amr, phi_nodes_line_amr, 24, 's', 'MarkerEdgeColor','k','MarkerFaceColor','none');
        end
        xlabel(axCrossPot,'x'); ylabel(axCrossPot,'\phi'); xlim(axCrossPot,[min(line_x) max(line_x)]);
        ylim(axCrossPot,[0 max_phi]); legend(axCrossPot, {'Analytical','Uniform','AMR','Uni nodes','AMR nodes'}, 'Location','best'); grid(axCrossPot,'on');
        % --- Cross-section error (log)
        axes(axCrossErr); cla(axCrossErr);
        semilogy(axCrossErr, line_x, max(err_line_uni, eps), 'b-', 'LineWidth', 1.2); hold(axCrossErr,'on');
        semilogy(axCrossErr, line_x, max(err_line_amr, eps), 'g--', 'LineWidth', 1.2);
        xlabel(axCrossErr,'x'); ylabel(axCrossErr,'|phi_{num}-phi_{ana}|'); grid(axCrossErr,'on');
        legend(axCrossErr, {sprintf('Uniform (L2=%.1e)', L2_line_uni), sprintf('AMR (L2=%.1e)', L2_line_amr)}, 'Location','best');
        % --- AMR Mesh plot
        axes(axMesh); cla(axMesh);
        triplot(elem_amr, nodes_amr(:,1), nodes_amr(:,2), 'k-');
        hold on;
        theta = linspace(0,2*pi,100);
        plot(axMesh, proton_pos(1) + 3*R0*cos(theta), proton_pos(2) + 3*R0*sin(theta), 'r-');
        text(0.05*L, 0.95*L, sprintf('#nodes: %d', size(nodes_amr,1)));
        text(0.05*L, 0.85*L, sprintf('#elements: %d', size(elem_amr,1)));
        text(0.05*L, 0.75*L, sprintf('new nodes: %d', new_nodes));
        axis(axMesh,'equal'); xlim(axMesh,[0 L]); ylim(axMesh,[0 L]);
        drawnow;
        pause(0.02); % small pause to allow GUI update
    end
end
%% ------------------ final summary plots & trajectory -------------------
figure('Name','Proton Trajectory & Velocity','Units','normalized','Position',[0.2 0.15 0.6 0.5]);
subplot(1,2,1);
plot(proton_traj(:,1), proton_traj(:,2), '-o'); axis equal; xlim([0 L]); ylim([0 L]);
xlabel('x'); ylabel('y'); title('Proton trajectory');
subplot(1,2,2);
plot(1:Nt, sqrt(sum(proton_vel.^2,2)), '-k'); xlabel('time step'); ylabel('speed'); title('Proton speed (constant)');
figure('Name','Errors & Node Counts over time','Units','normalized','Position',[0.15 0.05 0.7 0.4]);
yyaxis left;
plot(1:Nt, nodecount_amr_history, '-k', 'LineWidth', 1.6); hold on;
plot(1:Nt, repmat(N_uni, Nt,1), ':b', 'LineWidth', 1);
ylabel('AMR node count');
yyaxis right;
plot(1:Nt, L2_history_uni, 'b--', 'LineWidth', 1.2); hold on;
plot(1:Nt, L2_history_amr, 'g-', 'LineWidth', 1.6);
ylabel('L2 error (away from singular)');
xlabel('time step');
legend({'AMR N','Uniform N','L2 uni','L2 amr'}, 'Location','best');
grid on;
fprintf('\n--- Simulation finished ---\n');
fprintf('Final node counts: uniform=%d AMR=%d\n', size(nodes_uni,1), size(nodes_amr,1));
fprintf('Final L2 errors (uni/amr) = %.3e / %.3e\n', L2_history_uni(end), L2_history_amr(end));
%% ================================ HELPERS ==============================
% All helper functions are local below
function [nodes, elem] = generate_initial_mesh(L, Nx, Ny)
    x_coords = linspace(0, L, Nx+1);
    y_coords = linspace(0, L, Ny+1);
    [Xg, Yg] = meshgrid(x_coords, y_coords);
    nodes = [Xg(:), Yg(:)];
    elem = zeros(2*Nx*Ny, 3);
    N_x = Nx + 1;
    idx = 1;
    for j = 1:Ny
        for i = 1:Nx
            n1 = (j-1)*N_x + i;
            n2 = n1 + 1;
            n3 = n1 + N_x;
            n4 = n2 + N_x;
            elem(idx,:) = [n1 n2 n3]; idx = idx + 1;
            elem(idx,:) = [n2 n4 n3]; idx = idx + 1;
        end
    end
    elem = elem(1:idx-1,:);
end
function [A, b, dN_dx_list, dN_dy_list, Area_list] = assemble_system(nodes, elem, N, Ne, center, rho_func)
    A = sparse(N, N); b = zeros(N,1);
    dN_dx_list = zeros(3, Ne); dN_dy_list = zeros(3, Ne); Area_list = zeros(Ne,1);
    for e = 1:Ne
        ids = elem(e,:);
        x = nodes(ids,1); y = nodes(ids,2);
        Area = 0.5 * abs( x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)) );
        Area_list(e) = Area;
        if Area < eps, continue; end
        dN_dx = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / (2*Area);
        dN_dy = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / (2*Area);
        dN_dx_list(:, e) = dN_dx; dN_dy_list(:, e) = dN_dy;
        Ke = Area * (dN_dx * dN_dx' + dN_dy * dN_dy');
        A(ids, ids) = A(ids, ids) + Ke;
        % load vector using rho at element centroid
        xc = mean(x); yc = mean(y);
        rho = rho_func(xc, yc);
        b(ids) = b(ids) + (Area * rho / 3) * ones(3,1);
    end
end
function rho = point_rho(x,y,center,R0,coef)
    % Correct rho = coef * (2 R0^2 - r^2) / (r^2 + R0^2)^{5/2} (positive at center)
    r2 = (x-center(1)).^2 + (y-center(2)).^2;
    s2 = r2 + R0^2;
    rho = coef * (2 * R0^2 - r2) ./ (s2 .^(5/2));
end
function eta = error_estimation(nodes, elem, phi, dN_dx_list, dN_dy_list, Area_list)
    % Edge-jump estimator (residual-like). Returns per-element indicator.
    Ne = size(elem,1);
    phi_x = zeros(Ne,1); phi_y = zeros(Ne,1);
    for e = 1:Ne
        ids = elem(e,:);
        phi_e = phi(ids);
        phi_x(e) = dN_dx_list(:,e)' * phi_e;
        phi_y(e) = dN_dy_list(:,e)' * phi_e;
    end
    TR = triangulation(elem, nodes);
    E_adj = edges(TR);
    T_adj = edgeAttachments(TR, E_adj);
    eta_sq = zeros(Ne,1);
    for k = 1:size(E_adj,1)
        T_ids = T_adj{k};
        if numel(T_ids) ~= 2, continue; end
        e1 = T_ids(1); e2 = T_ids(2);
        grad_jump_sq = (phi_x(e1) - phi_x(e2))^2 + (phi_y(e1) - phi_y(e2))^2;
        jump_mag = sqrt(grad_jump_sq);
        n_ids = E_adj(k,:);
        h_E = norm(nodes(n_ids(1),:) - nodes(n_ids(2),:));
        indicator = jump_mag * h_E;
        eta_sq(e1) = eta_sq(e1) + indicator^2;
        eta_sq(e2) = eta_sq(e2) + indicator^2;
    end
    eta = sqrt(eta_sq);
    if max(eta) > 0
        eta = eta / max(eta);
    end
end
function marked_elements = mark_elements(eta, theta)
    % Dörfler marking: select minimal set whose squared-sum >= theta * total
    if isempty(eta)
        marked_elements = [];
        return;
    end
    eta_sq = eta.^2;
    total = sum(eta_sq);
    if total == 0
        marked_elements = [];
        return;
    end
    [~, idx] = sort(eta_sq, 'descend');
    csum = 0; k = 0;
    while csum < theta*total && k < numel(idx)
        k = k + 1;
        csum = csum + eta_sq(idx(k));
    end
    marked_elements = idx(1:k);
end
function marked_closed = closure_neighbours(marked_elements, elem)
    % make a simple 1-ring closure around marked elements for conformity
    Ne = size(elem,1);
    neighbors = build_neighbors(elem);
    marked = false(Ne,1);
    if ~isempty(marked_elements)
        marked(marked_elements) = true;
    end
    cur = find(marked);
    for e = cur'
        marked(neighbors{e}) = true;
    end
    marked_closed = find(marked);
end
function neighbors = build_neighbors(elem)
    Ne = size(elem,1);
    edge_map = containers.Map('KeyType','char','ValueType','any');
    for e = 1:Ne
        tri = elem(e,:);
        edges = [tri([1 2]); tri([2 3]); tri([3 1])];
        for k = 1:3
            ed = sort(edges(k,:));
            key = sprintf('%d-%d', ed(1), ed(2));
            if isKey(edge_map, key)
                edge_map(key) = [edge_map(key), e];
            else
                edge_map(key) = e;
            end
        end
    end
    neighbors = cell(Ne,1);
    vals = values(edge_map);
    for i = 1:length(vals)
        elems_sharing_edge = vals{i};
        if numel(elems_sharing_edge) < 2, continue; end
        for j = 1:numel(elems_sharing_edge)
            e1 = elems_sharing_edge(j);
            others = setdiff(elems_sharing_edge, e1);
            neighbors{e1} = unique([neighbors{e1}, others]);
        end
    end
end
function [new_nodes, new_elem, edge_mid_map] = refine_mesh(nodes, elem, marked_elements, edge_mid_map)
    % Midpoint refinement of marked triangles; reuses midpoints in edge_mid_map
    if nargin < 4 || isempty(edge_mid_map)
        edge_mid_map = containers.Map('KeyType','char','ValueType','int32');
    end
    N_old = size(nodes,1);
    Ne_old = size(elem,1);
    new_nodes = nodes;
    next_node_id = N_old + 1;
    new_elem = zeros(4*length(marked_elements) + (Ne_old-length(marked_elements)), 3);
    idx = 1;
    marked_mask = false(Ne_old,1);
    if ~isempty(marked_elements)
        marked_mask(marked_elements) = true;
    end
    for e = 1:Ne_old
        tri = elem(e,:);
        n1 = tri(1); n2 = tri(2); n3 = tri(3);
        if ~marked_mask(e)
            new_elem(idx,:) = tri; idx = idx + 1;
            continue;
        end
        edges = [n1 n2; n2 n3; n3 n1];
        mids = zeros(3,1);
        for k = 1:3
            a = edges(k,1); b = edges(k,2);
            key = sprintf('%d-%d', min(a,b), max(a,b));
            if isKey(edge_mid_map, key)
                mids(k) = edge_mid_map(key);
            else
                new_nodes(next_node_id,:) = 0.5*(nodes(a,:) + nodes(b,:));
                edge_mid_map(key) = next_node_id;
                mids(k) = next_node_id;
                next_node_id = next_node_id + 1;
            end
        end
        m12 = mids(1); m23 = mids(2); m31 = mids(3);
        new_elem(idx,:) = [n1, m12, m31]; idx = idx + 1;
        new_elem(idx,:) = [n2, m23, m12]; idx = idx + 1;
        new_elem(idx,:) = [n3, m31, m23]; idx = idx + 1;
        new_elem(idx,:) = [m12, m23, m31]; idx = idx + 1;
    end
    new_elem = new_elem(1:idx-1,:);
    maxIdx = max(new_elem(:));
    if maxIdx > size(new_nodes,1)
        error('REFINE_MESH: invalid mesh: element references node %d but only %d nodes exist', maxIdx, size(new_nodes,1));
    end
end
function Vq = triinterp(TR, V, Pq)
    % TRIINTERP - barycentric interpolation on TR at Pq (M x 2)
    nNodesTR = size(TR.Points,1);
    if length(V) ~= nNodesTR
        error('TRIINTERP: length(V) incompatible with TR.');
    end
    elem_idx = pointLocation(TR, Pq);
    Vq = nan(size(Pq,1),1);
    for i = 1:size(Pq,1)
        e = elem_idx(i);
        if isnan(e), continue; end
        ids = TR.ConnectivityList(e,:);
        nodes_e = TR.Points(ids,:);
        x = nodes_e(:,1); y = nodes_e(:,2);
        xq = Pq(i,1); yq = Pq(i,2);
        Area_e = 0.5 * abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)));
        if Area_e == 0
            Vq(i) = NaN; continue;
        end
        L1 = (xq*(y(2)-y(3)) + x(2)*(y(3)-yq) + x(3)*(yq-y(2))) / (2*Area_e);
        L2 = (xq*(y(3)-y(1)) + x(3)*(y(1)-yq) + x(1)*(yq-y(3))) / (2*Area_e);
        L3 = 1 - L1 - L2;
        Vq(i) = L1*V(ids(1)) + L2*V(ids(2)) + L3*V(ids(3));
    end
end
function v_filled = fill_line_nans(xq, v)
    v_filled = v;
    nans = isnan(v);
    if any(nans)
        xi = xq(~nans);
        vi = v(~nans);
        if ~isempty(xi)
            v_interp = interp1(xi, vi, xq(nans), 'pchip', 'extrap');
            v_filled(nans) = v_interp;
        end
    end
end