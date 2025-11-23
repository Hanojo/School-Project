% run_static_proton_AMR_relativistic_fixed.m
clear; close all; clc;

%% ------------------------- physical / domain params --------------------
L = 1; 
R0 = 0.02; 
coef = 1; 
epsilon0 = 1; % For Poisson equation: ∇²φ = -ρ/ε₀

%% ------------------------- simulation params ---------------------------
Nx = 60; Ny = 60;
MAX_AMR_ITER = 6; 
THETA = 0.25; 
tol_eta_stop = 0.05; 
Nt = 1; 
beta = 0.9; % Increased for clearer relativistic effects
gamma = 1 / sqrt(1 - beta^2);

% Motion direction (x-direction for relativistic contraction)
v_dir = [1, 0]; % velocity direction (unit vector)

% CORRECTED: Relativistic analytical potential
% For a moving charge, the potential is contracted in the direction of motion
phi_ana = @(x,y,xc,yc) coef * gamma ./ sqrt(gamma^2*(x-xc).^2 + (y-yc).^2 + R0^2);

% CORRECTED: Relativistic charge density
% The charge density transforms as ρ = γρ₀ in the lab frame
rho_ana = @(x,y,xc,yc) coef * gamma * (3*R0^2) ./ (gamma^2*(x-xc).^2 + (y-yc).^2 + R0^2).^(5/2);

proton_pos = [0.5*L, 0.5*L];

fprintf('Relativistic parameters: beta=%.3f, gamma=%.3f\n', beta, gamma);
fprintf('Expected contraction factor in x-direction: %.3f\n', 1/gamma);

%% ----------------------- meshes initialization -------------------------
[nodes_uni, elem_uni] = generate_initial_mesh(L, Nx, Ny);

%% ----------------------- main computation ------------------------------
fprintf('\n=== Starting computation ===\n');

% 1) UNIFORM SOLVE - WITH RELATIVISTIC CORRECTIONS
fprintf('Solving uniform mesh (with relativistic corrections)...\n');
[A_u, b_u] = assemble_system_relativistic(nodes_uni, elem_uni, proton_pos, rho_ana, gamma, v_dir);

% Boundary conditions
tol_bc = 1e-8;
onB = (nodes_uni(:,1) <= tol_bc) | (nodes_uni(:,1) >= L-tol_bc) | ...
      (nodes_uni(:,2) <= tol_bc) | (nodes_uni(:,2) >= L-tol_bc);
onB_idx = find(onB);
free_u = find(~onB);

phi_ana_uni = phi_ana(nodes_uni(:,1), nodes_uni(:,2), proton_pos(1), proton_pos(2));

% Apply boundary conditions
A_ff = A_u(free_u, free_u);
b_f = b_u(free_u) - A_u(free_u, onB_idx) * phi_ana_uni(onB_idx);

% Solve
phi_uni = zeros(size(nodes_uni,1),1);
phi_uni(free_u) = A_ff \ b_f;
phi_uni(onB) = phi_ana_uni(onB);

fprintf('Uniform: max phi = %.3e, min phi = %.3e\n', max(phi_uni), min(phi_uni));

% 2) AMR SOLVE - WITH RELATIVISTIC CORRECTIONS
fprintf('Starting AMR iterations (with relativistic corrections)...\n');
nodes_amr = nodes_uni;
elem_amr = elem_uni;
phi_amr_final = [];

for amr_it = 1:MAX_AMR_ITER
    fprintf('AMR iteration %d: nodes=%d, elements=%d\n', amr_it, size(nodes_amr,1), size(elem_amr,1));
    
    [A_a, b_a, dN_dx_list, dN_dy_list, Area_list] = assemble_system_relativistic(nodes_amr, elem_amr, proton_pos, rho_ana, gamma, v_dir);
    
    phi_ana_amr = phi_ana(nodes_amr(:,1), nodes_amr(:,2), proton_pos(1), proton_pos(2));
    
    onB_a = (nodes_amr(:,1) <= tol_bc) | (nodes_amr(:,1) >= L-tol_bc) | ...
           (nodes_amr(:,2) <= tol_bc) | (nodes_amr(:,2) >= L-tol_bc);
    onB_a_idx = find(onB_a);
    free_a = find(~onB_a);
    
    if isempty(free_a)
        warning('No free nodes in AMR mesh');
        phi_amr_final = phi_ana_amr;
        break;
    end
    
    A_ff_a = A_a(free_a, free_a);
    b_f_a = b_a(free_a) - A_a(free_a, onB_a_idx) * phi_ana_amr(onB_a_idx);
    
    phi_amr_current = zeros(size(nodes_amr,1),1);
    phi_amr_current(free_a) = A_ff_a \ b_f_a;
    phi_amr_current(onB_a) = phi_ana_amr(onB_a);
    
    phi_amr_final = phi_amr_current;
    
    % Error estimation
    eta = error_estimation_relativistic(nodes_amr, elem_amr, phi_amr_current, dN_dx_list, dN_dy_list, Area_list, gamma, v_dir);
    
    if isempty(eta) || max(eta) < tol_eta_stop
        fprintf('AMR converged after %d iterations\n', amr_it);
        break;
    end
    
    marked = mark_elements(eta, THETA);
    if isempty(marked)
        fprintf('No elements marked for refinement\n');
        break;
    end
    
    marked_closed = closure_neighbours(marked, elem_amr);
    [nodes_amr, elem_amr] = refine_mesh_fixed(nodes_amr, elem_amr, marked_closed);
    
    if amr_it == MAX_AMR_ITER
        fprintf('Reached maximum AMR iterations\n');
        % Final solve on refined mesh
        [A_a, b_a, ~, ~, ~] = assemble_system_relativistic(nodes_amr, elem_amr, proton_pos, rho_ana, gamma, v_dir);
        phi_ana_amr_final = phi_ana(nodes_amr(:,1), nodes_amr(:,2), proton_pos(1), proton_pos(2));
        onB_a_final = (nodes_amr(:,1) <= tol_bc) | (nodes_amr(:,1) >= L-tol_bc) | ...
                     (nodes_amr(:,2) <= tol_bc) | (nodes_amr(:,2) >= L-tol_bc);
        onB_a_idx_final = find(onB_a_final);
        free_a_final = find(~onB_a_final);
        
        A_ff_a_final = A_a(free_a_final, free_a_final);
        b_f_a_final = b_a(free_a_final) - A_a(free_a_final, onB_a_idx_final) * phi_ana_amr_final(onB_a_idx_final);
        
        phi_amr_final = zeros(size(nodes_amr,1),1);
        phi_amr_final(free_a_final) = A_ff_a_final \ b_f_a_final;
        phi_amr_final(onB_a_final) = phi_ana_amr_final(onB_a_final);
    end
end

% Ensure final solution
if isempty(phi_amr_final) || length(phi_amr_final) ~= size(nodes_amr,1)
    fprintf('Final AMR solution invalid, recomputing...\n');
    [A_a, b_a, ~, ~, ~] = assemble_system_relativistic(nodes_amr, elem_amr, proton_pos, rho_ana, gamma, v_dir);
    phi_ana_amr_final = phi_ana(nodes_amr(:,1), nodes_amr(:,2), proton_pos(1), proton_pos(2));
    onB_a_final = (nodes_amr(:,1) <= tol_bc) | (nodes_amr(:,1) >= L-tol_bc) | ...
                 (nodes_amr(:,2) <= tol_bc) | (nodes_amr(:,2) >= L-tol_bc);
    onB_a_idx_final = find(onB_a_final);
    free_a_final = find(~onB_a_final);
    
    A_ff_a_final = A_a(free_a_final, free_a_final);
    b_f_a_final = b_a(free_a_final) - A_a(free_a_final, onB_a_idx_final) * phi_ana_amr_final(onB_a_idx_final);
    
    phi_amr_final = zeros(size(nodes_amr,1),1);
    phi_amr_final(free_a_final) = A_ff_a_final \ b_f_a_final;
    phi_amr_final(onB_a_final) = phi_ana_amr_final(onB_a_final);
end

fprintf('AMR final: nodes=%d, elements=%d\n', size(nodes_amr,1), size(elem_amr,1));

%% ----------------------- analysis and visualization --------------------
phi_amr = phi_amr_final;

% Create comparison plots
figure('Position',[100 100 1400 800], 'Name', 'Relativistic Effects Comparison');

% 1) Potential surfaces comparison
subplot(2,3,1);
if size(nodes_uni,1) == length(phi_uni)
    trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), zeros(size(nodes_uni,1),1), phi_uni, 'EdgeColor','none');
else
    triplot(elem_uni, nodes_uni(:,1), nodes_uni(:,2), 'k-');
end
title('Uniform Numerical');
view(0,90); axis equal; colorbar;
xlabel('x'); ylabel('y');

subplot(2,3,2);
if size(nodes_amr,1) == length(phi_amr)
    trisurf(elem_amr, nodes_amr(:,1), nodes_amr(:,2), zeros(size(nodes_amr,1),1), phi_amr, 'EdgeColor','none');
else
    triplot(elem_amr, nodes_amr(:,1), nodes_amr(:,2), 'k-');
end
title('AMR Numerical');
view(0,90); axis equal; colorbar;
xlabel('x'); ylabel('y');

subplot(2,3,3);
% Analytical solution with proper sampling to show relativistic effects
Ng = 150;
gx = linspace(0,L,Ng); gy = linspace(0,L,Ng);
[GX,GY] = meshgrid(gx,gy);
GA = phi_ana(GX, GY, proton_pos(1), proton_pos(2));
surf(GX, GY, zeros(size(GX)), GA, 'EdgeColor','none');
title('Analytical Relativistic');
view(0,90); axis equal; colorbar;
xlabel('x'); ylabel('y');

% 2) Cross-sections to show relativistic distortion
subplot(2,3,4);
x_line = linspace(0.1*L, 0.9*L, 300);
y_line = proton_pos(2) * ones(size(x_line));

% Interpolate numerical solutions
valid_uni = ~isnan(phi_uni) & isfinite(phi_uni);
valid_amr = ~isnan(phi_amr) & isfinite(phi_amr);

if any(valid_uni)
    phi_line_uni = griddata(nodes_uni(valid_uni,1), nodes_uni(valid_uni,2), phi_uni(valid_uni), x_line, y_line, 'linear');
else
    phi_line_uni = nan(size(x_line));
end

if any(valid_amr)
    phi_line_amr = griddata(nodes_amr(valid_amr,1), nodes_amr(valid_amr,2), phi_amr(valid_amr), x_line, y_line, 'linear');
else
    phi_line_amr = nan(size(x_line));
end

phi_line_ana = phi_ana(x_line, y_line, proton_pos(1), proton_pos(2));

plot(x_line, phi_line_ana, 'k-', 'LineWidth', 3, 'DisplayName','Analytical');
hold on;
plot(x_line, phi_line_uni, 'b-', 'LineWidth', 2, 'DisplayName','Uniform');
plot(x_line, phi_line_amr, 'r--', 'LineWidth', 2, 'DisplayName','AMR');
xline(proton_pos(1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName','Proton Center');
xlabel('x'); ylabel('\phi');
title('Horizontal cut (y=0.5) - SHOULD SHOW CONTRACTION');
legend('Location','best'); grid on;

% 3) Vertical cross-section for comparison
subplot(2,3,5);
y_line_vert = linspace(0.1*L, 0.9*L, 300);
x_line_vert = proton_pos(1) * ones(size(y_line_vert));

phi_line_ana_vert = phi_ana(x_line_vert, y_line_vert, proton_pos(1), proton_pos(2));

if any(valid_uni)
    phi_line_uni_vert = griddata(nodes_uni(valid_uni,1), nodes_uni(valid_uni,2), phi_uni(valid_uni), x_line_vert, y_line_vert, 'linear');
else
    phi_line_uni_vert = nan(size(y_line_vert));
end

if any(valid_amr)
    phi_line_amr_vert = griddata(nodes_amr(valid_amr,1), nodes_amr(valid_amr,2), phi_amr(valid_amr), x_line_vert, y_line_vert, 'linear');
else
    phi_line_amr_vert = nan(size(y_line_vert));
end

plot(y_line_vert, phi_line_ana_vert, 'k-', 'LineWidth', 3, 'DisplayName','Analytical');
hold on;
plot(y_line_vert, phi_line_uni_vert, 'b-', 'LineWidth', 2, 'DisplayName','Uniform');
plot(y_line_vert, phi_line_amr_vert, 'r--', 'LineWidth', 2, 'DisplayName','AMR');
yline(proton_pos(2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName','Proton Center');
xlabel('y'); ylabel('\phi');
title('Vertical cut (x=0.5)');
legend('Location','best'); grid on;

% 4) Ratio plot to show relativistic contraction
subplot(2,3,6);
% Calculate the full width at half maximum (FWHM) to quantify contraction
[max_ana, max_idx] = max(phi_line_ana);
half_max = max_ana / 2;

% Find FWHM for horizontal cut
left_idx = find(phi_line_ana(1:max_idx) <= half_max, 1, 'last');
right_idx = find(phi_line_ana(max_idx:end) <= half_max, 1, 'first') + max_idx - 1;

if ~isempty(left_idx) && ~isempty(right_idx)
    fwhm_x = x_line(right_idx) - x_line(left_idx);
    text(0.1, 0.9, sprintf('FWHM_x = %.4f', fwhm_x), 'Units','normalized', 'FontSize', 12);
end

% Find FWHM for vertical cut  
[max_ana_vert, max_idx_vert] = max(phi_line_ana_vert);
half_max_vert = max_ana_vert / 2;

left_idx_vert = find(phi_line_ana_vert(1:max_idx_vert) <= half_max_vert, 1, 'last');
right_idx_vert = find(phi_line_ana_vert(max_idx_vert:end) <= half_max_vert, 1, 'first') + max_idx_vert - 1;

if ~isempty(left_idx_vert) && ~isempty(right_idx_vert)
    fwhm_y = y_line_vert(right_idx_vert) - y_line_vert(left_idx_vert);
    text(0.1, 0.8, sprintf('FWHM_y = %.4f', fwhm_y), 'Units','normalized', 'FontSize', 12);
    
    % Calculate contraction ratio
    if fwhm_x > 0
        contraction_ratio = fwhm_x / fwhm_y;
        expected_ratio = 1/gamma;
        text(0.1, 0.7, sprintf('Ratio (x/y) = %.3f', contraction_ratio), 'Units','normalized', 'FontSize', 12);
        text(0.1, 0.6, sprintf('Expected = %.3f', expected_ratio), 'Units','normalized', 'FontSize', 12);
    end
end

axis off;
title('Relativistic Contraction Analysis');

%% ====================== RELATIVISTIC FINITE ELEMENT FUNCTIONS ==========

function [A, b, dN_dx_list, dN_dy_list, Area_list] = assemble_system_relativistic(nodes, elem, center, rho_func, gamma, v_dir)
    % RELATIVISTIC finite element assembly
    % The relativistic Poisson equation becomes anisotropic
    % ∇·(ε∇φ) = -ρ with anisotropic ε due to relativistic effects
    
    N = size(nodes,1);
    Ne = size(elem,1);
    A = sparse(N, N);
    b = zeros(N,1);
    
    dN_dx_list = zeros(3, Ne);
    dN_dy_list = zeros(3, Ne);
    Area_list = zeros(Ne,1);
    
    % Relativistic permittivity tensor (simplified)
    % For motion along x: ε_xx = 1/γ², ε_yy = 1, ε_xy = ε_yx = 0
    epsilon = [1/gamma^2, 0; 0, 1];
    
    for e = 1:Ne
        ids = elem(e,:);
        x = nodes(ids,1); y = nodes(ids,2);
        
        Area = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
        Area_list(e) = Area;
        
        if Area < 1e-14
            continue;
        end
        
        % Shape function gradients
        dN_dx = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / (2*Area);
        dN_dy = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / (2*Area);
        
        dN_dx_list(:, e) = dN_dx;
        dN_dy_list(:, e) = dN_dy;
        
        % RELATIVISTIC stiffness matrix: K_e = ∫(∇N·ε∇N) dΩ
        % For anisotropic ε, this becomes:
        % K_e = Area * (dN_dx * epsilon(1,1) * dN_dx' + dN_dx * epsilon(1,2) * dN_dy' + ...
        %               dN_dy * epsilon(2,1) * dN_dx' + dN_dy * epsilon(2,2) * dN_dy')
        
        Ke = Area * (epsilon(1,1) * (dN_dx * dN_dx') + ...
                    epsilon(1,2) * (dN_dx * dN_dy') + ...
                    epsilon(2,1) * (dN_dy * dN_dx') + ...
                    epsilon(2,2) * (dN_dy * dN_dy'));
        
        A(ids, ids) = A(ids, ids) + Ke;
        
        % Load vector - relativistic charge density
        xc = mean(x); yc = mean(y);
        rho_val = rho_func(xc, yc, center(1), center(2));
        be = (Area * rho_val / 3) * ones(3,1);
        b(ids) = b(ids) + be;
    end
end

function eta = error_estimation_relativistic(nodes, elem, phi, dN_dx_list, dN_dy_list, Area_list, gamma, v_dir)
    % Relativistic error estimation - accounts for anisotropic gradients
    
    Ne = size(elem,1);
    eta = zeros(Ne,1);
    
    % Relativistic metric for error estimation
    epsilon = [1/gamma^2, 0; 0, 1];
    
    for e = 1:Ne
        ids = elem(e,:);
        phi_e = phi(ids);
        
        % Compute gradient in element
        grad_x = dN_dx_list(:,e)' * phi_e;
        grad_y = dN_dy_list(:,e)' * phi_e;
        
        % Relativistic gradient norm
        grad_norm_sq = epsilon(1,1)*grad_x^2 + 2*epsilon(1,2)*grad_x*grad_y + epsilon(2,2)*grad_y^2;
        grad_norm = sqrt(grad_norm_sq);
        
        eta(e) = Area_list(e) * grad_norm;
    end
    
    if max(eta) > 0
        eta = eta / max(eta);
    end
end

%% ====================== STANDARD HELPER FUNCTIONS ======================

function [nodes, elem] = generate_initial_mesh(L, Nx, Ny)
    x_coords = linspace(0, L, Nx+1);
    y_coords = linspace(0, L, Ny+1);
    [Xg, Yg] = meshgrid(x_coords, y_coords);
    nodes = [Xg(:), Yg(:)];
    
    elem = [];
    for j = 1:Ny
        for i = 1:Nx
            n1 = (j-1)*(Nx+1) + i;
            n2 = n1 + 1;
            n3 = n1 + (Nx+1);
            n4 = n2 + (Nx+1);
            
            elem = [elem; n1 n2 n3];
            elem = [elem; n2 n4 n3];
        end
    end
end

function [new_nodes, new_elem] = refine_mesh_fixed(nodes, elem, marked_elements)
    new_nodes = nodes;
    next_node_id = size(nodes,1) + 1;
    
    marked_mask = false(size(elem,1),1);
    marked_mask(marked_elements) = true;
    
    new_elem = [];
    edge_mid_map = containers.Map();
    
    for e = 1:size(elem,1)
        if ~marked_mask(e)
            new_elem = [new_elem; elem(e,:)];
            continue;
        end
        
        n1 = elem(e,1); n2 = elem(e,2); n3 = elem(e,3);
        
        edges = [n1 n2; n2 n3; n3 n1];
        mids = zeros(1,3);
        
        for k = 1:3
            a = min(edges(k,:)); b = max(edges(k,:));
            key = sprintf('%d-%d', a, b);
            
            if isKey(edge_mid_map, key)
                mids(k) = edge_mid_map(key);
            else
                new_nodes(next_node_id,:) = 0.5 * (nodes(a,:) + nodes(b,:));
                edge_mid_map(key) = next_node_id;
                mids(k) = next_node_id;
                next_node_id = next_node_id + 1;
            end
        end
        
        m12 = mids(1); m23 = mids(2); m31 = mids(3);
        
        new_elem = [new_elem;
                   n1, m12, m31;
                   n2, m23, m12;
                   n3, m31, m23;
                   m12, m23, m31];
    end
end

function marked_elements = mark_elements(eta, theta)
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
    Ne = size(elem,1);
    marked = false(Ne,1);
    if ~isempty(marked_elements)
        marked(marked_elements) = true;
    end
    
    neighbors = build_neighbors(elem);
    for e = marked_elements'
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
    keys = edge_map.keys();
    for i = 1:length(keys)
        elems = edge_map(keys{i});
        if length(elems) == 2
            neighbors{elems(1)} = unique([neighbors{elems(1)}, elems(2)]);
            neighbors{elems(2)} = unique([neighbors{elems(2)}, elems(1)]);
        end
    end
end