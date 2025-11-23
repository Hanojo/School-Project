    % ========================================================================
    % PROJECT 3 (Corrected): FEM for Stationary Proton Potential - Uniform Mesh
    % ========================================================================
    % This script solves:  -Δφ = ρ  on a square domain with φ = 0 on the boundary.
    % The charge density ρ represents a stationary proton modeled as a Gaussian.
    % The code:
    %   1. Builds a uniform triangular mesh
    %   2. Assembles the FEM stiffness matrix and load vector
    %   3. Solves the system for φ
    %   4. Compares to analytical φ = 1 / r
    %   5. Visualizes FEM vs analytical and computes error metrics
    %
    % AUTHOR: Johan Lorenzi (corrected version)
    % COMMENTS: Each section is heavily documented for clarity.
    %
    % ========================================================================
    
    clear; close all; clc;
    
    %% =====================================================================
    % 1. PARAMETERS (Non-dimensionalized for numerical stability)
    % =====================================================================
    L       = 1.0;              % Domain size (non-dimensional)
    R0      = 0.01;             % Proton radius (used to regularize 1/r)
    center  = [0.5*L, 0.5*L];  % Proton position
    epsilon = 3 * R0;           % Gaussian width to approximate delta function
    
    % Mesh resolution
    Nx = 50;                    % Number of intervals in x (→ Nx+1 nodes)
    Ny = 50;                    % Number of intervals in y
    
    %% =====================================================================
    % 2. MESH GENERATION: Structured grid → Triangular mesh
    % =====================================================================
    % Create equally spaced node coordinates
    x_coords = linspace(0, L, Nx+1);
    y_coords = linspace(0, L, Ny+1);
    
    % Generate full 2D grid (Xg(i,j), Yg(i,j))
    [Xg, Yg] = meshgrid(x_coords, y_coords);
    
    % Store all node coordinates as a list of (x,y) pairs
    nodes = [Xg(:), Yg(:)];
    N = size(nodes, 1);         % Total number of nodes
    
    % Build triangular elements (each square → 2 triangles)
    elem = [];
    for j = 1:Ny
        for i = 1:Nx
            % Identify the 4 nodes of the current square cell
            n1 = (j-1)*(Nx+1) + i;      % bottom-left
            n2 = n1 + 1;                % bottom-right
            n3 = n1 + (Nx+1);           % top-left
            n4 = n2 + (Nx+1);           % top-right
            
            % Split square into two triangles
            elem = [elem;
                    n1 n2 n3;            % lower-left triangle
                    n2 n4 n3];           % upper-right triangle
        end
    end
    Ne = size(elem, 1);          % Total number of triangular elements
    
    fprintf('Mesh generated: %d nodes, %d triangular elements\n', N, Ne);
    
    %% =====================================================================
    % 3. ASSEMBLY OF GLOBAL STIFFNESS MATRIX (A) AND LOAD VECTOR (b)
    %    Equation:  A φ = b  with A_ij = ∫ ∇Ni·∇Nj dΩ  and  b_i = ∫ Ni ρ dΩ
    % =====================================================================
    A = sparse(N, N);            % Sparse stiffness matrix
    b = zeros(N, 1);             % Load vector
    
    % Define Gaussian charge distribution parameters
    sigma = epsilon / 3;         % std. deviation of Gaussian
    norm_factor = 1 / (2*pi*sigma^2);   % Ensures ∫ρ dxdy = 1
    
    % Loop over all elements to assemble local matrices
    for e = 1:Ne
        % Node indices for this triangle
        ids = elem(e, :);
        
        % Coordinates of the triangle's 3 vertices
        x = nodes(ids, 1);
        y = nodes(ids, 2);
        
        % Compute triangle area (shoelace formula)
        Area = 0.5 * abs( x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)) );
        
        % Skip degenerate elements (none expected in structured mesh)
        if Area < 1e-12
            continue;
        end
        
        % Compute gradients of shape functions (P1 linear elements)
        dN_dx = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / (2*Area);
        dN_dy = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / (2*Area);
        
        % Local stiffness matrix (3x3)
        Ke = Area * (dN_dx * dN_dx' + dN_dy * dN_dy');
        
        % Assemble into global matrix
        A(ids, ids) = A(ids, ids) + Ke;
        
        % ---- LOAD VECTOR: Gaussian charge ----
        % Evaluate ρ at centroid (1-point quadrature)
        xc = mean(x);
        yc = mean(y);
        dx = xc - center(1);
        dy = yc - center(2);
        rho = norm_factor * exp(-(dx^2 + dy^2)/(sigma^2));
        
        % Local load vector (each node gets Area/3 * ρ)
        be = (Area * rho / 3) * ones(3,1);
        b(ids) = b(ids) + be;
    end
    
    %% =====================================================================
    % 4. APPLY DIRICHLET BOUNDARY CONDITIONS (φ = 0 on boundary)
    % =====================================================================
    % Identify boundary nodes
    bottom = 1:(Nx+1);
    top    = (N - Nx):N;
    left   = 1:(Nx+1):N;
    right  = (Nx+1):(Nx+1):N;
    
    boundary_nodes = unique([bottom, top, left, right]);
    free_nodes = setdiff(1:N, boundary_nodes);
    
    fprintf('Boundary nodes: %d, Free nodes: %d\n', ...
            length(boundary_nodes), length(free_nodes));
    
    % Extract system for free DOFs
    A_ff = A(free_nodes, free_nodes);
    b_f  = b(free_nodes);
    
    % Check matrix conditioning
    cond_A = condest(A_ff);
    fprintf('Condition number of A_ff: %.2e\n', cond_A);
    
    %% =====================================================================
    % 5. SOLVE THE LINEAR SYSTEM
    % =====================================================================
    phi = zeros(N, 1);            % Initialize full potential vector
    phi(free_nodes) = A_ff \ b_f; % Solve only for interior nodes
    
    %% =====================================================================
    % 6. ANALYTICAL SOLUTION (φ = 1 / max(|r - r0|, R0))
    % =====================================================================% === Analytical solution ===
    % === Analytical solution ===
    % === Analytical potential ===
    r = sqrt((nodes(:,1)-center(1)).^2 + (nodes(:,2)-center(2)).^2);
    phi_ana = 1 ./ sqrt(r.^2 + R0^2);
    
    Phi_ana = reshape(phi_ana, Ny+1, Nx+1);
    
    % === Clipping only for visualization ===
    phi_display = Phi_ana;
    phi_clip_value = 1e13;   % you can tune this; try 1e13–1e14 range
    phi_display(phi_display > phi_clip_value) = phi_clip_value;
    
    
    
    
    
    
    
    %% =====================================================================
    % 7. RESHAPE FOR VISUALIZATION (match meshgrid structure)
    % =====================================================================
    
    % Reshape to match FEM mesh orientation
    % === Reshape for plotting ===
    % Both FEM and analytical fields must match Xg, Yg orientation
    Phi     = reshape(phi,     Ny+1, Nx+1);
    Phi_ana = reshape(phi_ana, Ny+1, Nx+1);
    
    
    
    
    %% =====================================================================
    % 8. VISUALIZATION: FEM, Analytical, and Error
    % =====================================================================
    figure('Position', [100, 100, 1400, 500]);
    
    % FEM solution
    subplot(1,3,1);
    surf(Xg, Yg, Phi, 'EdgeColor', 'none');
    colorbar; view(3); axis equal tight;
    title('FEM Solution (Uniform Mesh)');
    xlabel('x'); ylabel('y'); zlabel('\phi');
    hold on;
    plot3(center(1), center(2), max(phi)*1.3, 'rp', ...
        'MarkerFaceColor','r','MarkerSize',10);
    text(center(1), center(2), max(phi)*1.4, 'H^+', ...
        'Color','r','FontSize',12,'HorizontalAlignment','center');
    
    subplot(1,3,2);
    surf(Xg, Yg, phi_display, 'EdgeColor', 'none');
    colorbar;
    view(35,40);
    axis equal tight;
    zlim([0, phi_clip_value]);
    title(['Analytical (clipped at ', sprintf('%.1e',phi_clip_value), ')']);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('\phi');
    
    % Error
    subplot(1,3,3);
    Error = abs(Phi - Phi_ana);
    surf(Xg, Yg, log10(Error+1e-12), 'EdgeColor','none');
    title('log_{10}(Error)');
    xlabel('x'); ylabel('y');
    axis equal tight; view(30,35);
    colorbar;
    %% =====================================================================
    % 9. CROSS-SECTION THROUGH PROTON CENTER (aligned with mesh)
    % =====================================================================
    % =====================================================================
    % 9. CROSS-SECTION THROUGH PROTON CENTER
    % =====================================================================
    line_x = linspace(0, L, 500);
    line_y = center(2) * ones(size(line_x));
    
    line_fem = interp2(Xg, Yg, Phi, line_x, line_y, 'linear', 0);
    r = sqrt((line_x - center(1)).^2 + (line_y - center(2)).^2);
    line_ana = min(1 ./ max(r, R0), 1/R0); % clipped analytical
    
    figure;
    plot(line_x, line_fem, 'b-', 'LineWidth', 2); hold on;
    plot(line_x, line_ana, 'r--', 'LineWidth', 2);
    xlim([0 L]); ylim([0, 10]); % clip y-axis for visibility
    xlabel('x'); ylabel('\phi');
    legend('FEM', 'Analytical (clipped 1/r)', 'Location', 'best');
    title('Cross-section through Proton Center');
    grid on;
    
    
    
    %% =====================================================================
    % 10. ERROR ANALYSIS (away from singularity)
    % =====================================================================
    valid = (nodes(:,1) > 0.1*L) & (nodes(:,1) < 0.9*L) & ...
            (nodes(:,2) > 0.1*L) & (nodes(:,2) < 0.9*L) & ...
            (sqrt((nodes(:,1)-center(1)).^2 + (nodes(:,2)-center(2)).^2) > 5*R0);
    
    err = abs(phi(valid) - phi_ana(valid));
    L2_err  = sqrt(mean(err.^2));
    Linf_err = max(err);
    
    fprintf('\nERROR METRICS (interior, away from singularity):\n');
    fprintf('  L2 error   = %.3e\n', L2_err);
    fprintf('  L∞ error   = %.3e\n', Linf_err);
    fprintf('\nComputation completed successfully.\n');
