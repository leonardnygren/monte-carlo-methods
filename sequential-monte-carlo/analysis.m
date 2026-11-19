%% Question 3
max_n = 10;    
N = 100000;      
estimated_cn = zeros(max_n, 1);
for n = 1:max_n
    N_SA = 0;  % Set Counter for SAWs
    
    parfor i = 1:N  
        % Generate a random walk of n steps
        pos = zeros(n+1, 2);  
        current = [0, 0];
        pos(1, :) = current;
        is_SAW = true;
        
        for k = 1:n
            % Random direction (1=up, 2=down, 3=left, 4=right)
            dir = randi(4);
            switch dir
                case 1  
                    step = [0, 1];
                case 2
                    step = [0, -1];
                case 3
                    step = [-1, 0];
                case 4
                    step = [1, 0];
            end
            current = current + step;
            pos(k+1, :) = current;
        end
        
     
        if size(unique(pos, 'rows'), 1) == n + 1   % Check if all positions are unique
            N_SA = N_SA + 1;
        end
    end
    
    % Estimate c_n(2) as (N_SA / N) * 4^n
    estimated_cn(n) = (N_SA / N) * 4^n;
    fprintf('n = %d, Estimated c_n = %.2f\n', n, estimated_cn(n));
end
figure;
plot(1:max_n, estimated_cn, 'o-', 'LineWidth', 2);
xlabel('Walk length (n)');
ylabel('Estimated c_n(2)');
title('SIS Estimate of Self-Avoiding Walk Counts');
grid on;

%% Question 4
max_n = 18;
N = 100000;

weights = zeros(N, max_n);

parfor i = 1:N

    % Temporary weight vector for this iteration
    local_weights = zeros(1, max_n);

    current = [0, 0];

    visited = containers.Map('KeyType', 'char', ...
        'ValueType', 'logical');

    visited(sprintf('%d,%d', current(1), current(2))) = true;

    omega = 1;  % omega_0 = 1

    for k = 1:max_n

        neighbors = current + [0,1; 0,-1; 1,0; -1,0];
        available = [];

        % Check unvisited neighbors
        for j = 1:4
            pos = neighbors(j, :);
            key = sprintf('%d,%d', pos(1), pos(2));
            if ~isKey(visited, key)
                available = [available; pos];
            end
        end

        d_k = size(available, 1);
        if d_k == 0
            % Trapped: remaining weights stay zero
            break;
        else
            omega = omega * d_k;
            local_weights(k) = omega;

            % Randomly select an available neighbor
            chosen_idx = randi(d_k);
            new_pos = available(chosen_idx, :);

            key = sprintf('%d,%d', new_pos(1), new_pos(2));
            visited(key) = true;

            current = new_pos;
        end
    end

    % Assign the complete row (valid parfor slicing)
    weights(i,:) = local_weights;
end

% Calculate estimated c_n(2)
estimated_cn = mean(weights, 1);

% Plot results
figure;
semilogy(1:max_n, estimated_cn, 'o-', 'LineWidth', 2);
xlabel('Walk Length (n)');
ylabel('Estimated c_n(2)');
title('SIS for SAWs (Validated Implementation)');
grid on;
%% Question 5: Sequential Importance Sampling with Resampling (SISR)

clear;
clc;

% Parameters
max_n = 18;
N = 100000;
resample_threshold = 0.8;

rng(2024);
% Initialize particles
% Each particle has:
% pos     - current position
% visited - value-based matrix of previously visited positions
% weight  - importance weight

particles = repmat(struct( ...
    'pos', [0, 0], ...
    'visited', [0, 0], ...
    'weight', 1), N, 1);

% Estimated self-avoiding walk counts
estimated_cn = zeros(1, max_n);

% Accumulated normalization constant
scale = 1;

% Number of resampling events
resampling_events = 0;

% Main SISR loop

for k = 1:max_n

    % Mutation step

    for i = 1:N

        if particles(i).weight == 0
            continue;
        end

        current = particles(i).pos;
        visited = particles(i).visited;

        % Four possible neighboring positions
        neighbors = [current + [0, 1];
                     current + [0, -1];
                     current + [1, 0];
                     current + [-1, 0]];

        % Find unvisited neighboring positions
        available = neighbors( ...
            ~ismember(neighbors, visited, 'rows'), :);

        d_k = size(available, 1);

        if d_k == 0

            % Particle is trapped
            particles(i).weight = 0;

        else

            % Select a uniformly random available neighbor
            chosen_idx = randi(d_k);
            new_pos = available(chosen_idx, :);

            % Update particle position and history
            particles(i).pos = new_pos;
            particles(i).visited = [visited; new_pos];

            % Update importance weight
            particles(i).weight = ...
                particles(i).weight * d_k;

        end
    end

    % Estimate c_n(2)

    weights = [particles.weight]';

    % Current normalizing constant estimate
    estimated_cn(k) = scale * mean(weights);

    % Check for extinction
    sum_weights = sum(weights);

    if sum_weights == 0
        fprintf('All particles extinct at n = %d\n', k);
        estimated_cn(k:end) = 0;
        break;
    end

    % Calculate Effective Sample Size (ESS)

    normalized_weights = weights / sum_weights;

    ESS = 1 / sum(normalized_weights.^2);

    ESS_ratio = ESS / N;

    % Resampling step

    if ESS < resample_threshold * N

        % Preserve accumulated normalizing constant
        scale = estimated_cn(k);

        % Systematic resampling
        cumulative_weights = cumsum(normalized_weights);
        cumulative_weights(end) = 1;

        % Generate systematic sampling points
        u = ((0:N-1)' + rand) / N;

        % Efficient one-pass resampling
        idx = zeros(N, 1);

        i = 1;

        for j = 1:N

            while u(j) > cumulative_weights(i)
                i = i + 1;
            end

            idx(j) = i;

        end

        % Select new particles
        particles = particles(idx);

        % Reset importance weights
        [particles.weight] = deal(1);

        % Update resampling counter
        resampling_events = resampling_events + 1;

    end

    % Display results

    fprintf(['n = %2d | Estimated c_n = %.4f | ' ...
             'ESS/N = %.4f | Resampling events = %d\n'], ...
             k, estimated_cn(k), ESS_ratio, resampling_events);

end

% Results

fprintf('\nTotal resampling events: %d\n', resampling_events);

% Plot estimated self-avoiding walk counts

figure;

semilogy(1:max_n, estimated_cn, 'o-', ...
    'LineWidth', 2, 'MarkerSize', 6);

xlabel('Walk Length (n)');
ylabel('Estimated c_n(2) (log scale)');

title('Sequential Importance Sampling with Resampling');

grid on;



%% Question 6: Estimate A_2, mu_2 and gamma_2 using SISR
% Repeat the SISR estimation independently to assess Monte Carlo variation.
% The SISR method is the same adaptive-ESS approach as in Question 5,
% generalized in the sisr_saw_counts function at the end of this file.
max_n = 18;
N = 10000;
num_iterations = 3;
resample_threshold = 0.8;
rng(2025);  % Reproducible independent replications

all_estimated_cn = zeros(num_iterations, max_n);
estimated_params = zeros(num_iterations, 3);
n = (1:max_n)';
Xreg = [ones(max_n, 1), n, log(n)];

for iter = 1:num_iterations
    % Sequential importance sampling WITH resampling, not ordinary SIS.
    [cn, resampling_events] = sisr_saw_counts( ...
        2, max_n, N, resample_threshold);
    all_estimated_cn(iter, :) = cn;

    % Fit log(c_n) = log(A_2) + n*log(mu_2) + (gamma_2-1)*log(n).
    beta = Xreg \ log(cn(:));
    A2 = exp(beta(1));
    mu2 = exp(beta(2));
    gamma2 = beta(3) + 1;
    estimated_params(iter, :) = [A2, mu2, gamma2];

    fprintf('Question 6, iteration %d (resampling events: %d):\n', ...
        iter, resampling_events);
    fprintf('  A_2 = %.6f, mu_2 = %.6f, gamma_2 = %.6f\n', ...
        A2, mu2, gamma2);
end

fprintf('Question 6 mean [A_2, mu_2, gamma_2]:\n');
disp(mean(estimated_params, 1));
fprintf('Question 6 sample SD [A_2, mu_2, gamma_2]:\n');
disp(std(estimated_params, 0, 1));


%% Question 9: Estimate asymptotic parameters using SISR in higher dimensions
% Use one generalized implementation rather than separate SIS code blocks.
max_n = 18;
N = 10000;
resample_threshold = 0.8;
dimensions = [5, 10];
rng(2026);

n = (1:max_n)';
Xreg = [ones(max_n, 1), n, log(n)];
q9_estimated_params = zeros(numel(dimensions), 3);

for dim_idx = 1:numel(dimensions)
    d = dimensions(dim_idx);
    [cn, resampling_events] = sisr_saw_counts( ...
        d, max_n, N, resample_threshold);

    beta = Xreg \ log(cn(:));
    Ad = exp(beta(1));
    mud = exp(beta(2));
    gammad = beta(3) + 1;
    q9_estimated_params(dim_idx, :) = [Ad, mud, gammad];

    % Large-dimensional asymptotic approximation from the assignment.
    mu_asymptotic = (2*d - 1) - 1/(2*d) ...
        - 3/(2*d)^2 - 16/(2*d)^3;

    fprintf('\nQuestion 9, d = %d (resampling events: %d):\n', ...
        d, resampling_events);
    fprintf('  A_%d = %.6f\n', d, Ad);
    fprintf('  mu_%d = %.6f (theoretical bounds: %d <= mu_%d <= %d)\n', ...
        d, mud, d, d, 2*d-1);
    fprintf('  gamma_%d = %.6f\n', d, gammad);
    fprintf('  Large-d approximation to mu_%d: %.6f\n', ...
        d, mu_asymptotic);
end


%% Task 10

load('population_2024.mat');
%parameters
A = 0.8; B = 3.8; C = 0.6; D = 0.99; G = 0.8; H = 1.25;
n = 100; % generations
N_vals = [500, 1000, 10000]; % particle sizes

% (uniform) observation density
p = @(x, y) (y >= G * x & y <= H * x) / (H * x - G * x);

% Loop over different particle sizes (N)
for N_idx = 1:length(N_vals)
    N = N_vals(N_idx);
    
    %initialize vectors for filter estimates, intervals and weights
    tau = zeros(1, n + 1);
    taulower = zeros(1, n + 1);
    tauupper = zeros(1, n + 1); 
    w = zeros(N, 1);
    part = C + (D - C) * rand(N, 1); % particles
    
    
    w = arrayfun(@(x) p(x, Y(1)), part); % weights (for first observation)
    
    if sum(w) == 0
        error('All particle weights are zero.');
    end

    tau(1) = sum(part .* w) / sum(w); % estimate (for first generation)
    
    % confidence interval for generation 0
    [xx, I] = sort(part); % sort data (particles)
    cw = cumsum(w(I)) / sum(w); % cumulative normalized weightsum
    Ilower = find(cw >= 0.025, 1); % index for lower 2.5% quantile
    Iupper = find(cw >= 0.975, 1); % index for upper 2.5% quantile
    taulower(1) = xx(Ilower); % lower 2.5% quantile
    tauupper(1) = xx(Iupper); % upper 2.5% quantile
    
    % Resample particles using weights from the first observation
    ind = randsample(N, N, true, w);
    part = part(ind);


    % Particle filter main loop
    for k = 1:n
        Rk = A + (B - A) * rand(N, 1);
        part = Rk .* part .* (1 - part);
        
        w = arrayfun(@(x) p(x, Y(k + 1)), part); % weights for the current observation
        
        if sum(w) == 0
            error('All particle weights are zero.');
        end

        tau(k + 1) = sum(part .* w) / sum(w);  % estimate filtered expectation 
        
        %confidence intervals
        [xx, I] = sort(part); % Sort particles
        cw = cumsum(w(I)) / sum(w); % cumulative normalized weightsum
                                    % for sorted data
        Ilower = find(cw >= 0.025, 1); % index for lower 2.5% quantile
        Iupper = find(cw >= 0.975, 1); % index for upper 2.5% quantile
        taulower(k + 1) = xx(Ilower); % lower 2.5% quantile
        tauupper(k + 1) = xx(Iupper); % upper 2.5% quantile
        
        %resampling indices and particles
        ind = randsample(N, N, true, w); 
        part = part(ind); 
    end
    
    % save results for each N in a different tau variable
    if N == 500
        tau_500 = tau;
        taulower_500 = taulower;
        tauupper_500 = tauupper;
    elseif N == 1000
        tau_1000 = tau;
        taulower_1000 = taulower;
        tauupper_1000 = tauupper;
    elseif N == 10000
        tau_10000 = tau;
        taulower_10000 = taulower;
        tauupper_10000 = tauupper;
    end
    
    % Display results (estimate, confidence intervals, and if true value is inside the interval)
    disp(['Results for N = ', num2str(N), ' particles:']);
    total_diff = 0; % variable accumulated total absolute difference
    for k = 1:n + 1
        diff = abs(tau(k) - X(k));  % Absolute difference
        % Check if true value is inside confidence interval
        inside_interval = (X(k) >= taulower(k) && X(k) <= tauupper(k)); % Check if true value is inside the CI
        
        % display current generation results
        fprintf('Generation %d: Estimate = %.4f, True Value = %.4f, CI = [%.4f, %.4f], Inside CI = %d\n', ...
            k, tau(k), X(k), taulower(k), tauupper(k), inside_interval);
        
        % get total difference
        total_diff = total_diff + diff;
    end
    
    % display total difference
    fprintf('Total absolute difference for N = %d: %.4f\n\n', N, total_diff);
end
%%
% plotting (one plot for each N) of estimates vs real values
figure()
hold on
plot(0:n, tau_500, '--o','LineWidth', 1);
plot(0:n, X, '--x', 'LineWidth', 1);
hold off
xlabel('generation')
ylabel('relative population size')
legend('Estimation','X')
title('Filter estimate')

figure()
hold on
plot(0:n, tau_1000, '--o','LineWidth', 1);
plot(0:n, X, '--x', 'LineWidth', 1);
hold off
xlabel('generation')
ylabel('relative population size')
legend('Estimation','X')
title('Filter estimate')

figure()
hold on
plot(0:n, tau_10000, '--o','LineWidth', 1);
plot(0:n, X, '--x', 'LineWidth', 1);
hold off
xlabel('generation')
ylabel('relative population size')
legend('Estimation','X')
title('Filter estimate')



%% Shared local functions (Questions 6 and 9)
% Local functions must appear after all script code in a MATLAB .m script.
function [estimated_cn, resampling_events] = sisr_saw_counts(d, max_n, N, resample_threshold)
    % Adaptive SISR for counts of self-avoiding walks on Z^d.
    % Each particle stores its complete history as a numeric matrix.
    % Weights accumulate until ESS-based resampling occurs; after resampling
    % the normalizing-constant estimate is preserved in 'scale'.

    directions = [eye(d); -eye(d)];
    paths = repmat({zeros(1, d)}, N, 1);
    weights = ones(N, 1);
    estimated_cn = zeros(1, max_n);
    scale = 1;
    resampling_events = 0;

    for k = 1:max_n
        % Mutation: extend each live walk by a uniformly chosen free neighbor.
        for i = 1:N
            if weights(i) == 0
                continue;
            end

            path = paths{i};
            neighbors = path(end, :) + directions;
            free = ~ismember(neighbors, path, 'rows');
            available = neighbors(free, :);
            dk = size(available, 1);

            if dk == 0
                weights(i) = 0;
            else
                next_pos = available(randi(dk), :);
                paths{i} = [path; next_pos];
                weights(i) = weights(i) * dk;
            end
        end

        total_weight = sum(weights);
        if total_weight == 0
            error('All SISR particles trapped at step %d; increase N.', k);
        end

        % Normalizing-constant estimator, valid across resampling events.
        estimated_cn(k) = scale * (total_weight / N);

        % Resample only if ESS drops below the chosen threshold.
        % No resampling is needed after the last step.
        if k < max_n
            normalized_weights = weights / total_weight;
            ESS = 1 / sum(normalized_weights.^2);
            if ESS < resample_threshold * N
                cumulative_weights = cumsum(normalized_weights);
                cumulative_weights(end) = 1;
                u = ((0:N-1)' + rand) / N;

                % One-pass systematic resampling.
                idx = zeros(N, 1);
                j = 1;
                for i = 1:N
                    while u(i) > cumulative_weights(j)
                        j = j + 1;
                    end
                    idx(i) = j;
                end

                % Cell arrays contain value data; duplicate paths remain
                % independent when subsequently extended (copy-on-write).
                paths = paths(idx);
                weights = ones(N, 1);
                scale = estimated_cn(k);
                resampling_events = resampling_events + 1;
            end
        end
    end
end
