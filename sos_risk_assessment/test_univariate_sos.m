load test.mat


n_cases = size(collision_moments, 1);

% Result arrays.
risk_bounds = zeros(1, n_cases);
flags = zeros(1, n_cases);
solve_times = zeros(1, n_cases);

parfor i = 1:n_cases
    [rb, f, sol]  = univariate_sos_risk_bound(collision_moments(i, :), "sedumi");
    risk_bounds(i) = rb;
    flags(i) = f;
    solve_times(i) = sol.solvertime;
end

chebyshev_risk_bounds = zeros(1, n_cases);
for i = 1:n_cases
    if collision_moments(i, 1) <= 0
        chebyshev_risk_bounds(i) = nan;
    else
        chebyshev_risk_bounds(i) = (collision_moments(i, 2) - collision_moments(i, 1)^2)/collision_moments(i,2);
    end
end

plot(risk_bounds);
hold on;
plot(chebyshev_risk_bounds);

legend("SOS Order 4", "Chebyshev")