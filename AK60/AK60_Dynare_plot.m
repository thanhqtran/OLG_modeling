% Extract steady-state values from Dynare
k_vals = oo_.steady_state(1:59); % k2 to k60 (first 59 variables)
n_vals = oo_.steady_state(60:99); % n1 to n40 (variable 60th to 99th)

% Define age vector
ages = 1:60; % Age 1 to 60

% Prepend k1 = 0
k_vals = [0; k_vals]; 

% Append n41 to n60 = 0
n_vals = [n_vals; zeros(20,1)]; 

% Plot capital stock (k1 to k60)
figure;
subplot(2,1,1); % Upper plot
plot(ages, k_vals, 'LineWidth', 1.5);
xlabel('Age');
ylabel('Capital Stock (k)');
title('Capital Stock Over the Life Cycle');
grid on;

% Plot labor supply (n1 to n60)
subplot(2,1,2); % Lower plot
plot(ages, n_vals, 'LineWidth', 1.5);
xlabel('Age');
ylabel('Labor Supply (n)');
title('Labor Supply Over the Life Cycle');
grid on;
