# The AK60 model

This note is based on Heer (2009)'s DGEM. You should buy the book and read more in-depth analysis of this problem. The note here is very preliminary, prone to errors, and mainly used to explain the code.

## Households

No population growth. The representative household's life span is 60, such that

$$T + TR = 40 + 20 = 60,$$

where $T$ is the working length and $TR$ is the retirement length. Labor supply $n^s_t$ follows

$$l^s_t = 1 - n^s_t \text{ for } t \in \{1, 2, \dots, 40\},$$

$$l^s_t = 1 \text{ for } t \in \{41, 42, \dots, 60\},$$

where $l^s_t$ is leisure. After T years, retirement is mandatory. 

The agent's maximization problem is
```math
\sum^{T+TR}_{s=1} \beta^{s-1} u(c^s_{s+t-1}, l^s_{t+s-1}),
```

where $\beta$ is the discount factor. The instantaneous utility function:
$$u(c,l) = \frac{((c+\psi)l^\gamma)^{1-\eta} - 1}{1-\eta}.$$

An agent is born without wealth and leaves no bequests upon death, thus $k^1_t = k^{61}_t = 0$. 

The real budget constraint of the working agent is given by

$$k^{s+1}_{t+1} = (1+r_t) k^s_t + (1-\tau_t) w_t n^s_t - c^s_t \text{ for } s = 1, \dots, T, $$

where $r_t, w_t$ are the interest and wage rates, while $\tau$ is the social security contribution tax. 

Once retired, the agents receive public pensions $b$ and no labor earnings. The budget constraint for a retiree is

$$k^{s+1}_{t+1} = (1+r_t) k^s_t + b - c^s_t \text{ for } s=T+1, \dots, TR$$

The consumption-leisure FOC is

$$\frac{u_l(c^s_t, l^s_t)}{u_c(c^s_t, l^s_t)} = \gamma \frac{c^s_t + \psi}{l^s_t} = (1-\tau_t) w_t.$$

Euler equation:

```math
\frac{1}{\beta} = \frac{u_c(c^{s+1}_{t+1}, l^{s+1}_{t+1})}{u_c(c^s_t, l^s_t)}(1+r_{t+1}).
```

To derive the FOCs, we can solve by forming a Lagrangian equation or Bellman's Recursive equation.

## Production

Cobb-Douglas technology


$$Y_t = N_t^{1-\alpha} K_t^{\alpha}.$$


Under a perfectly competitive market, the factor prices are
$$w_t = (1-\alpha) K_t^\alpha N_t^{1-\alpha},$$

$$r_t = \alpha K_t^{\alpha-1} N_t^{-\alpha} - \delta.$$

## Government

The government funds a social security program via taxation. Its budget is balanced every period such that

$$\tau w_t N_t = \frac{TR}{T+TR}b.$$

## Equilibrium

Given policy $b$ and initial distribution of capital $\{k_0\}^{T+TR}_{s=1}$ , a collection of policy rules $c^s (k^s_t, K_t, N_t)$ $n^s (k^s_t, K_t, N_t)$, $k^{s+1} (k^s_t, K_t, N_t)$, relative prices $\{ w_t, r_t \}$ such that

- Individual and aggregate behavior are consistent

```math
N_t = \sum^T_{s=1} \frac{n^s_t}{T+TR},
```

```math
K_t = \sum^T_{s=1} \frac{k^s_t}{T+TR}.
```

- Relative prices solve the firm's optimization problem.

- Given input prices and government policy, $c^s_t, n^s_t, k^s_{t+1}$ solve the households' problem.

- Goods market clears:

```math
Y_t = N_t^{1-\alpha} K_t^\alpha = \sum^{T+TR}_{s=1} \frac{c^s_t}{T+TR} + K_{t+1} - (1-\delta) K_t.
```

- The government's budget is balanced.

## Calibration

Parameters: $\eta = 2, \beta = 0.99, \alpha = 0.3, \delta = 0.1, \gamma=2, \psi=0,001$​

Replacement ratio:

```math
\xi =\frac{b}{(1-\tau)w \bar{n}} = 0.3
```

with $\bar{n}$ is the average labor supply. This function is used to calculate $\tau$.

## Steady State Computation

In the steady state, distribution is constant over generations

```math
\{ k^s_t\}^{60}_{s=1} = \{ k^s_{t+1}\}^{60}_{s=1} = \{ \bar{k}^s\}^{60}_{s=1}.
```

Aggregate capital and labor are also constant
```math
K_t = K_{t+1} = K, \\
N_t = N_{t+1} = N.
```
Prices are also constant: $w, r, \tau$.

**Algorithm (steps)**

1. Make initial guesses of the steady state values of $K$ and $N$
2. Compute $w, r,\tau$ that solve the firm's and government's problems.
3. Compute the decision rule for $c^s, k^s, n^s$ by backward iteration, given the initial and terminal capital stock $k^1 = k^{61}=0$.
4. Compute $K$ and $N$ again and calculate the error from the initial guesses.
5. Update $K$ and $N$ , return to step 2 until convergence.

**More on Step 3**:

For backward iteration, since $k^{61}=0 = k^1$, we only need to calculate the decision rule for $k^{59}, k^{58}, \dots, k^2$. 

For $t=59, \dots, 41$: Since the retiree does not work, we only need to deal with the Euler equation to determine $k$. To calculate the optimal $k^{59}$, we need the information of $k^{60}$ and $k^{61}$. Although $k^{61}$ is known, $k^{60}$ is not. Here, we provide the first 2 guesses for the first 2 iterations $k^{60}_1 = 0.15$, and $k^{60}_2 = 0.2$. Note: the subscript denotes the $i^{th}$ iteration.


For $t=40$, you need to find the optimal values of $n^{40} (n^{41}, k^{41})$, $k^{40}(k^{41}, k^{42}, n^{41}, n^{42})$. As this is the last period with a positive labor supply, we have $n^{41} = n^{42} = 0$. Capital holdings $k^{41}, k^{42}$ have been solved in the retiree's problem. Hence, you have 2 equations of 2 unknowns. We can use Julia's `NLSolve`.


For $t=39,\dots, 1$, we will do things similar to the previous steps. Note that when we calculate $k^1$, if it is different from the tolerance, we must update $K$ and $N$.


Updating $k^{60}$ using the Secant Method:

```math
k^{60}_{i} = k^{60}_{i-1} - \frac{k^{60}_{i-1} - k^{60}_{i-2}}{k^{1}_{i-1} - k^{1}_{i-2}} k^1_{i-1}
```

**More on steps 1 and 5:**

It is actually sufficient to make only a guess of $N$ since $K$ must scale up with it due to Cobb-Douglas technology. Since there is a target for the steady state $N$, it is a more appropriate approach than guessing randomly. In the code, $nbar$ is the aggregate $N$.

The update algorithm is simple:

```math
N^{guess}_{j} = \phi N^{guess}_{j-1} + (1-\phi) N^{out}_{j-1}
```

with $\phi$ is the learning rate, $j$ is the $j^{th}$ iteration. 

$N^{out}$ is the outcome of an iteration taking $N^{guess}$ as an input. 

The algorithm should update and converge at iteration $\hat{j}$ such that 

```math
N^{guess}_{\hat{j}} = N^{out}_{\hat{j}}.
```

