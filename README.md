# Large Scale OLG
Learning how to model Large-Scale OLG models (life cycle) in Julia and MATLAB.

- [Auerbach & Kotlikoff (1987)](https://kotlikoff.net/wp-content/uploads/2019/03/Dynamic-Fiscal-Policy_1.pdf)
- [Huggett (1996)](http://drphilipshaw.com/Huggett%201996.pdf)
- [Aiyagari (1994)](http://drphilipshaw.com/AyagariQJE94.pdf)
- [Krusell & Smith (1998)](http://www.econ.yale.edu/smith/250034.pdf)
- [French(2005)]
- [De Nardi et al (2017)]

# Working with Julia
- [Cheat Sheet](https://cheatsheet.juliadocs.org/)
- After testing the functionality of the code in .ipynb, the entire code must be wrapped inside a `main()` function so that the script can be called and run in the terminal
  ```
  using Packages

  function main()
    the code
  end

  @time main()

  ```
  - Nonlinear Solver Example

  ```
  # solving for x,y -> z, with parameters a,b,c
  function system_eq!(F, z, a, b, c)
    x, y = z
    F[1] = x^2 + y^2 - a + 2b
    F[2] = (x * y)^c - a - b^2
  end
  #specify values of parameters
  a = 1.0  
  b = 2.0
  c = 2.0
  z_guess = [0.0, 0.0]  # Initial guess for x and y
  # Solve the system of equations
  result = nlsolve((F, z) -> system_eq!(F, z, a, b, c), z_guess)
  # Extract the solution
  x_solution = result.zero[1]
  y_solution = result.zero[2]
  ```

  When applying this solver, the correct guess is very important. With backward iteration, remember to take the previously solved known value as the guess. For example, in Auerbach-Kotlikoff, the guess to solve `k[39]` should take the initial guess of `k`,`n` in `k[40]`, `k[41]`, and `n[40]`.
