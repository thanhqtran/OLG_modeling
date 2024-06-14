# Some Large-Scale OLG
Learning how to model Large-Scale OLG models (life cycle) in Julia. Some papers:

- [Auerbach & Kotlikoff (1987)](https://kotlikoff.net/wp-content/uploads/2019/03/Dynamic-Fiscal-Policy_1.pdf)
- [Aiyagari (1994)](http://drphilipshaw.com/AyagariQJE94.pdf)
- [Hubbard, Skinner & Zeldes (1995)](https://doi.org/10.1086/261987)
- [Imrohoroglu, Imorohoroglu & Joines (1995)](https://doi.org/10.1007/BF01213942)
- [Huggett (1996)](http://drphilipshaw.com/Huggett%201996.pdf)
- [Carroll (1997)](https://doi.org/10.1162/003355397555109)
- [Krusell & Smith (1998)](http://www.econ.yale.edu/smith/250034.pdf)
- [Conesa & Krueger (1999)](https://doi.org/10.1006/redy.1998.0039)
- [Attanasio, Low & Sanchez-Marcos (2006)](https://doi.org/10.1257/aer.98.4.1517)
- [Huggett, Ventura & Yaron (2006)](https://doi.org/10.1016/j.jmoneco.2005.10.013)
- [Braun and Joines (2014)](https://www.sciencedirect.com/science/article/pii/S0165188915000780)

**Books**
- [Heer & Maussner (2009), Dynamic General Equilibrium Modeling Computational Methods and Applications](https://www.uni-augsburg.de/en/fakultaet/wiwi/prof/vwl/heer/dsge-book/).
- [Heer (2019), Public Economics: The Macroeconomic Perspective](https://www.uni-augsburg.de/en/fakultaet/wiwi/prof/vwl/heer/pubec-book/)

More: [https://github.com/robertdkirkby/LifeCycleOLGReadingList](https://github.com/robertdkirkby/LifeCycleOLGReadingList) and [https://www.robertdkirkby.com/life-cycle-and-olg-models/](https://www.robertdkirkby.com/life-cycle-and-olg-models/)

# Running Julia Script
- [Cheat Sheet](https://cheatsheet.juliadocs.org/)
- After testing the functionality of the code in `.ipynb`, the entire code must be wrapped inside a `main()` function. Then, the script can be called and run in the terminal.
  ```julia
  using Packages

  function main()
    the code
  end

  @time main()

  ```
  Run Julia from the terminal
  ```shell
  cd "path_to_script_containing_folder"
  julia script.jl
  ```

# Running Dynare

## in Julia
(*) does not work all the time
Documentation [https://juliapackages.com/p/dynare](https://juliapackages.com/p/dynare)
- Installation
```julia
using Pkg
pkg"add Dynare"
```
- Running
```julia
using Dynare
```
(if there is an error related to 'OpenBLAS32' )

```julia
import LinearAlgebra, OpenBLAS32_jll
LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
```
- Invoke `.mod'` file
```julia
context = @dynare "path/model.mod";
```
The results are stored in the `context` structure.
- View results
The context structure is saved in the directory `<path to modfile>/<modfilenane>/output/<modfilename>.jld2`. It can be loaded with
```julia
using JLD2
DD = load("<path to modfile>/<modefilename>/output/<modefilename>.jld2")``
```
The IRF graphs are saved in `<path to modfile>/<modfilenane>/graphs`.

## in Matlab
Quick Start: [here](https://www.dynare.org/resources/quick_start/)
```matlab
addpath /Applications/Dynare/x.y/matlab
cd '/Users/USERNAME/work'
dynare example1
```
to edit
```matlab
edit example1.mod
```

# Snippets
- Nonlinear Solver Example

  ```julia
  # Pkg.add("NLsolve")
  using NLsolve
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
- Loop
  - while `loop`
  ```julia
  max_iter = 10
  i = 1
  tol = 1e-6
  err = 0.1
  while i <= max_iter && abs(err) > tol
    do something
    if abs(err) > tol
      i += 1
    else
      break
    end
  end
  ```
  - for `loop` forward
  ```julia
  for s in 1:1:60
    k[s] = solve(k[s-1)
  end
  ```
  - for `loop` backward
  ```julia
  for s in 60:-1:1
    k[s] = solve(k[s+1])
  end
  ```
  
