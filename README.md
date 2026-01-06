# Exploring Horizontal Discretizations of the Primitive Equations on Triangular and Hexagonal Grids
Data and notebook to explore Fourier symbols and instabilities of the linearized primitive equations in the Eady configuration discretized on triangular and hexagonal grids.

## Description
The linearized multilayer primitive equations in the Eady configuration can be solved by the superposition of plane waves. A Fourier transform of the system in the horizontal coordinates yields an eigenvalue problem for the system matrix. The system matrix may be decomposed into the Fourier symbols of the respective linearized (differential) operators. Here the equations are discretized on a triangular or hexagonal grid. The Fourier symbols of the linearized operators are provided as expressions in the _Julia_ programming language in the directory [data/symbols](https://github.com/smaasz/TriEadyInstabilities/tree/main/data/symbols) for different discretization schemes. The [Pluto](https://plutojl.org/) notebook at [notebooks/symbols.jl](https://github.com/smaasz/TriEadyInstabilities/blob/main/notebooks/symbols.jl) may be used to explore the Fourier symbols, assemble the system matrix from the Fourier symbols, and the growth rates of the baroclinically and symmetrically unstable modes. The maximal growth rates as functions of the wavenumber have been computed for different discretization schemes and parameter settings. The results are saved in the CSV-file [data/simspub.csv](https://github.com/smaasz/TriEadyInstabilities/blob/main/data/simspub.csv) which can be further explored in the notebook. Note that there are two versions of the notebook. The difference is in the packages used for plotting. The notebook [notebooks/symbols.jl](https://github.com/smaasz/TriEadyInstabilities/blob/main/notebooks/symbols.jl) relies on [Plots.jl](https://github.com/JuliaPlots/Plots.jl) whereas [notebooks/symbols_makie.jl](https://github.com/smaasz/TriEadyInstabilities/blob/main/notebooks/symbols.jl) relies on [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl/tree/master/CairoMakie). The latter requires longer precompilation.

## Getting Started
### Dependencies
- [Julia](https://julialang.org/) (version 1.11.5 used to obtain the package dependency resolution in the notebooks)
- [Pluto.jl](https://plutojl.org/) package (version 0.20.6)
- several other Julia packages downloaded in the notebooks

### Installation and 
1. (*Only necessary once!*) __Install Julia__. It is recommended to use the cross-platform installer [juliaup](https://github.com/JuliaLang/juliaup) and follow the instructions in the [README](https://github.com/JuliaLang/juliaup/blob/main/README.md). In order to reproduce the same results please install the Julia version 1.11.5. Using `juliaup` this is down by running the following command in your terminal: ```juliaup add 1.11.5```.
2. __Run Julia__. If you used `juliaup` for the installation you can run the correct version by executing `julia +1.11.5` in your terminal.
3. (*Only necessary once!*) __Install Pluto__. 
```julia
julia> Pkg.add("Pluto")
```
4. __Run Pluto__. This will open the Pluto starting page in your default browser.
```julia
julia> using Pluto; Pluto.run()
```
5. __Open Notebook__. Paste the URL (`https://github.com/smaasz/TriEadyInstabilities/blob/main/notebooks/symbols.jl` or `https://github.com/smaasz/TriEadyInstabilities/blob/main/notebooks/symbols_makie.jl`) into the respective text box of the Pluto starting page.  

## Version History

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
