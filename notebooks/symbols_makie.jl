### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c365fb96-df3d-11f0-305f-cd90975c054e
begin
	using CSV
	import DataFrames: subset, DataFrame, select, filter, Not, ByRow, sort!, innerjoin
	import Downloads
	import DrWatson: @unpack, @dict, savename, wload
	import HTTP
	using HypertextLiteral: @htl, @html_str
	import Latexify
	import LinearAlgebra: eigvals, kron!, I, Diagonal
	import NaNMath
	using CairoMakie
	using PlutoUI
	using Printf
	import Symbolics: @variables, simplify, expand, get_variables, value, Num, substitute, simplify_fractions, taylor, expand_derivatives, Differential
	import Symbolics
	import SymbolicUtils
	using UUIDs: uuid1
	using RuntimeGeneratedFunctions
end;

# ╔═╡ 26b93993-390f-4278-b7ee-67350d6970be
TableOfContents(title="📚 Table of Contents", indent=true, depth=3, aside=true)

# ╔═╡ 6e2e12ab-3a7c-411e-b74a-2f6967162fa4
md"""
## Allow Downloads
"""

# ╔═╡ cb581f1d-ab4f-4e28-a4eb-747b77315ad6
md"""
The data needed to run the notebook is contained in the git repository `InstabilityOnGrids` hosted on GitHub at: [smaasz/TriEadyInstabilities/data](https://github.com/smaasz/TriEadyInstabilities/tree/main/data).
"""

# ╔═╡ dc12b2fc-26d7-4a3b-8c72-ab24a1fe4d41
md"""
Download data: $(@bind isdownload CheckBox(default=false))
"""

# ╔═╡ eb07fe92-8237-4c5f-ac9c-8ebeb30ff1ab
md"## Grid"

# ╔═╡ d906b120-1595-4bb9-b40a-6e17900e7539
md"""Grid: $(@bind grid_t Select([:TriA => "triangular A grid", :TriB => "triangular B grid", :TriC => "triangular C grid", :HexC => "hexagonal C grid"]))"""

# ╔═╡ 7b13f90d-7a51-45bd-acf2-b02545e05bde
function scalarloc(grid)
	if grid == :TriA
		:vertex
	elseif grid == :TriB
		:vertex
	elseif grid == :TriC
		:cell
	elseif grid == :HexC
		:vertex
	end
end

# ╔═╡ eeeaca1b-7861-4e66-9d48-62430f9e49c4
function vectorloc(grid)
	if grid == :TriA
		:vertex
	elseif grid == :TriB
		:cell
	elseif grid == :TriC
		:edge
	elseif grid == :HexC
		:edge
	end
end

# ╔═╡ 34ec4989-0add-4f87-b679-393a286077dd
function nlattices(colpt)
	if colpt == :vertex
		1
	elseif colpt == :cell
		2
	elseif colpt == :edge
		3
	end
end

# ╔═╡ 0c4d6003-d3b1-4fdf-acc4-b5348b606825
md"""
$(Resource("https://github.com/smaasz/InstabilityOnGrids/blob/main/plots/hexlattices.jpeg?raw=true", :height=>500))
"""

# ╔═╡ 6054481d-0bba-40c7-baef-c72cfe02e632
md"""
__Left__: Primitive unit cell of a hexagonal lattice. The red circle is the representative of the hexagonal lattice in the unit cell and located at a vertex, $v$, of a triangle. The unit cell contains two triangle centers, $c^u$ and $c^d$, that form two hexagonal lattices represented by the blue and purple diamond respectively. Moreover, the unit cell contains three edge, $e^{(1)}$, $e^{(2)}$, and $e^{(3)}$, centers that form three hexagonal lattices which are represented by the orange, yellow, and green square respectively. __Right__: Excerpt of the triangular/hexagonal grid centered at a vertex. The hexagon surrounding the vertex is called the _dual control volume_ corresponding to the vertex and its corners are the three centers of the upward pointing triangles (blue) and the three centers of the downward pointing triangles (red). The labels of the corners denote the _phase shift_ of a plane wave with fixed wavevector, $\boldsymbol{k}$, with respect to the central vertex. The dependence on the wavevector is neglected in the notation.
"""

# ╔═╡ 80e53723-b027-49ec-a178-aefc6d563672
md"""
__scalar quantities__: $(scalarloc(grid_t) == :vertex ? "🔴" : "🔷, ♦️")

__vector quantities__: $(vectorloc(grid_t) == :cell ? "🔷, ♦️" : (vectorloc(grid_t) == :edge ? "🟩, 🟨, 🟧" : "🔴"))
"""

# ╔═╡ 11dbe7db-a8ac-4286-8b43-e352f857842f
begin
	db = nlattices(scalarloc(grid_t))
	du = nlattices(vectorloc(grid_t))
end;

# ╔═╡ 44db3ae2-fb94-48ad-a306-9a95627618e7
begin
	hmt_schemes = Dict(
		:TriA => [
			:standard => "Standard",
		],
		:TriB => [
			:asc => "advective form, streamline derivative on cells",
			:avi => "advective form, vector-invariant", 
			:fdv => "flux form, divergence on vertices", 
			:fdcre => "flux form, diverence on cells with reconstruction on edges"
		],
		:TriC => [
			:ICON => "ICON",
			:MICON => "MICON"
		],
		:HexC => [
			:MPAS => "MPAS"
		]
	)

	hst_schemes = Dict(
        :TriA => [
            :low => "low",
            :high => "high",
        ],
        :TriB => [
            :low => "low",
            :high => "high",
        ],
        :TriC => [
            :low => "low",
        ],
        :HexC => [
            :low => "low",
            :high => "high",
        ],
    )
	"hmt & hst schemes"
end

# ╔═╡ fe40d664-3137-47ea-90a0-03bf4058be5c
md"Scheme: $(@bind hmt_scheme Select(hmt_schemes[grid_t]))"

# ╔═╡ 4589f32c-f6dc-4805-94fc-424655fa5289
md"## Fourier Symols"

# ╔═╡ 4dcd3ed3-6fa5-4f54-a9fa-2f109095292f
@variables f₀ g N² Ri M² β θU 𝕂ᵘ 𝕂ᵇ H Nz k l a h z α₁ α₂ α₃

# ╔═╡ 97c50791-1c10-452b-8f7c-a2e70f71872b
md"""
Operators between two tensor products of a discrete horizontal with the vertical space will be denoted by an underlining.
For notational convenience the underlining of an operator between two discrete horizontal spaces denotes the tensor product of the operator with the identity on the vertical space (i.e. Kronecker product of the matrices).

The fully assembled system matrices for the considered discretizations on a triangular/hexagonal grid are summarized in the following equation.
The system matrix acts on the vector of amplitudes specifying the horizontal velocity fields on the different vertical layers, the buoyancy field on the different vertical layers, and the surface elevation.
"""

# ╔═╡ 5136efea-d24b-4016-8942-581a91552796
WideCell(md"""
```math
\mathsf{S}=
\begin{bmatrix}
\mathsf{G_x}\otimes\underline{U} + \mathsf{G_y}\otimes\underline{V} + \left( \mathsf{A_x^u}\otimes\underline{U_z} + \mathsf{A_y^u}\otimes\underline{V_z}\right) \mathsf{\underline{W}} + f_0 \mathsf{\underline{M}} + \mathsf{\underline{D^u}} & \mathsf{\underline{G}}~\mathsf{\underline{P}} & g  \mathsf{G\otimes 1_V}\\
N^2 \mathsf{\underline{W}} + \mathsf{A_x^b}\otimes\underline{B_x} + \mathsf{A_y^b}\otimes\underline{B_y}& \mathsf{\Gamma_x}\otimes\underline{U} + \mathsf{\Gamma_y}\otimes\underline{V} +  \mathsf{\underline{D^b}}& \underline{0}\\
\Delta_z \mathsf{D} \otimes \mathsf{1_V^T} & \mathsf{0 \otimes 1_V^T} & \mathsf{0}
\end{bmatrix}
```
""")

# ╔═╡ 5b815179-5799-4544-81f7-37ff80bb1442
md"""__Small wavenumber approximation__: $(@bind doapprox PlutoUI.CheckBox(default=false))"""

# ╔═╡ f7cec82a-4606-4463-8ce3-25db7683df02
md"""__Phase substitutions__: $(@bind dophasesubs PlutoUI.CheckBox(default=false))"""

# ╔═╡ 55e52182-cbc1-4af3-9141-c02580e77195
md"""
#### ``\mathsf{G_x}``
"""

# ╔═╡ b2ff6746-af69-4a09-8ef5-46dcc1776c05
md"""
#### ``\mathsf{G_y}``
"""

# ╔═╡ 3877ac7e-142d-4c3b-8b11-3574cd4b2f9c
md"""
#### ``\mathsf{A^u_x}``
"""

# ╔═╡ 31a8215d-64a0-44af-b54c-97d4de167ccd
md"""
#### ``\mathsf{A^u_y}``
"""

# ╔═╡ d44d24eb-e829-42a1-9372-1f54e2a56d0f
md"""
#### ``\mathsf{M}``
"""

# ╔═╡ 8609a240-a029-48ea-b29d-8e0587b1338c
md"""
#### ``\mathsf{D^u}``
"""

# ╔═╡ 1607ff67-b476-400a-b013-0a5a3111de7a
md"""
#### ``\mathsf{G}``
"""

# ╔═╡ 37a62d52-1154-4820-80fc-eb64439bec68
md"""
#### ``\mathsf{A^b_x}``
"""

# ╔═╡ 6f28619c-df35-4942-a66b-76a7d5ee90da
md"""
#### ``\mathsf{A^b_y}``
"""

# ╔═╡ 82b95f73-f209-4302-ac27-896d463393aa
md"""
#### ``\mathsf{Γ_x}``
"""

# ╔═╡ 1e33e3aa-9338-472d-b6c3-45b94a304f21
md"""
#### ``\mathsf{Γ_y}``
"""

# ╔═╡ 70a0b93e-a878-4d6a-9a42-e24ed96abe97
md"""
#### ``\mathsf{D^b}``
"""

# ╔═╡ bc57b52f-deba-48d4-a3ae-589c19748249
md"""
#### ``\mathsf{D}``
"""

# ╔═╡ 655ba582-8c5e-44eb-b711-6ec8369326a5
md"## Symbolic Rewriting"

# ╔═╡ 0b6025ec-0e21-4a6a-a604-6aa98eef2350
md"""
### Small wavenumber approximations
by Taylor expansion of the trigonometric functions around 0.
"""

# ╔═╡ d75cbf05-6c02-4622-93d5-9231231cdb3c
nantrigexpand = let
	function rcosexpand(n)
	    x = Symbolics.variable(:x)
	    Symbolics.@rule NaNMath.cos(~x) => substitute(taylor(cos(x), x, 0, 0:n), Dict([x => ~x]))
	end

	function rsinexpand(n)
	    x = Symbolics.variable(:x)
	    Symbolics.@rule NaNMath.sin(~x) => substitute(taylor(sin(x), x, 0, 0:n), Dict([x => ~x]))
	end

	n -> 	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcosexpand(n), rsinexpand(n)])))
end

# ╔═╡ 38313589-8fda-4c9e-bcf5-7ce2b24360a4
function smallkapprox(fsymbol)
	simplify.(fsymbol; rewriter=nantrigexpand(4)) .|> expand
end

# ╔═╡ 7c38c4dc-2711-48f7-a31f-14b972b0e447
md"""
### Cancelation of rational numbers in fractions
must be executed explicitly due to a bug(?) in the Symbolics.jl package.
"""

# ╔═╡ 1540c9d9-25ac-48d6-9c2c-fe7650d1b872
function simfrac(f)
	e = reim(f)
	ne = []
	for pe in e
		npe = 0
		for t in Symbolics.terms(numerator(pe))
			for f in Symbolics.factors(denominator(pe))
				t = t / f
			end
			npe += t
		end
		push!(ne, npe)
	end
	ne[1] + im * ne[2]
end

# ╔═╡ c8446c65-9992-40e9-be38-7fb6ec9b0f62
function displayfs(fs, tvarout, tvarin)
	s = if tvarout == :vector && tvarin == :scalar
		if vectorloc(grid_t) == :edge
			fs
		else
			s = permutedims(fs, (2, 1, 3))
			reshape(s, 2 * du, db)
		end
	elseif tvarout == :scalar && tvarin == :vector
		if vectorloc(grid_t) == :edge
			fs
		else
			s = permutedims(fs, (1, 3, 2))
			reshape(s, db, 2 * du)
		end
	elseif tvarout == :vector && tvarin == :vector
		if vectorloc(grid_t) == :edge
			fs
		else
			s = permutedims(fs, (2, 1, 4, 3))
			reshape(s, 2 * du, 2 * du)
		end
	else
		fs
	end
	s = expand.(s)
	s = substitute.(s, Ref(Dict(
		a^2 => 4//3 * h^2, 
		a^4 => (4//3 * h^2)^2,
		a^6 => (4//3 * h^2)^3,
		a^8 => (4//3 * h^2)^4,
		a^10 => (4//3 * h^2)^5,
		a^12 => (4//3 * h^2)^6,
		a^14 => (4//3 * h^2)^7,
		a^16 => (4//3 * h^2)^8,
	)))
	s = try # issue in Symbolics
		simplify.(s; expand=true, rewriter=simfrac)
	catch
		simplify.(s; expand=true, simplify_fractions=false)
	end
end

# ╔═╡ dfc6fcc3-7357-4e5e-9061-f2571934d2d2
md"""
### Phase substitutions
in order to simplify notation.
"""

# ╔═╡ 92b811da-ccc3-4781-a085-5c556ce889ee
function phasesubs(a, h, k, l)
    _α₁ = a*k/2 + h*l/3
    _α₂ = -a*k/2 + h*l/3
    _α₃ = -2//3*h*l
    αs = [(_α₁, α₁), (_α₂, α₂), (_α₃, α₃)]
    cs = [(1,0), (2,0), (1,-1), (3,0), (2,-1), (4,0), (3,-1), (2,-2), (1//2, 0), (1//2, -1//2), (1, -1//2), (3//2, 0)]
    phase_subs = Dict{Num, Num}()
    for (c₁, c₂) in cs
		for (_ϕ₁, ϕ₁) in αs
		    for (_ϕ₂, ϕ₂) in αs
				t = simplify(c₁*_ϕ₁ + c₂*_ϕ₂; expand=true)
				phase_subs[t] = c₁*ϕ₁ + c₂*ϕ₂
				t = simplify(-c₁*_ϕ₁ - c₂*_ϕ₂; expand=true)
				phase_subs[t] = -c₁*ϕ₁ - c₂*ϕ₂
		    end
		end
    end
    phase_subs
end

# ╔═╡ 689e54d9-2ecd-46cb-9082-efb76a0ad99d
function subphases(fsymbol, phases)
	rsin = Symbolics.@rule NaNMath.sin(~x) => NaNMath.sin(get(phases, ~x, ~x))
	rcos = Symbolics.@rule NaNMath.cos(~x) => NaNMath.cos(get(phases, ~x, ~x))
	rewriter = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcos, rsin])))
	simplify.(fsymbol; rewriter) .|> expand .|> simplify
end

# ╔═╡ daa9e2bf-a1eb-4628-a974-024c7a0d61a6
phases = phasesubs(a,h,k,l);

# ╔═╡ 7891b7da-ba36-417d-bdde-c84187463636
md"""
## Instability Analysis
"""

# ╔═╡ 1de0d073-bc11-4e4f-832c-f3b0e01b58f7
md"""
The maximal growth rates as functions of wavenumber have been computed for different parameter values. These parameters can be varied in the box at the bottom right corner to explore different settings.
"""

# ╔═╡ e8dafe17-c065-426c-a928-dcec43e962b0
resp_instab = !isdownload ? Missing : Downloads.download("https://raw.githubusercontent.com/smaasz/TriEadyInstabilities/refs/heads/main/data/simspub.csv");

# ╔═╡ 2ea23935-697b-46be-886c-2f166475a590
ddf = let
	df = CSV.read(resp_instab, DataFrame)
	df = select(df, Not(:Ks, :iωs), :Ks => ByRow(x->eval(Meta.parse(x))) => :Ks, :iωs => ByRow(x-> eval(Meta.parse(x))) => :iωs)
	df = filter(row -> row.hmt_scheme ≠ "MICON", df)
	sort!(df, [:grid_t, :hmt_scheme])
	df
end

# ╔═╡ e45490c6-59ea-4538-9aeb-25e1398a0e23
colordf = let
	colors = RGBAf.(Makie.wong_colors())
	push!(colors, RGBAf(Makie.to_color(:purple)))
	grid_ts = vcat([[grid_t for i in 1:length(hmt_schemes[grid_t])] for grid_t in [:HexC, :TriA, :TriB, :TriC]]...)
	push!(grid_ts)
	hmt_schemes = vcat([first.(hmt_schemes[grid_t]) for grid_t in [:HexC, :TriA, :TriB, :TriC]]...)
	push!(hmt_schemes)
	DataFrame("grid_t" => String.(grid_ts), "hmt_scheme" => String.(hmt_schemes), "color" => colors[1:length(grid_ts)])
end

# ╔═╡ 7499fa65-cb31-45d1-aa43-d248dc1cd066
md"## Assembling the system matrix"

# ╔═╡ 5fde0451-f79c-472e-a431-52b72eb2b703
function vertical_ops()
    @variables f₀ N² M² Ri θU β H
    Ū  = (z + H * (1//2 + β)) * -M²/f₀
    Ūz = expand_derivatives(Differential(z)(Ū))

    vops = Dict(
        :Ū  => Ū, 
        :U  => Ū * cos(θU),
        :V  => Ū * sin(θU),
        :Uz => Ūz * cos(θU),
        :Vz => Ūz * sin(θU),
        :Bx => M² * -sin(θU),
        :By => M² * cos(θU),
    )

    @assert isequal(vops[:Bx], f₀ * vops[:Vz])
    @assert isequal(vops[:By], -f₀ * vops[:Uz])

	M²sub = Dict(M² => √(f₀^2 * N²/ Ri))
	params = (z, f₀, N², Ri, θU, β, H)
    (;
     [Symbol("$(name)f") => Symbolics.build_function(substitute(vop, M²sub), params...; expression=Val{false}) for (name, vop) in pairs(vops)]...
    )
end

# ╔═╡ 02be2a80-0fa3-42bd-ba7b-e0a9ec893c3c
function createbuffer(::Val{grid_t}) where {grid_t}
    if grid_t == :TriA
        Dict(
            :Gx    => zeros(ComplexF64, 2, du, 2, du),
            :Gy    => zeros(ComplexF64, 2, du, 2, du),
            :M     => zeros(ComplexF64, 2, du, 2, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, db),
            :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, db),
            :G     => zeros(ComplexF64, 2, du, db),
            :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
            :I     => zeros(ComplexF64, db, db),
            :Av⁽ˣ⁾ => zeros(ComplexF64, db, 2, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, db, 2, du),
            :Γx    => zeros(ComplexF64, db, db),
            :Γy    => zeros(ComplexF64, db, db),
            :Dᵇ    => zeros(ComplexF64, db, db),
            :D     => zeros(ComplexF64, db, 2, du),
        )
    elseif grid_t == :TriB
        Dict(
            :Gx    => zeros(ComplexF64, 2, du, 2, du),
            :Gy    => zeros(ComplexF64, 2, du, 2, du),
            :M     => zeros(ComplexF64, 2, du, 2, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, db),
            :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, db),
            :G     => zeros(ComplexF64, 2, du, db),
            :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
            :I     => zeros(ComplexF64, db, db),
            :Av⁽ˣ⁾ => zeros(ComplexF64, db, 2, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, db, 2, du),
            :Γx    => zeros(ComplexF64, db, db),
            :Γy    => zeros(ComplexF64, db, db),
            :Dᵇ    => zeros(ComplexF64, db, db),
            :D     => zeros(ComplexF64, db, 2, du),
        )
    elseif grid_t == :TriC
        Dict(
            :Gx    => zeros(ComplexF64, du, du),
            :Gy    => zeros(ComplexF64, du, du),
            :M     => zeros(ComplexF64, du, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, du, db),
            :A⁽ʸ⁾  => zeros(ComplexF64, du, db),
            :G     => zeros(ComplexF64, du, db),
            :Dᵘ    => zeros(ComplexF64, du, du),
            :I     => zeros(ComplexF64, db, db),
            :Av⁽ˣ⁾ => zeros(ComplexF64, db, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, db, du),
            :Γx    => zeros(ComplexF64, db, db),
            :Γy    => zeros(ComplexF64, db, db),
            :Dᵇ    => zeros(ComplexF64, db, db),
            :D     => zeros(ComplexF64, db, du),
        )
    else
        Dict(
            :Gx    => zeros(ComplexF64, du, du),
            :Gy    => zeros(ComplexF64, du, du),
            :M     => zeros(ComplexF64, du, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, du, db),
            :A⁽ʸ⁾  => zeros(ComplexF64, du, db),
            :G     => zeros(ComplexF64, du, db),
            :Dᵘ    => zeros(ComplexF64, du, du),
            :I     => zeros(ComplexF64, db, db),
            :Av⁽ˣ⁾ => zeros(ComplexF64, db, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, db, du),
            :Γx    => zeros(ComplexF64, db, db),
            :Γy    => zeros(ComplexF64, db, db),
            :Dᵇ    => zeros(ComplexF64, db, db),
            :D     => zeros(ComplexF64, db, du),
        )
    end
end

# ╔═╡ 13394a42-764d-44d4-a730-66c5109647ef
b = createbuffer(Val(grid_t));

# ╔═╡ b05f634b-56a8-49c7-8a2b-d410d72582fe
function systemmat(grid_t::Union{Val{:TriA}, Val{:TriB}}, fsyms, b, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme, useidealized=Dict(), vfops=vertical_ops())
    Δz  = H / Nz
    h   = a * √3/2

    # conversion to dissipation parameters
    if dissip_scheme == :biharmonic
		𝕂ᵘ = Vᵘ * a^3
		𝕂ᵇ = Vᵇ * a^3
    else
		𝕂ᵘ = Vᵘ * a
		𝕂ᵇ = Vᵇ * a
    end

    # compute fourier symbols and store in buffer
    for (name, fsym) in pairs(fsyms)
        useideal = get(useidealized, name, false)
	fsym[2](b[name], 0.0, f₀, N², Ri, θU, β, k, l, useideal ? 1e-20 : a, useideal ? √3/2*1e-20 : h)
    end

    # vertical operators
    @unpack Ūf, Uf, Vf, Uzf, Vzf, Bxf, Byf = vfops
    Ū̲ = Diagonal([Ūf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲ = Diagonal([Uf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲ = Diagonal([Vf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗ = (U̲, V̲)
    U̲z = Diagonal([Uzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲z = Diagonal([Vzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗z = (U̲z, V̲z)
    B̲x = Diagonal([Bxf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲y = Diagonal([Byf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲⃗ = (B̲x, B̲y)
    
    W̲ = let
		M = Δz * 1/2 * [iV ≤ iVi ≤ iV+1 for iV=1:Nz, iVi=1:Nz+1] * [iV < iVi ? 1 : 0 for iVi=1:Nz+1, iV=1:Nz]
		[kron(-b[:D][:,jTH,:], M) for jTH=1:2]
    end
    
    P̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ? 1 : 0 for iV=1:Nz, iVi=1:Nz-1] * [iVi ≤ iV ≤ iVi+1 for iVi=1:Nz-1, iV=1:Nz]
		kron(I(db), -M)
    end

    # assembling	
    S̲ = zeros(ComplexF64, (2*du+db)*Nz+db, (2*du+db)*Nz+db)

    ru⃗ = [(iTH-1)*du*Nz+1:iTH*du*Nz for iTH=1:2]
    rb = 2*du*Nz+1:(2*du+db)*Nz
    rη = (2*du+db)*Nz+1:(2*du+db)*Nz+db
    
    # U⃗
    @views for iTH = 1:2
		for jTH = 1:2
		    kron!(S̲[ru⃗[iTH], ru⃗[jTH]], b[:Gx][iTH,:,jTH,:], U̲)
	        #kron!(S̲[ru⃗[iTH], ru⃗[jTH]], b[:Gx][iTH,:,jTH,:], Ū̲)
		    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:Gy][iTH,:,jTH,:], V̲)
		    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:A⁽ˣ⁾][iTH,:,:], U̲z) * W̲[jTH]
		    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:A⁽ʸ⁾][iTH,:,:], V̲z) * W̲[jTH]
		    S̲[ru⃗[iTH], ru⃗[jTH]] += f₀ * kron(b[:M][iTH,:,jTH,:], I(Nz))
		    S̲[ru⃗[iTH], ru⃗[jTH]] += 𝕂ᵘ * kron(b[:Dᵘ][iTH,:,jTH,:], I(Nz))
		end
	
		kron!(S̲[ru⃗[iTH], rb], b[:G][iTH,:,:], I(Nz))
		S̲[ru⃗[iTH], rb] *= P̲
	
		kron!(S̲[ru⃗[iTH], rη], b[:G][iTH,:,:], g*ones(Nz,1))
    end

    # b
    @views for jTH = 1:2
		S̲[rb, ru⃗[jTH]] += N² * kron(b[:I], I(Nz)) * W̲[jTH]
		S̲[rb, ru⃗[jTH]] += kron(b[:Av⁽ˣ⁾][:,jTH,:], B̲x)
		S̲[rb, ru⃗[jTH]] += kron(b[:Av⁽ʸ⁾][:,jTH,:], B̲y)
    end
    @views kron!(S̲[rb, rb], b[:Γx], U̲)
    @views S̲[rb, rb] += kron(b[:Γy], V̲)
    @views S̲[rb, rb] += 𝕂ᵇ * kron(b[:Dᵇ], I(Nz))

    # η
    @views for jTH = 1:2
		kron!(S̲[rη, ru⃗[jTH]], b[:D][:,jTH,:], Δz*ones(1,Nz))
    end

    S̲[rη, rη] .= U̲[end] * b[:Γx] + V̲[end] * b[:Γy]
    
    S̲
end

# ╔═╡ 22ad33e5-0641-40aa-9aec-fc4e9be6b4b8
function systemmat(grid_t::Union{Val{:TriC}, Val{:HexC}}, fsyms, b, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme, useidealized=Dict(), vfops=vertical_ops())
    Δz  = H / Nz
    h   = a * √3/2

    # conversion to dissipation parameters
    if dissip_scheme == :biharmonic
		𝕂ᵘ = Vᵘ * a^3
		𝕂ᵇ = Vᵇ * a^3
	    else
		𝕂ᵘ = Vᵘ * a
		𝕂ᵇ = Vᵇ * a
    end

    # compute fourier symbols
    for (name, fsym) in pairs(fsyms)
        useideal = get(useidealized, name, false)
	fsym[2](b[name], 0.0, f₀, N², Ri, θU, β, k, l, useideal ? 1e-20 : a, useideal ? √3/2*1e-20 : h)
    end

    # vertical operators
    @unpack Uf, Vf, Uzf, Vzf, Bxf, Byf = vfops
    U̲ = Diagonal([Uf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲ = Diagonal([Vf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗ = (U̲, V̲)
    U̲z = Diagonal([Uzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲z = Diagonal([Vzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗z = (U̲z, V̲z)
    B̲x = Diagonal([Bxf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲y = Diagonal([Byf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲⃗ = (B̲x, B̲y)

    W̲ = let
		M = Δz * 1/2 * [iV ≤ iVi ≤ iV+1 for iV=1:Nz, iVi=1:Nz+1] * [iV < iVi ? 1 : 0 for iVi=1:Nz+1, iV=1:Nz]
		kron(-b[:D], M)
    end
    
    P̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ? 1 : 0 for iV=1:Nz, iVi=1:Nz-1] * [iVi ≤ iV ≤ iVi+1 for iVi=1:Nz-1, iV=1:Nz]
		kron(I(db), -M)
    end

    # assembling	
    S̲ = zeros(ComplexF64, (du+db)*Nz+db, (du+db)*Nz+db)

    ru⃗, rb, rη = (1:du*Nz, du*Nz+1:(du+db)*Nz, (du+db)*Nz+1:(du+db)*Nz+db)
    # U⃗

    @views begin
		kron!(S̲[ru⃗, ru⃗], b[:Gx], U̲)
		S̲[ru⃗, ru⃗] += kron(b[:Gy], V̲)
		S̲[ru⃗, ru⃗] += kron(b[:A⁽ˣ⁾], U̲z) * W̲
		S̲[ru⃗, ru⃗] += kron(b[:A⁽ʸ⁾], V̲z) * W̲
		S̲[ru⃗, ru⃗] += f₀ * kron(b[:M], I(Nz))
		S̲[ru⃗, ru⃗] += 𝕂ᵘ * kron(b[:Dᵘ], I(Nz))
		
		kron!(S̲[ru⃗, rb], b[:G], I(Nz))
		S̲[ru⃗, rb] *= P̲
		
		kron!(S̲[ru⃗, rη], b[:G], g * ones(Nz,1))
    end

    # b
    @views begin
		S̲[rb, ru⃗] += N² * kron(b[:I], I(Nz)) * W̲
		S̲[rb, ru⃗] += kron(b[:Av⁽ˣ⁾], B̲x)
		S̲[rb, ru⃗] += kron(b[:Av⁽ʸ⁾], B̲y)
	
		kron!(S̲[rb, rb], b[:Γx], U̲)
		S̲[rb, rb] += kron(b[:Γy], V̲)
		S̲[rb, rb] += 𝕂ᵇ * kron(b[:Dᵇ], I(Nz))
    end

    # η
    @views begin
		kron!(S̲[rη, ru⃗], b[:D], Δz*ones(1,Nz))
		S̲[rη, rη] .+= U̲[end] * b[:Γx] + V̲[end] * b[:Γy]
    end
    
    S̲
end


# ╔═╡ 84865a4a-9752-44e4-9746-41674be387c7
html"""<hr>"""

# ╔═╡ 2c658919-9be3-43f5-934b-9dd11a415eec
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ 8e1bcfe0-e9c8-40b2-a671-f0795d4c4b52
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;

# ╔═╡ 2474940b-0b25-44fe-995e-7e6e1dc94536
begin
	floataside(md"""
Number of layers: $(@bind sNz PlutoUI.Slider([32, 64]; show_value=true))
			   
Flow direction θU: $(@bind sθU PlutoUI.Slider([0.0, π/6]; show_value=true))

Flow shift β (in -HM²/f₀): $(@bind sβ PlutoUI.Slider(0.0:0.5:0.5; default=0.0, show_value=true))

Grid parameter a (in m): $(@bind sa PlutoUI.Slider([6.25e3, 12.5e3]; default=6.25e3, show_value=true))

Richardson number Ri: $(@bind sRi Select([100.0, 0.5]; default=100//1))
			   
horiz. scalar transport scheme: $(@bind shst_scheme Select([:low => "low-order accurate", :high => "high-order accurate"]))

Biharmonic friction V (in m/s): $(@bind sV PlutoUI.Slider([0, 1e-3, 1e-2]; show_value=true))
""", top=340)
end

# ╔═╡ 41f3ff09-cd17-4691-bd7a-dce52f1a7e97
fsyms = !isdownload ? Missing : let
	RuntimeGeneratedFunctions.init(@__MODULE__)
	config = @dict(grid_t, hmt_scheme, hst_scheme=shst_scheme, dissip_scheme=:biharmonic)
	path   = joinpath(
		"https://github.com/smaasz/TriEadyInstabilities/raw/refs/heads/main/data/symbols/", 
		savename(config, "jld2")
	)
	resp_symbols = Downloads.download(path)
	data   = wload(resp_symbols)
	
	@unpack fsyms, = data
	for (name, fsym) in pairs(fsyms) # mistake in symbols
		fsyms[name] = (
			Meta.parse(replace(string(fsym[1]), "BigFloat" => "Real")),
			Meta.parse(replace(string(fsym[2]), "BigFloat" => "Real"))
		)
	end
	fsyms_generated = Dict([
		name => (@RuntimeGeneratedFunction(fsym[1]), @RuntimeGeneratedFunction(fsym[2])) for (name, fsym) in pairs(fsyms)
	])
end

# ╔═╡ c1570e2e-e7c7-4cce-80a4-bd0c7efaa0c2
fsymbols = !isdownload ? Missing : let
	params = (z, f₀, N², Ri, θU, β, k, l, a, h)
	fsymbols = Dict{Symbol, Array{Complex{Num}}}([
	    name => fsym[1](params...) for (name, fsym) in pairs(fsyms) 
	])
	for (name, fsymbol) in pairs(fsymbols)
		fsymbols[name] = try # issue in Symbolics
			simplify.(fsymbol; expand=true, rewriter=simfrac)
		catch
			simplify.(fsymbol; expand=true, simplify_fractions=false)
		end
	end
	if doapprox
		for (name, fsymbol) in pairs(fsymbols)
			fsymbols[name] = smallkapprox(fsymbol)
		end
	elseif dophasesubs
		for (name, fsymbol) in pairs(fsymbols)
			fsymbols[name] = subphases(fsymbol, phases)
		end
	end
	fsymbols
end;

# ╔═╡ 456e8985-3197-47c6-8658-e7f7bb8253e5
let
	fs = displayfs(fsymbols[:Gx], :vector, :vector)
end

# ╔═╡ 33004562-7af3-4fb8-98d2-f2aed9d1203e
let
	fs = displayfs(fsymbols[:Gy], :vector, :vector)
end

# ╔═╡ a0341bc4-7240-4d8c-a6c4-9929241851de
let
	fs = displayfs(fsymbols[:A⁽ˣ⁾], :vector, :scalar)
end

# ╔═╡ 9b901c1a-b28f-4cd8-a800-d5d8b38b478e
let
	fs = displayfs(fsymbols[:A⁽ʸ⁾], :vector, :scalar)
end

# ╔═╡ 854b6598-20b2-4aa2-9410-112cf703d96e
let
	fs = displayfs(fsymbols[:M], :vector, :vector)
end

# ╔═╡ 58ef8cec-6e7a-4491-aa00-e4516094925e
let
	fs = displayfs(fsymbols[:Dᵘ], :vector, :vector)
end

# ╔═╡ 159467de-d9aa-4718-a8e1-012a90fa93a1
let
	fs = displayfs(fsymbols[:G], :vector, :scalar)
end

# ╔═╡ a40fb713-9ed6-4ec1-8f16-928d2b9a65cd
let
	fs = displayfs(fsymbols[:Av⁽ˣ⁾], :scalar, :vector)
end

# ╔═╡ 7de0eff1-04f9-4b11-a201-def19335c853
let
	fs = displayfs(fsymbols[:Av⁽ʸ⁾], :scalar, :vector)
end

# ╔═╡ 51b7e139-b647-444e-bd4b-534ff7a28e78
let
	fs = displayfs(fsymbols[:Γx], :scalar, :scalar)
end

# ╔═╡ 177d8297-0064-4b4f-a04e-cddf3cc30786
let
	fs = displayfs(fsymbols[:Γy], :scalar, :scalar)
end

# ╔═╡ 772c653a-5ca8-4672-ae48-2086fff8c639
let
	fs = displayfs(fsymbols[:Dᵇ], :scalar, :scalar)
end

# ╔═╡ 007d3be0-7cc0-4aeb-b46f-6f894a4480ef
let
	fs = displayfs(fsymbols[:D], :scalar, :vector)
end

# ╔═╡ 89277803-a204-4371-b475-8a8ff7233d60
function plotinstabilities(df, Ri, grid_t, hmt_scheme; Nz, θU, β, a, Vᵘ, Vᵇ, hst_scheme, scatterpts=nothing)
	size   = Ri > 1 ? (1400, 610) : (1400, 610)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 2.5 : 2.5
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)

	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 title  = "Max. Growth Rate on" * (Ri > 1 ? "Baroclinic" : "Symmetric") * " Axis",
			 xlabel = "wavenumber / kₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==Nz, :β=>x->x.==β, :θU=>x->x.==θU, :a=>x->x.==a, :Vᵘ=>x->x.≈Vᵘ, :Vᵇ=>x->x.≈Vᵇ, :Ri=>x->x.==Ri, :grid_t=>x->Symbol.(x).==grid_t, :hmt_scheme=>x->Symbol.(x).==hmt_scheme, :hst_scheme=>x->Symbol.(x).==hst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
		scatter!(ax, scatterpts[1] ./ fₛ, scatterpts[2] .* (sqrt(N²) / abs(M²) ); color, label=nothing, strokecolor=:black, strokewidth=2, markersize=16)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈Ri))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 0f4ecc8e-d314-4ed4-bc87-bc28b0dba078
function plotinstabilities(df, Ri; Nz, θU, β, a, Vᵘ, Vᵇ, hst_scheme)
	size   = Ri > 1 ? (1400, 610) : (1400, 610)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 2.5 : 2.5
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)

f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 title  = "Max. Growth Rate on" * (Ri > 1 ? "Baroclinic" : "Symmetric") * " Axis",
			 xlabel = "wavenumber / kₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==Nz, :β=>x->x.==β, :θU=>x->x.==θU, :a=>x->x.==a, :Vᵘ=>x->x.≈Vᵘ, :Vᵇ=>x->x.≈Vᵇ, :Ri=>x->x.==Ri, :hst_scheme=>x->Symbol.(x).==hst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈Ri))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 7271061e-092f-4454-b90b-7d0d9295f815
WideCell(
	plotinstabilities(ddf, sRi; Nz=sNz, θU=sθU, β=sβ, a=sa, Vᵘ=sV, Vᵇ=sV, hst_scheme=shst_scheme)
)

# ╔═╡ 4af8a93f-8a2a-4531-8c54-495bf97b1f4f
md"""
__Wavenumber K__: $(@bind sK PlutoUI.Slider(range(1e-10, min(2/√3*π/6.25e3, 2/√3*π/sa)*1.1, 500); show_value=true))
"""

# ╔═╡ b240a2e2-f574-4275-b471-2313428b18eb
S̲ = let
	a   = sa
	f₀  = -1e-4
	g   = 1e9
	N²  = 1e-6
	Ri  = sRi
	θU  = sθU
	β   = sβ
	Vᵘ  = -sV
	Vᵇ  = -sV
	H   = 4000
	Nz  = sNz
	θ = (Ri > 1 ? 0 : π/2) + θU
	k   = sK * cos(θ)
	l   = sK * sin(θ)
	systemmat(Val(grid_t), fsyms, b, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme=:biharmonic)
end

# ╔═╡ 9a3027be-a9d9-47b6-81c3-2f2ea23776b7
begin
	mgr = real(eigvals(S̲)[end])
	md"""__Max. growth rate__: $(@sprintf("%2.4e", mgr))"""
end

# ╔═╡ b78217cf-652a-403c-bb65-3721f963b95a
WideCell(
	plotinstabilities(ddf, sRi, grid_t, hmt_scheme; Nz=sNz, θU=sθU, β=sβ, a=sa, Vᵘ=sV, Vᵇ=sV, hst_scheme=shst_scheme, scatterpts=([sK], [mgr]))
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
DrWatson = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
RuntimeGeneratedFunctions = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
UUIDs = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[compat]
CSV = "~0.10.15"
CairoMakie = "~0.15.8"
DataFrames = "~1.8.1"
DrWatson = "~2.19.1"
HTTP = "~1.10.19"
HypertextLiteral = "~0.9.5"
Latexify = "~0.16.10"
NaNMath = "~1.1.3"
PlutoUI = "~0.7.77"
RuntimeGeneratedFunctions = "~0.5.16"
SymbolicUtils = "~4.10.0"
Symbolics = "~7.4.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "9974ef1deae23559a9ec33e742c28117df9a445d"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "5017d6849aff775febd36049f7d926a5fb6677ec"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.8"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkCodecCore]]
git-tree-sha1 = "1a3ad7e16a321667698a19e77362b35a1e94c544"
uuid = "0b6fb165-00bc-4d37-ab8b-79f91016dbe1"
version = "1.0.1"

[[deps.ChunkCodecLibZlib]]
deps = ["ChunkCodecCore", "Zlib_jll"]
git-tree-sha1 = "cee8104904c53d39eb94fd06cbe60cb5acde7177"
uuid = "4c0bbee4-addc-4d73-81a0-b6caacae83c8"
version = "1.0.0"

[[deps.ChunkCodecLibZstd]]
deps = ["ChunkCodecCore", "Zstd_jll"]
git-tree-sha1 = "34d9873079e4cb3d0c62926a225136824677073f"
uuid = "55437552-ac27-4d47-9aa3-63184e8fd398"
version = "1.0.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "76dab592fa553e378f9dd8adea16fe2591aa3daa"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.6"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"
weakdeps = ["Makie", "Random"]

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DrWatson]]
deps = ["Dates", "FileIO", "JLD2", "LibGit2", "MacroTools", "Pkg", "Random", "Requires", "Scratch", "UnPack"]
git-tree-sha1 = "5b6632df14cf24fc2cdb805aab24147001463336"
uuid = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
version = "2.19.1"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "d60eb76f37d7e5a40cc2e7c36974d864b82dc802"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.1"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bfcd42851cf2f1b303f51525a54dc5e98d408a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.15.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "6b4d2dc81736fe3980ff0e8879a9fc7c33c44ddf"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.2+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "02b61501dbe6da3b927cc25dacd7ce32390ee970"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.2"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["ChunkCodecLibZlib", "ChunkCodecLibZstd", "FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues"]
git-tree-sha1 = "8f8ff711442d1f4cfc0d86133e7ee03d62ec9b98"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.6.3"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "5b6bb73f555bc753a6153deec3717b8904f5551c"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.3.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "d1b974f376c24dad02c873e951c5cd4e351cd7c2"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.8"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f19301ae653233bc88b1810ae908194f07f8db9d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39a11854f0cba27aa41efaedf43c77c5daa6be51"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.0+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "6ed167db158c7c1031abf3bd67f8e689c8bdf2b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.77"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.ReadOnlyArrays]]
git-tree-sha1 = "e6f7ddf48cf141cb312b078ca21cb2d29d0dc11d"
uuid = "988b38a3-91fc-5605-94a2-ee2116b3bd83"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "be5733d4a2b03341bdcab91cea6caa7e31ced14b"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.9"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "79529b493a44927dd5b13dde1c7ce957c2d049e4"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.0"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils", "TermInterface"]
git-tree-sha1 = "49201c2793ce02f141c6f8b5194ce34e8012cd68"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.4"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "EnumX", "ExproniconLite", "LinearAlgebra", "MacroTools", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "ReadOnlyArrays", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "WeakCacheSets"]
git-tree-sha1 = "5e0bf32048585519709211f22cf7f375b86bb3a8"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "4.10.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsChainRulesCoreExt = "ChainRulesCore"
    SymbolicUtilsDistributionsExt = "Distributions"
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "AbstractPlutoDingetjes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "Moshi", "MultivariatePolynomials", "MutableArithmetics", "NaNMath", "PrecompileTools", "Preferences", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "ec82ed16711fde5790fdf284b6eb0a4c7e1aae69"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "7.4.1"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsDistributionsExt = "Distributions"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLatexifyExt = ["Latexify", "LaTeXStrings"]
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c25751629f5baaa27fef307f96536db62e1d754e"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.27.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WeakCacheSets]]
git-tree-sha1 = "386050ae4353310d8ff9c228f83b1affca2f7f38"
uuid = "d30d5f5c-d141-4870-aa07-aabb0f5fe7d5"
version = "0.1.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "de8ab4f01cb2d8b41702bab9eaad9e8b7d352f73"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.53+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╠═c365fb96-df3d-11f0-305f-cd90975c054e
# ╠═26b93993-390f-4278-b7ee-67350d6970be
# ╟─6e2e12ab-3a7c-411e-b74a-2f6967162fa4
# ╟─cb581f1d-ab4f-4e28-a4eb-747b77315ad6
# ╟─dc12b2fc-26d7-4a3b-8c72-ab24a1fe4d41
# ╟─eb07fe92-8237-4c5f-ac9c-8ebeb30ff1ab
# ╟─d906b120-1595-4bb9-b40a-6e17900e7539
# ╟─7b13f90d-7a51-45bd-acf2-b02545e05bde
# ╟─eeeaca1b-7861-4e66-9d48-62430f9e49c4
# ╟─34ec4989-0add-4f87-b679-393a286077dd
# ╟─0c4d6003-d3b1-4fdf-acc4-b5348b606825
# ╟─6054481d-0bba-40c7-baef-c72cfe02e632
# ╟─80e53723-b027-49ec-a178-aefc6d563672
# ╠═11dbe7db-a8ac-4286-8b43-e352f857842f
# ╟─44db3ae2-fb94-48ad-a306-9a95627618e7
# ╟─fe40d664-3137-47ea-90a0-03bf4058be5c
# ╟─4589f32c-f6dc-4805-94fc-424655fa5289
# ╠═4dcd3ed3-6fa5-4f54-a9fa-2f109095292f
# ╟─97c50791-1c10-452b-8f7c-a2e70f71872b
# ╟─5136efea-d24b-4016-8942-581a91552796
# ╟─5b815179-5799-4544-81f7-37ff80bb1442
# ╟─f7cec82a-4606-4463-8ce3-25db7683df02
# ╟─41f3ff09-cd17-4691-bd7a-dce52f1a7e97
# ╠═c1570e2e-e7c7-4cce-80a4-bd0c7efaa0c2
# ╟─c8446c65-9992-40e9-be38-7fb6ec9b0f62
# ╟─55e52182-cbc1-4af3-9141-c02580e77195
# ╟─456e8985-3197-47c6-8658-e7f7bb8253e5
# ╟─b2ff6746-af69-4a09-8ef5-46dcc1776c05
# ╟─33004562-7af3-4fb8-98d2-f2aed9d1203e
# ╟─3877ac7e-142d-4c3b-8b11-3574cd4b2f9c
# ╟─a0341bc4-7240-4d8c-a6c4-9929241851de
# ╟─31a8215d-64a0-44af-b54c-97d4de167ccd
# ╟─9b901c1a-b28f-4cd8-a800-d5d8b38b478e
# ╟─d44d24eb-e829-42a1-9372-1f54e2a56d0f
# ╟─854b6598-20b2-4aa2-9410-112cf703d96e
# ╟─8609a240-a029-48ea-b29d-8e0587b1338c
# ╟─58ef8cec-6e7a-4491-aa00-e4516094925e
# ╟─1607ff67-b476-400a-b013-0a5a3111de7a
# ╟─159467de-d9aa-4718-a8e1-012a90fa93a1
# ╟─37a62d52-1154-4820-80fc-eb64439bec68
# ╟─a40fb713-9ed6-4ec1-8f16-928d2b9a65cd
# ╟─6f28619c-df35-4942-a66b-76a7d5ee90da
# ╟─7de0eff1-04f9-4b11-a201-def19335c853
# ╟─82b95f73-f209-4302-ac27-896d463393aa
# ╟─51b7e139-b647-444e-bd4b-534ff7a28e78
# ╟─1e33e3aa-9338-472d-b6c3-45b94a304f21
# ╟─177d8297-0064-4b4f-a04e-cddf3cc30786
# ╟─70a0b93e-a878-4d6a-9a42-e24ed96abe97
# ╟─772c653a-5ca8-4672-ae48-2086fff8c639
# ╟─bc57b52f-deba-48d4-a3ae-589c19748249
# ╟─007d3be0-7cc0-4aeb-b46f-6f894a4480ef
# ╟─655ba582-8c5e-44eb-b711-6ec8369326a5
# ╟─0b6025ec-0e21-4a6a-a604-6aa98eef2350
# ╟─d75cbf05-6c02-4622-93d5-9231231cdb3c
# ╠═38313589-8fda-4c9e-bcf5-7ce2b24360a4
# ╟─7c38c4dc-2711-48f7-a31f-14b972b0e447
# ╟─1540c9d9-25ac-48d6-9c2c-fe7650d1b872
# ╟─dfc6fcc3-7357-4e5e-9061-f2571934d2d2
# ╟─92b811da-ccc3-4781-a085-5c556ce889ee
# ╟─689e54d9-2ecd-46cb-9082-efb76a0ad99d
# ╠═daa9e2bf-a1eb-4628-a974-024c7a0d61a6
# ╟─7891b7da-ba36-417d-bdde-c84187463636
# ╟─1de0d073-bc11-4e4f-832c-f3b0e01b58f7
# ╟─2474940b-0b25-44fe-995e-7e6e1dc94536
# ╠═e8dafe17-c065-426c-a928-dcec43e962b0
# ╟─2ea23935-697b-46be-886c-2f166475a590
# ╟─e45490c6-59ea-4538-9aeb-25e1398a0e23
# ╠═89277803-a204-4371-b475-8a8ff7233d60
# ╟─0f4ecc8e-d314-4ed4-bc87-bc28b0dba078
# ╠═7271061e-092f-4454-b90b-7d0d9295f815
# ╟─7499fa65-cb31-45d1-aa43-d248dc1cd066
# ╠═5fde0451-f79c-472e-a431-52b72eb2b703
# ╟─02be2a80-0fa3-42bd-ba7b-e0a9ec893c3c
# ╠═13394a42-764d-44d4-a730-66c5109647ef
# ╠═b05f634b-56a8-49c7-8a2b-d410d72582fe
# ╠═22ad33e5-0641-40aa-9aec-fc4e9be6b4b8
# ╠═b240a2e2-f574-4275-b471-2313428b18eb
# ╟─9a3027be-a9d9-47b6-81c3-2f2ea23776b7
# ╟─4af8a93f-8a2a-4531-8c54-495bf97b1f4f
# ╠═b78217cf-652a-403c-bb65-3721f963b95a
# ╟─84865a4a-9752-44e4-9746-41674be387c7
# ╟─2c658919-9be3-43f5-934b-9dd11a415eec
# ╟─8e1bcfe0-e9c8-40b2-a671-f0795d4c4b52
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
