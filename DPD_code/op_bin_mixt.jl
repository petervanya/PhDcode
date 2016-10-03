#!/usr/bin/julia
"""
Order parameter in Julia

04/09/16
"""

function gen_ordered_box(N1, N2)
  types = vcat(ones(Int, N1), 2*ones(Int, N2))
  box = Float64[L/2, L, L]'
  xyz1 = rand(N1, 3).*box
  xyz2 = rand(N2, 3).*box
  xyz2[:, 1] += L/2
  xyz = vcat(xyz1, xyz2)
  types, xyz
end


function gen_disordered_box(N1, N2)
  types = vcat(ones(Int, N1), 2*ones(Int, N2))
  xyz = rand(N1+N2, 3) * L
  types, xyz
end


function dist_vec(xyz, L)
  N = size(xyz)[1]
  Np = Int(N*(N-1)/2)
  dv = zeros(Np)
  tv = zeros(Int, Np, 2)

  cell = L * eye(3)
  inv_cell = pinv(cell)
  G, Gn = zeros(3), zeros(3)
  dr, drn = zeros(3), zeros(3)

  cnt = 1
  for i = 1:N
    for j in 1:i-1
      dr = xyz[i, :] - xyz[j, :]
      G = inv_cell * dr
      Gn = G - round(G)
      dr = cell * Gn
      dv[cnt] = norm(dr)
      tv[cnt] = [types[i], types[j]]
      cnt += 1
    end
  end
  dv, tv
end


function order_param_naive(types, xyz, L, rc)
  N = size(xyz)[1]
  cell = L * eye(3)
  inv_cell = pinv(cell)
  G, Gn = zeros(3), zeros(3)
  dr, drn = zeros(3), zeros(3)
  op = Float64[]
  d = 0.0
  phi = 0.0

  for i = 1:N
    n1, n2 = 0, 0
    for j = 1:N
      dr = (xyz[i, :] - xyz[j, :])'
      G = inv_cell * dr
      Gn = G - round(G)
      dr = cell * Gn
      d = norm(dr)
      if i == j
        continue
      end
      if d < rc
        if types[j] == 1
          n1 += 1
        elseif types[j] == 2
          n2 += 1
        end
      end
    end
    if n1 + n2 != 0
      phi = (n1 - n2)^2 / (n1 + n2)^2
      push!(op, phi)
    end
  end
  sum(op) / length(op)
end


function order_param_lc(types, xyz, lc, L, rc)
  N = size(xyz)[1]
  cell = L * eye(3)
  inv_cell = pinv(cell)
  G, Gn = zeros(3), zeros(3)
  dr, drn = zeros(3), zeros(3)
  op = Float64[]
  d = 0.0
  phi = 0.0

  for i = 1:N
    n1, n2 = 0, 0
   
    id_i = id_from_coord(div(xyz[i, :], L/Nx), Nx)
    neigh_cells = neighbour_cells(id_i, Nx)
    neigh_atoms = Int[]
    for nc in neigh_cells
      append!(neigh_atoms, lc[nc])
    end

    for j in neigh_atoms
      dr = (xyz[i, :] - xyz[j, :])'
      G = inv_cell * dr
      Gn = G - round(G)
      dr = cell * Gn
      d = norm(dr)
      if i == j
        continue
      end
      if d < rc
        if types[j] == 1
          n1 += 1
        elseif types[j] == 2
          n2 += 1
        end
      end
    end
    if n1 + n2 != 0
      phi = (n1 - n2)^2 / (n1 + n2)^2
      push!(op, phi)
    end
  end
  sum(op) / length(op)
end


# link cells
function cell_coord(id, Nx)
  if id > Nx^3
    error("ID must be less than the number of cells Nx^3.")
  end
  nx = Int(div(id-1, Nx^2))
  ny = Int(div(id-1 - nx*Nx^2, Nx))
  nz = Int(id-1 - nx * Nx^2 - ny * Nx)
  [nx ny nz]
end


function id_from_coord(n, Nx)
  n = round(Int, n)
  if all(n.>Nx)
    error("Each coordinate must be less than Nx.")
  end
  id = Int(n[1] * Nx^2 + n[2] * Nx + n[3]) + 1   # id numbering from 1
end


function atom_to_cell(r, Lx)
  n = round(Int, div(r, Lx))
  id_from_coord(n)
end


function product(arr, repeat=3)
  res = zeros(Int, length(arr)^3, 3)
  cnt = 1
  for i in arr
    for j in arr
      for k in arr
        res[cnt, :] = [i j k]
        cnt += 1
      end
    end
  end
  res
end


function neighbour_cells(id, Nx)
  r = cell_coord(id, Nx)
  nn = 3^3              # number of neighbours
  neighs = zeros(Int, nn, 3)
  p = product([-1 0 1], 3)
  for i in 1:nn
    neighs[i, :] = mod(r + p[i, :], Nx)
  end
  ids = [id_from_coord(neighs[i, :], Nx) for i = 1:nn]
end


function init_lc(Lx, Nx)
  lc = Dict([i => Int[] for i = 1:Nx^3])
end


function populate_lc(lc, xyz, Lx, Nx)
  N = size(xyz)[1]
  for i in 1:N
    id = id_from_coord(div(xyz[i, :], Lx), Nx)
    push!(lc[id], i)
  end
end


function test()
  Nx = 10
  @printf "Testing link cells Nx = %i\n" Nx
  n = [0 0 0]
  @printf "ID from coord [0 0 0] (should be 1) %i\n" id_from_coord(n, Nx)
  n = [1 1 1]
  @printf "ID from coord [1 1 1] (should be 112) %i\n" id_from_coord(n, Nx)

  id = 1000
  println("coord from ID = 1000 (should be [9 9 9]): ", cell_coord(id, Nx))

  println("Neighbour cells of ID = 1:")
  println(neighbour_cells(414, Nx))
end


# ===== main
srand(1234)
L = 10.0
Nx = 6
f = 0.5
N = 1000
rc = 1.3

if length(ARGS) == 2
  println("2 args")
  N = parse(Int, ARGS[1])
  rc = parse(Float64, ARGS[2])
elseif length(ARGS) == 1
  if ARGS[1] == "-h"
    @printf "Usage:\n    op_bin_mixt.jl <N> <rc>\n"
    exit(0)
  else
    N = parse(Int, ARGS[1])
  end
end

N1, N2 = Int(f * N), Int((1-f) * N)
@printf "N: %i | rc: %.2f | L: %.1f | f: %.2f | Nx: %i\n" N rc L f Nx

# ===== one cell
@printf "Testing ordered box.\n"
types, xyz = gen_ordered_box(N1, N2)
tic()
op = order_param_naive(types, xyz, L, rc)
toc()
@printf "=== op: %.3f\n" op

@printf "Testing disordered box.\n"
types, xyz = gen_disordered_box(N1, N2)
tic()
op = order_param_naive(types, xyz, L, rc)
toc()
@printf "=== op: %.3f\n" op

# ===== link cells
#Lx = L / Nx
#lc = init_lc(Lx, Nx)
##println(lc)
#if Lx < rc
#  println("WARNING: Cell size Lx smaller than cutoff rc.")
#end
#
#types, xyz = gen_ordered_box(N1, N2)
#populate_lc(lc, xyz, Lx, Nx)
##test()
##for i = 1:Nx^3
##  println(i, "  ", lc[i])
##end
#
#tic()
#op = order_param_lc(types, xyz, lc, L, rc)
#toc()
#@printf "=== op: %.3f\n" op
#
#
