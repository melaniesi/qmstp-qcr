# ==============================================================================
#   QMSTP-QCR -- Code to compute the quadratic minimum spanning tree of a graph.
# ------------------------------------------------------------------------------
#   Copyright (C) 2024 Melanie Siebenhofer <melaniesi@edu.aau.at>
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see https://www.gnu.org/licenses/.
# ==============================================================================

# see below for an example!

using Mosek, MosekTools, JuMP
using LinearAlgebra
using Gurobi

include("data-io.jl");

"""
    backtracking(i, j, parent, edges)

Returns a list of edges in the cycle that has been
detected in the DFS algorithm by finding another path
from `i` to its neighbor `j`, where `j` has already been
visited before in the algorithm.

In `parent` on position k is the parent/predecessor (in the depth
first traversal) of vertex k if it has already been
visited. `edges` is a list of all edges in the graph. 
"""
function backtracking(i, j, parent, edges)
    cycle = [(j,i)]
    while !((parent[i],j) in edges) && !((j,parent[i]) in edges)
        push!(cycle, (i, parent[i]))
        i = parent[i]
    end
    push!(cycle,(i, parent[i]))
    push!(cycle,(parent[i], j))
    return(cycle)
end

"""
    find_cycles(n, edges)

DFS algorithm to find all cycles in a graph
with `n` vertices and given `edges`.

Returns a list of cycles, represented by lists
of edges in the respictive cycle.
"""
function find_cycles(n, edges)
    edges_copy = copy(edges)
    parent = [-1 for i=1:n]
    visited = [false for i=1:n]
    cycles = []

    while !isempty(edges_copy)
        to_check = [edges_copy[1][1]]
        visited[to_check[1]] = true
        while !isempty(to_check)
            i = pop!(to_check)
            Ni = [j for j=1:n if (i,j) in edges_copy || (j,i) in edges_copy]
            for j in Ni
                if visited[j]
                    push!(cycles, backtracking(i, j, parent, edges))
                else
                    visited[j] = true
                    parent[j] = i
                    push!(to_check, j)
                end
            end
            edges_copy = [e for e in edges_copy if !(i in e)]
        end
    end
    return cycles
end

"""
    vertices_in_cycles(cycles)

Returns a list of subsets of vertices forming `cycles`.
The `cycles` are returned as a list of edges in
the corresponding cycles. 
"""
function vertices_in_cycles(cycles)
    subsets = []
    for c in cycles
        S = unique(collect(Iterators.flatten(c)))
        push!(subsets, S)
    end
    return subsets
end

"""
    find_connected_components(n, edges)

DFS algorithm to find all connected components of a graph
with `n` vertices and given `edges`.

Returns a list of the vertices in the components.
"""
function find_connected_components(n, edges)
    edges_copy = copy(edges)
    visited = [false for i=1:n]
    components = []

    while !isempty(edges_copy)
        to_visit = [edges_copy[1][1]]
        visited[to_visit[1]] = true
        comp = []
        while !isempty(to_visit)
            i = pop!(to_visit)
            push!(comp, i)
            Ni = [j for j=1:n if (i,j) in edges_copy || (j,i) in edges_copy]
            append!(to_visit, [j for j in Ni if !visited[j]])
            visited[Ni] .= true
            edges_copy = [e for e in edges_copy if !(i in e)]
        end
        push!(components,comp)
    end
    for i in [j for j=1:n if !visited[j]] # isolated vertices
        push!(components, [i])
    end
    return components
end

"""
    find_cycles_and_components(n, edges)

DFS algorithm to find all cycles and connected
components in a graph with `n` vertices and
given `edges`.

Returns a list of edges in the cycles and a list of
the vertices in the components.
"""
function find_cycles_and_components(n, edges)
    edges_copy = copy(edges)
    parent = [-1 for i=1:n]
    visited = [false for i=1:n]
    cycles = []
    components = []

    while !isempty(edges_copy)
        to_check = [edges_copy[1][1]]
        visited[to_check[1]] = true
        comp = []
        while !isempty(to_check)
            i = pop!(to_check)
            append!(comp, i)
            Ni = [j for j=1:n if (i,j) in edges_copy || (j,i) in edges_copy]
            for j in Ni
                if visited[j]
                    push!(cycles, backtracking(i,j,parent, edges))
                else
                    visited[j] = true
                    parent[j] = i
                    push!(to_check, j)
                end
            end
            edges_copy = [e for e in edges_copy if !(i in e)]
        end
        push!(components, comp)
    end
    for i in [j for j=1:n if !visited[j]] # isolated vertices
        push!(components, [i])
    end
    return (cycles, components)
end


"""
    solve_SDPRelax_primal_withDuals(n, Q)

Solves SDPQ from the paper "Improving the
performance of standard solvers for quadratic
0 - 1 programs by a tight convex reformulation:
QCR method" by Billionet, Elloumi, Plateau (2008)
for the BQP
    
    min y^T * Q * y
    st. e^T * y = (n - 1)
        y ∈ {0,1}^n
    
using JuMP and Mosek.

Returns a dictionary containing
    * a: the optimal values α* - dual of the constraints sum_e Y_ef = (n - 1) * y_f
    * u: the optimal values u* - dual of the constraints Y_ee = y_e
    * opt: the optimal value of the SDP (relaxation) solved
    * Y: the optimum of the SDP (Y_{m+1,m+1} = 1)
    * time: the time needed to solve the SDP
"""
function solve_SDPRelax_primal_withDuals(n, Q)
    m = size(Q, 1)
    model = Model(Mosek.Optimizer)
    Y = @variable(model, Y[1:m+1,1:m+1], PSD)
    @objective(model, Min, LinearAlgebra.dot(Q - diagm(diag(Q)), Y[1:m,1:m]) + sum( 1/2 * (Y[i,m+1] + Y[m+1,i]) * Q[i,i] for i=1:m)) # ⟨ Q - Diag(diag(Q) , Y) + y^T * diag(Q)
    r = @constraint(model, Y[m+1,m+1] == 1)
    u = @constraint(model, [e=1:m], 1/2*(Y[e,m+1]+Y[m+1,e]) - Y[e,e] == 0) # Y_ee = y_e
    b = @constraint(model, (n-1)*Y[m+1,m+1] - 1/2*(sum(Y[1:m,m+1]) + sum(Y[m+1,1:m])) == 0) # sum_e y_e = (n - 1)
    #@constraint(model, sum(Y[1:m,1:m]) == (n-1)^2) -> constraint is not in SDP for QCR
    a = @constraint(model,[i=1:m], (n-1)/2*(Y[i,m+1]+Y[m+1,i]) - 1/2*(sum(Y[i,1:m]) + sum(Y[1:m,i])) == 0) # sum_e Y_ef = (n - 1) y_f
    optimize!(model)
    return Dict("a" => dual.(a), "u" => dual.(u), "opt" => objective_value(model), "Y" => value.(Y), "time" => solve_time(model))
end


"""
    convexify_qmstp(n::Int, Q::Matrix)

Applies the quadratic convex reformulation (QCR)
method of Billionet, Elloumi, Plateau (2008) to
the BQP

    min y^T * Q * y
    st. e^T * y = (n - 1)
        y ∈ {0,1}^n.

Returns Q_conv, c_conv and the time needed to
solve the SDP.
The reformulated convex QP is then
    min y^T * Q_conv * y + c_conv^T * y
    st. e^T * y = (n - 1)
        y ∈ {0,1}^n
and the continuous relaxation of this CQP is
equal to the SDP relaxation bound.
"""
function convexify_qmstp(n::Int, Q::Matrix)
    result = solve_SDPRelax_primal_withDuals(n, Q)
    m = size(Q, 1)
    u = result["u"]
    a = result["a"]
    c_conv = diag(Q)
    Q_conv = Q - diagm(c_conv - u) + 1 / 2 * ( repeat(a, 1, m) + repeat(a', m, 1) )
    c_conv -= (n - 1) * a + u
    return Q_conv, c_conv, result["time"]
end


"""
    convexify_qmstp(result_sdp::Dict, n::Int, Q::Matrix)

Uses the result `result_sdp` of
solve_SDPRelax_primal_withDuals(n, Q) to construct
the matrices Q_conv and c_conv for the CQR.
The reformulated convex QP is then
    min y^T * Q_conv * y + c_conv^T * y
    st. e^T * y = (n - 1)
        y ∈ {0,1}^n
Further, the time needed to solve the SDP is
returned as well.
"""
function convexify_qmstp(result_sdp::Dict, n::Int, Q::Matrix)
    m = size(Q, 1)
    u = result_sdp["u"]
    a = result_sdp["a"]
    c_conv = diag(Q)
    Q_conv = Q - diagm(c_conv - u) + 1 / 2 * ( repeat(a, 1, m) + repeat(a', m, 1) )
    c_conv -= (n - 1) * a + u
    return Q_conv, c_conv, result_sdp["time"]
end

"""
    solve_SDPRelax_dual(n, Q)

Solves the dual of SDPQ from the paper "Improving the
performance of standard solvers for quadratic
0 - 1 programs by a tight convex reformulation:
QCR method" by Billionet, Elloumi, Plateau (2008)
for the BQP
    
    min y^T * Q * y
    st. e^T * y = (n - 1)
        y ∈ {0,1}^n
    
using JuMP and Mosek.

Returns a dictionary containing
    * a: the optimal values α* - dual of the constraints sum_e Y_ef = (n - 1) * y_f
    * u: the optimal values u* - dual of the constraints Y_ee = y_e
    * opt: the optimal value of the SDP (relaxation) solved
    * Y: the optimum of the SDP (Y_{m+1,m+1} = 1)
    * time: the time needed to solve the SDP
"""
function solve_SDPRelax_dual(n::Int, Q::Matrix)
    m = size(Q, 1)
    dsdp = Model(Mosek.Optimizer)
    r = @variable(dsdp, r)
    a = @variable(dsdp, a[1:m])
    u = @variable(dsdp, u[1:m])
    b = @variable(dsdp, b)

    @objective(dsdp, Max, r)

    Qmat = Q - diagm( diag(Q) ) + diagm(u) + 1 / 2 * ( repeat(a, 1, m) + repeat(a', m, 1) )
    cvec = diag(Q) - (n - 1) * a - u + b * ones(m)
    mat = Matrix{AffExpr}(undef, m + 1, m + 1)
    mat[1:m,1:m] .= Qmat
    mat[m+1,m+1] = - (n - 1) * b - r
    mat[1:m,m+1] .= 1 / 2 * cvec
    mat[m+1,1:m] .= 1 / 2 * cvec
    @constraint(dsdp, mat in PSDCone())

    optimize!(dsdp)

    return Dict("a" => value.(a), "u" => value.(u), "opt" => value.(r), "time" => solve_time(dsdp))
end


"""
    add_lazycon_subtour(model, edges_graph, X)

Adds the subtour elimination constraints to `model`
as lazy constraints for integer points.
For each detected cycle, a constraint gets added.
In `edges` a list of edges in the graph is provided and
`X` is the variable modelling the adjacency of the subgraph.
"""
function add_lazycon_subtour(model, edges_graph, X)
    function callback_function_subtourelim(cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            x_val = callback_value.(cb_data, model[:x])
            #x = model[:x]
            n = size(X, 1)
            edges = edges_graph[findall(e->abs(e-1)<1e-6, x_val)]
            cycles = find_cycles(n, edges)
            if !isempty(cycles)
                for S in vertices_in_cycles(cycles)
                    #println("add subtour constraint for cycle $S")
                    print("-")
                    subtour_con = @build_constraint(sum(X[S,S]) <= 2*(length(S)-1))
                    #println(subtour_con)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), subtour_con)
                end
                print("\n")
            else
                return
            end
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), callback_function_subtourelim)
end

"""
    add_lazycon_subtour(model, edges_graph, X)

Adds the subtour elimination and cutset constraints to `model`
as lazy constraints for integer points.

In `edges` a list of edges in the graph is provided and
`X` is the variable modelling the adjacency of the subgraph.
If cutgc1 = true the first type of GC-cuts derived from the LMI are added on top.
"""
function add_lazycon_subtour_cutset(model, edges_graph, X, countercuts; cutgc1=false)
    function callback_function_subtourcutset(cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            x_val = callback_value.(cb_data, model[:x])
            #x = model[:x]
            n = size(X, 1)
            edges = edges_graph[findall(e->abs(e-1)<1e-6, x_val)]
            cycles, components = find_cycles_and_components(n, edges)
            nsubtourcon = 0; ncutsetcon = 0; ncutgc1 = 0;
            if !isempty(cycles)
                for S in vertices_in_cycles(cycles)
                    #println("add subtour constraint for cycle $S")
                    subtour_con = @build_constraint(sum(X[S,S]) <= 2*(length(S)-1))
                    #println(subtour_con)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), subtour_con)
                    nsubtourcon += 1
                end
                beta = 2*(1-cos(pi/n))
                comp = components[1]
                notInComp = [i for i=1:n if !(i in comp)]
                #println("add cut set constraint for S = $(components[1])")
                cutset_con = @build_constraint(sum(X[comp, notInComp]) >= 1)
                MOI.submit(model, MOI.LazyConstraint(cb_data), cutset_con)
                ncutsetcon += 1
                if cutgc1
                    k = length(comp)  
                    v = -k*ones(n)
                    v[comp] .+= n
                    cutset_con = @build_constraint(LinearAlgebra.dot(v*v' - v.^2*ones(n)', X) <= floor(-beta*k*n*(n-k)))
                    MOI.submit(model, MOI.LazyConstraint(cb_data), cutset_con)
                    ncutgc1 += 1
                end
                if length(components) > 2
                    for comp in components[2:end]
                        notInComp = [i for i=1:n if !(i in comp)]
                        #println("add cut set constraint for S = $(comp)")
                        cutset_con = @build_constraint(sum(X[comp, notInComp]) >= 1)
                        MOI.submit(model, MOI.LazyConstraint(cb_data), cutset_con)
                        ncutsetcon += 1
                        if cutgc1
                            k = length(comp)  
                            v = -k*ones(n)
                            v[comp] .+= n
                            cutset_con = @build_constraint(LinearAlgebra.dot(v*v' - v.^2*ones(n)', X) <= floor(-beta*k*n*(n-k)))
                            MOI.submit(model, MOI.LazyConstraint(cb_data), cutset_con)
                            ncutgc1 += 1
                        end
                    end
                    print("\n")
                end
                println("Added $ncutsetcon cut-set constraints, $nsubtourcon subtour-elimination constraints, $ncutgc1 GC cuts of type 1.")
                countercuts[1] += ncutsetcon
                countercuts[2] += nsubtourcon
                countercuts[3] += ncutgc1
            else
                return
            end
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), callback_function_subtourcutset)
end

"""
    solve_qmstp_conv(n, Q, edges; cutgc1=false)

"""
function solve_qmstp_conv(n, Q, edges, timelimit=3600.0; cutgc1=false)
    @assert Q == Q' "Q is not symmetric!"
    m = length(edges)
    println("Convexify objective function, compute SDP relaxation.")
    result_sdp = solve_SDPRelax_primal_withDuals(n, Q) # compute CQR
    Q_conv, c_conv, time_sdp = convexify_qmstp(result_sdp, n, Q) # build cost matrix and vector
    delta = eigmin(Q_conv) # ensure Q_conv is psd
    Q_conv -= delta * I    # subtract vom Q_conv to make PSD
    c_conv += delta * ones(m) # add to c_conv

    println("Solve QMSTP with MIQP solver Gurobi.")
    qmodel = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(qmodel, "TimeLimit", timelimit)
    set_optimizer_attribute(qmodel, "LazyConstraints", 1)
    set_optimizer_attribute(qmodel, "MIPGapAbs", 1)
    x = @variable(qmodel, x[1:m], Bin) # vector of edge variables
    # X ... adjacency matrix built from edge vars
    #       (used for adding lazy constraints only)
    X = Matrix{Union{VariableRef, Int}}(fill(0, n, n)) 
    for i=1:m # this is faster than value.(X)
        X[CartesianIndex(edges[i])] = x[i]
        X[CartesianIndex(reverse(edges[i]))] = x[i]
    end
    @objective(qmodel, Min, x'*Q_conv*x + c_conv'*x);
    @constraint(qmodel, sum(x) == n - 1);
    countercuts = [0,0,0]
    add_lazycon_subtour_cutset(qmodel, edges, X, countercuts, cutgc1=cutgc1)

    optimize!(qmodel)

    res = Dict{String,Any}("xopt" => Int.(round.(value.(x))), # ensure datatype is integer
                           "objective-bound" => objective_bound(qmodel),
                           "time-sdp" => time_sdp,
                           "root-bound" => result_sdp["opt"],
                           "timelimit" => timelimit,
                           "time-gurobi" => solve_time(qmodel),
                           "bbnodes-gurobi" => node_count(qmodel),
                           "ncutsetconstraints" => countercuts[1],
                           "nsubtourconstraints" => countercuts[2],
                           "nconstraints-gctype1" => countercuts[3])
    res["objective_value"] = termination_status(qmodel) == OPTIMAL ? objective_value(qmodel) : "-";
    res["gap"] = termination_status(qmodel) == TIME_LIMIT ? relative_gap(qmodel) : "-";
        
    return res
end

# data from:
# n, m, edges, Q = readInput_qmstp("../../04_Data/QMSTPInstances/CP1/qmstp_CP10_33_10_10.dat");
n = 10; m = 14;
edges = [(1, 8), (1, 10), (2, 4), (2, 6), (2, 10), (3, 10), (4, 5),
         (4, 6), (4, 10), (6, 7), (6, 10), (7, 10), (8, 10), (9, 10)];
Q = [7 6 7 9 1 4 6 5 5 8 7 3 2 4;
     6 3 9 4 3 8 4 1 3 3 3 4 6 7; 
     7 9 6 3 5 6 9 7 2 5 8 5 5 8; 
     9 4 3 6 1 5 8 8 3 5 6 8 2 2; 
     1 3 5 1 8 5 5 4 5 1 2 6 4 5; 
     4 8 6 5 5 1 7 3 7 5 8 8 2 3; 
     6 4 9 8 5 7 9 3 1 2 9 6 6 2; 
     5 1 7 8 4 3 3 2 9 9 5 3 8 6; 
     5 3 2 3 5 7 1 9 9 5 2 7 3 8; 
     8 3 5 5 1 5 2 9 5 1 7 8 8 5; 
     7 3 8 6 2 8 9 5 2 7 7 3 1 5; 
     3 4 5 8 6 8 6 3 7 8 3 9 5 2; 
     2 6 5 2 4 2 6 8 3 8 1 5 2 2; 
     4 7 8 2 5 3 2 6 8 5 5 2 2 1]
timelimit = 100
res = solve_qmstp_conv(n, Q, edges, timelimit; cutgc1=true)
