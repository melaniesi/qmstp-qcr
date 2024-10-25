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

using JSON
include("QMSTP-QCR.jl")

paths = ["../../04_Data/QMSTPInstances/" * benchmarkclass * "/" for benchmarkclass in ["test_small", "SV", "CP_all"]]

function evaluate_qmstp_qcr(filepath, timelimit=3600.0; writelogfile=true, movetoprocessed=true)
    n, m, edges, Q = readInput_qmstp(filepath)
    result = solve_qmstp_conv(n, Q, edges, timelimit; cutgc1=false)
    if writelogfile || movetoprocessed
        filename = split(filepath, "/")[end]
        instancename = split(filename, ".", keepempty=false)[begin]
        pathdir = split(filepath, instancename)[begin]
        if writelogfile
            result["instance"] = instancename
            result["n"] = n; result["m"] = m;
            if !isdir(pathdir*"logs-qmstpqcr") mkdir(pathdir*"logs-qmstpqcr"); end
            io = open(pathdir * "logs-qmstpqcr/" * instancename* ".json", "w")
            JSON.print(io, result, 4)
            close(io)
        end
        if movetoprocessed
            if !isdir(pathdir*"processed") mkdir(pathdir*"processed"); end
            mv(filepath, pathdir*"processed/"*filename, force=true)
        end
    end
    return result
end

# sort graph files in directory in ascending order of the number of edges
function get_sortperm_graphFiles(graphFiles, path_dir, threshold=100000)
    nedges = similar(graphFiles, Int64)
    for (i, filename) in enumerate(graphFiles)
        filepath = path_dir * filename
        @show filepath
        _, m, _, _ = readInput_qmstp(filepath, threshold)
        nedges[i] = m
    end
    return sortperm(nedges)
end

function evaluate_allindir(pathdir)
    graphFiles = filter(x->(endswith(x,r".dat|txt")), readdir(pathdir, sort=false))
    perm = get_sortperm_graphFiles(graphFiles, pathdir)
    for graphfile in graphFiles[perm]
        evaluate_qmstp_qcr(pathdir*graphfile, writelogfile=true, movetoprocessed=true)
    end
end

function evaluate_all(paths)
    for pathdir in paths
        if !(pathdir[end] == "/") pathdir = pathdir * "/" end
        evaluate_allindir(pathdir)
    end
end

evaluate_all(paths0)

