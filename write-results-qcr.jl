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

paths = ["../../04_Data/QMSTPInstances/" * benchmarkclass * "/" for benchmarkclass in ["test_small", "SV", "CP_all"]]
paths = paths[1:1]

for path in paths
    logfilespath = path * "logs-qmstpqcr/"
    jsonfiles = filter(x->endswith(x, ".json"), readdir(logfilespath,sort=false))
    io = open(logfilespath * "summary.csv", "w")
    write(io, "instance,n,m,objval,objbound,gap,time,timegurobi,bbnodes,ncuts,ncutset,nsubtour,rootbound,timesdp,timelimit\n")
    for instance in jsonfiles
        result = JSON.parsefile(logfilespath * instance)
        instancename = result["instance"]
        n = result["n"]; m = result["m"];
        timesdp = result["time-sdp"]
        rootbound = result["root-bound"]
        timegurobi = result["time-gurobi"]
        time = timesdp + timegurobi
        timelimit = result["timelimit"]
        objval = result["objective_value"]
        objbd = result["objective-bound"]
        gap = result["gap"]
        bbnodes = result["bbnodes-gurobi"]
        ncutset = result["ncutsetconstraints"]
        nsubtour = result["nsubtourconstraints"]
        ncuts = ncutset + nsubtour

        write(io, "$instancename,$n,$m,$objval,$objbd,$gap,$time,$timegurobi,$bbnodes,$ncuts,$ncutset,$nsubtour,$rootbound,$timesdp,$timelimit\n")
    end
    close(io)
end
