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

include("data-io.jl") # functions to read files and extract cost matrix Q
include("QMSTP-QCR.jl") # functions to compute the QCR and solve integer formulation with lazy cuts

filepath = "../../04_Data/QMSTPInstances/CP1/qmstp_CP15_67_10_10.dat"
n, m, edges, Q = readInput_qmstp(filepath)
timelimit = 100
# adds subtour and cutset constraints (+ GC cuts we first derived from LMI if cutgc1 = true)
res= solve_qmstp_conv(n, Q, edges, timelimit, cutgc1=false) 
