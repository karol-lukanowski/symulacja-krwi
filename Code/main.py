#import analysis as An
import blood_oxygen as Bo
import decrease as Dc
import draw_net as Dr
import oxygen as Ox
import pressure as Pr
#import pruning as Prun
import save as Sv
import upstream as Up
import vegf as Ve

import numpy as np

from build import build
from utils import solve_equation, simAnalysisData, collect_data

from config import simInputData

### Loading config data
sid = simInputData()

### Building the network
sid, nxGraph, edges, arterioleVenuleMarker, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = build(sid)

### Object collecting data for analysis
simAnalysisData = simAnalysisData()

### Definition of the RHS of evolution equations (TO DO: MOVE OUT OF MAIN)
pressureFlowEqRHS = Pr.create_pressureFlowEqRHS(sid, in_nodes, out_nodes)           # pressure-flow eq. RHS

if sid.includeOxygenEvolution:                                                      # if the simulation is to include oxygen-related evolution mechanisms:
    arterioleOxygenEqRHS = np.zeros(sid.networkSizeSquared)                         #   create the RHS of the arteriole oxygen convection-reaction eq.
    arterioleOxygenEqRHS[in_nodes_ox] = 1                                           #   corr. to boundary conditions of unit oxygen concentration at in_nodes

iterations = sid.old_iters + sid.iterations

####################################################################################################################
# Main loop
####################################################################################################################
for iter in range(sid.old_iters, iterations):
    print(f'Iteration {iter + 1}/{iterations}')

    ####################################################################################################################
    ### a. Pressure and flow (Hagen-Poiseuille + continuity eq.)
    ####################################################################################################################
    pressureFlowEqMatrix = Pr.create_pressureFlowEqMatrix(sid, edges, in_nodes, out_nodes)
    pressureVector = solve_equation(pressureFlowEqMatrix, pressureFlowEqRHS)

    ####################################################################################################################
    ### b. Oxygen concentration in arterioles (convection-reaction eq.)
    ### c. Oxygen concentration in capillaries (diffusion-reaction eq.)
    ####################################################################################################################
    if sid.includeOxygenEvolution:
        arterioleOxygenEqMatrix = Bo.update_matrix(sid, pressureVector, arterioleVenuleMarker, arterioleOxygenEqRHS, edges)
        arterioleOxygenVector = solve_equation(arterioleOxygenEqMatrix, arterioleOxygenEqRHS)

        capillaryOxygenEqMatrix = Ox.update_matrix(sid, arterioleVenuleMarker, edges)
        capillaryOxygenVector = solve_equation(capillaryOxygenEqMatrix, arterioleOxygenVector)
        VEGFEqRHS = Ve.create_vector(sid, arterioleVenuleMarker, capillaryOxygenVector)
    else:
        VEGFEqRHS = Ve.create_vector(sid, arterioleVenuleMarker)

    ####################################################################################################################
    ### d. VEGF concentration (diffusion-reaction eq.)
    ####################################################################################################################
    VEGFEqMatrix = Ve.update_matrix(sid, VEGFEqRHS, edges)
    VEGFVector = solve_equation(VEGFEqMatrix, VEGFEqRHS)

    ####################################################################################################################
    ### e. Metabolic signal concentration
    ####################################################################################################################
    if sid.includeSignalEvolution:
        signalEqMatrix, signalEqRHS = Up.update_matrix_upstream(sid, VEGFVector, pressureVector, edges, in_nodes)
        signalVectorUpstream = solve_equation(signalEqMatrix, signalEqRHS)

        signalEqMatrix, signalEqRHS = Up.update_matrix_downstream(sid, VEGFVector, pressureVector, edges, out_nodes)
        signalVectorDownstream = solve_equation(signalEqMatrix, signalEqRHS)

        signalVector = signalVectorUpstream + signalVectorDownstream

    ####################################################################################################################
    # Plotting
    ####################################################################################################################
    if iter % sid.plot_every == 0:
        vnow2 = VEGFVector.copy()
        for node in out_nodes:
            vnow2[node] = 0

        nxGraph = Pr.update_network(nxGraph, sid, edges, pressureVector)

        #Dr.uniform_hist(sid, nxGraph, in_nodes, out_nodes, boundary_edges, oxresult, pressureVector, vnow, snow_upstream, snow, name=f'histogram{sid.old_iters // sid.plot_every:04d}.png')
        #Dr.draw(sid, nxGraph, in_nodes, out_nodes, boundary_edges, oxresult, name=f'd{sid.old_iters // sid.save_every:04d}.png', data='d')
#        Dr.drawq(name = f'q{(i+old_iters)//save_every:04d}.png', oxdraw = [])
#        Dr.drawq(name=f'veq{i // save_every:04d}.png', oxdraw=vnow2)
#        Dr.drawq(name=f'oxq{i // save_every:04d}.png', oxdraw=snow)
#        Dr.drawblood(sid, nxGraph, in_nodes, out_nodes, boundary_edges, name=f'q_blood{sid.old_iters // sid.save_every:04d}.png', oxresult=oxresult, oxdraw = vnow, data='q')
        Dr.drawvessels(sid, nxGraph, in_nodes, out_nodes, boundary_edges, name=f'q_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=arterioleVenuleMarker, oxdraw = capillaryOxygenVector, data='q')
        Dr.drawvessels(sid, nxGraph, in_nodes, out_nodes, boundary_edges, name=f'd_vessels{sid.old_iters // sid.plot_every:04d}.png', oxresult=arterioleVenuleMarker, oxdraw = signalVector, data='d')
    
    if sid.shear_d:
        edges = Pr.update_graph(sid, edges, pressureVector)
    if sid.vegf_d:
        edges = Ve.update_graph(sid, VEGFVector, arterioleVenuleMarker, edges)
    if sid.signal_d:
        edges = Up.update_graph_upstream(sid, signalVectorUpstream, pressureVector, arterioleVenuleMarker, edges)
        edges = Up.update_graph_downstream(sid, signalVectorDownstream, pressureVector, arterioleVenuleMarker, edges)
    if sid.gradp_d:
        edges = Pr.update_graph_gradp(sid, edges, pressureVector)
        #edges = Pr.update_graph_gradp_tips(sid, edges, pressureVector, oxresult)
    if sid.decrease_d:
        edges = Dc.update_graph(sid, edges)
    arterioleVenuleMarker = Ve.update_blood(sid, arterioleVenuleMarker, edges)

    if sid.data_collection:
        collect_data(simAnalysisData, sid, in_nodes, pressureVector, VEGFVector, pressureVector)

    if iter % sid.save_every == 0 and iter != 0:
        Sv.save(f'/save{sid.old_iters}.dill', sid, nxGraph, edges, arterioleVenuleMarker, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)

    sid.old_iters += 1

Sv.save('/save.dill', sid, nxGraph, edges, arterioleVenuleMarker, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
#Dr.plot_params(sid)

#An.getStrahlerHistogram(nxGraph, pressureVector, oxresult, in_nodes, dirname)


#DiG = An.getDiGraph(nxGraph, pressureVector, oxresult, in_nodes)
#DiG = An.strahlerOrder(DiG)
#An.plotStrahlerGraph(DiG, 'deown101deown_drabina/26')

#Prun.pruning(nxGraph, reg_reg_edges, reg_something_edges, in_edges, pressureVector, pressureFlowEqRHS, oxresult)