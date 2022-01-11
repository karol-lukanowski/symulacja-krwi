import delaunay as De
import save as Sv

from geometry import set_geometry, create_edgelist
from oxygen import create_vector
from utils import make_dir

from config import simInputData


def build(sid:simInputData):

    ### Starting a new simulation with data from config
    if sid.loadMode == 'new_simulation':
        nxGraph, boundary_edges = De.Build_delaunay_net(sid.networkSize, periodic = sid.periodic, diameterStats = sid.diameterStats)
        in_nodes, out_nodes, reg_nodes, boundary_nodes_out, boundary_nodes_in = \
            set_geometry(sid.networkSize, nxGraph, geo=sid.geo, R=sid.networkSize//2.5, R_s=sid.networkSize//20, in_nodes=sid.in_nodes_own, out_nodes=sid.out_nodes_own)

        edges = create_edgelist(nxGraph, in_nodes, out_nodes, reg_nodes, boundary_nodes_out, boundary_nodes_in)
        
        in_nodes_ox = in_nodes
        out_nodes_ox = out_nodes
        make_dir(sid)
        oxresult = create_vector(sid, in_nodes_ox, out_nodes_ox)

        Sv.save('/template.dill', sid, nxGraph, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges)
        Sv.save_config(sid)

    ### Loading an existing network to evolve it further
    elif sid.loadMode == 'load_network':
        sid, nxGraph, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = Sv.load(sid.loadDirectory+'/save.dill')

    ### Loading a template network
    elif sid.loadMode == 'load_template':
        sid2, nxGraph, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges = Sv.load(sid.loadDirectory+'/template.dill')
        sid.networkSize = sid2.networkSize
        sid.dirname = sid2.dirname + '/template'
        make_dir(sid)
        Sv.save_config(sid)

    return sid, nxGraph, edges, oxresult, in_nodes, out_nodes, in_nodes_ox, out_nodes_ox, boundary_edges