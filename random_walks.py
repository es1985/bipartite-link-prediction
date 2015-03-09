import util
import datetime
import networkx as nx
import numpy as np
from scipy import sparse
from dataset_maker import get_date


def run_random_walks(data_dir, weight_edges=False):
    print "Loading data and building transition matrix..."
    data_dir = 'train'
    weight_edges=False
    examples = util.load_json('./data/' + data_dir + '/examples.json')
    G = nx.read_edgelist('./data/' + data_dir + '/graph.txt', nodetype=int)
    if weight_edges:
        reviews = util.load_json('./data/' + data_dir + '/review.json')
        end_date = datetime.date(2012, 1, 1) if data_dir == 'train' else datetime.date(2013, 1, 1)
        edges = G.edges()
        for e in util.logged_loop(edges, util.LoopLogger(20000, len(edges), True)):
            n1, n2 = str(e[0]), str(e[1])
            if n1 not in reviews or n2 not in reviews[n1]:
                n1, n2 = n2, n1
            G[e[0]][e[1]]['weight'] = 1.0 / ((end_date - get_date(reviews[n1][n2][0])).days + 90)
        del reviews  # save some memory

    adjacency_matrix = nx.adjacency_matrix(G)
    inverse_degree_matrix = sparse.diags([[1.0 / adjacency_matrix.getrow(i).sum()
                                           for i in range(adjacency_matrix.shape[0])]], [0])
    transition_matrix = inverse_degree_matrix.dot(adjacency_matrix)

    print "Running random walks..."
    for u in util.logged_loop(examples, util.LoopLogger(10, len(examples), True)):
        p = run_random_walk(transition_matrix, int(u), 10,0.2)
        p = p.todense()
        for b in examples[u]:
            examples[u][b] = p[0, int(b)]

    util.write_json(examples, './data/' + data_dir
                    + ('/weighted_random_walks.json' if weight_edges else '/random_walks.json'))

    return    
    
def run_random_walk(transition_matrix, u, iterations, jump_p):
    p = [0] * transition_matrix.shape[0]
    p[int(u)] = 1.0
    p = sparse.csr_matrix(p)
    
    for i in range(iterations):
        p = p.dot(transition_matrix)
        p *= (1 - jump_p)
        #        p[0, u] += jump_p
    
    return p


if __name__ == '__main__':
    run_random_walks('train', False)
    run_random_walks('train', False)
    run_random_walks('test', False)
    run_random_walks('test', True)
