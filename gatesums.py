from simplify import *

def sum_gates_perm(n,d,top):
    res = 0
    for g in ['d','m']:
        res+= perm(n,d,g,top,False)
    return res 

def sum_idle_perm(n,d,top):
    res = 0
    for g in ['id','im','ic']:
        res+= perm(n,d,g,top,False)
    return res 

def sum_gates_FO(n,top):
    res = 0
    for g in ['d','m']:
        res+= fo(n,g,top,False)
    return res

def sum_idle_FO(n,top):
    res = 0
    for g in ['id','im','ic']:
        res+= fo(n,g,top,False)
    return res


def sum_gates_OR(n,top):
    res = 0
    for g in ['d','m']:
        res+= orr(n,g,top,False)
    return res

def sum_idle_OR(n,top):
    res = 0
    for g in ['id','im','ic']:
        res+= orr(n,g,top,False)
    return res


def sum_gates_AND(ANDfunc,ns):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
        res = 0
        for g in ['d','m']:
            res+= ANDfunc(n,g,"all",False)
        gates[i] = res
    return gates

def sum_idle_AND(ANDfunc,ns):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
        res = 0
        for g in ['id','im','ic']:
            res+= ANDfunc(n,g,"all",False)
        gates[i] = res
    return gates

def sum_gates_UCG_QSP(ucg,ns,ds,top):
    gates = np.zeros(len(ns))
    for i , (n,d) in enumerate(zip(ns,ds)):
        res = 0
        for g in ['d','m']:
            res+= sparse_UCG_QSP(ucg,n,d,g,top)
        gates[i] = res
    return gates

def sum_idle_UCG_QSP(ucg,ns,ds,top):
    gates = np.zeros(len(ns))
    for i , (n,d) in enumerate(zip(ns,ds)):
        res = 0
        for g in ['id','im','ic']:
            res+= sparse_UCG_QSP(ucg,n,d,g,top)
        gates[i] = res
    return gates


def sum_gates_unary_based(ns,ds,top, dense = False):
    gates = np.zeros(len(ns))
    if dense:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['d','m']:
                res+= unary_based_QSP(n,d,g,top,False)
            gates[i] =round(res)
    else:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['d','m']:
                res+= unary_based_QSP(n,d,g,top,False)
            gates[i] =round(res)
    return gates

def sum_idle_unary_based(ns,ds,top, dense = False):
    gates = np.zeros(len(ns))
    if dense:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['id','im','ic']:
                res+= unary_based_QSP(n,d,g,top,False)
            gates[i] =round(res)
    else:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['id','im','ic']:
                res+= unary_based_QSP(n,d,g,top,False)
            gates[i] =round(res)
    return gates

def sum_gates_dense_UCG_QSP(ucg,ns,top):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
         
        res = 0
        for g in ['d','m']:
            res+= dense_UCG_QSP(ucg,n,2**n,g,top)
             
        gates[i] = res
    return gates

def sum_idle_dense_UCG_QSP(ucg,ns,top):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
         
        res = 0
        for g in ['id','im','ic']:
             
            res+= dense_UCG_QSP(ucg,n,2**n,g,top)
        gates[i] = res
    return gates