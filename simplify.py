from sympy import *
import numpy as np
import re
import math
import matplotlib.pyplot as plt




n = symbols('n')


def decideBase(n): # must be called with n-1 (:
    if n == 3:
        return np.array([1,0,0])
    if n == 4:
        return np.array([0,1,0])
    if n == 5:
        return np.array([0,0,1])
    if n % 2 == 0:
        return 2*decideBase(n/2)
    else:
        return decideBase((n-1)/2) + decideBase((n+1)/2)
    

def and_recursive(n,g,top,idle):
    """Succ. prob of recursive implementation of AND-gate from Nie et al. (2024)"""
    k = floor(log(n,2))-2
    bases = [3,4,5]
    if not idle:
        base_occurrences = decideBase(n-1)
        left = (2**k-1)*(and_base(4,g,False)+and_base(5,g,False)) + (n*k-3*(2**k - 1))*(and_base(5,g,True))+ (n*k-2*(2**k - 1))*(and_base(4,g,True))
        right = (2**k-1)*(and_base(5,g,False)) + (n*k-3*(2**k - 1))*(and_base(5,g,True))
        for base, occurrence in zip(bases,base_occurrences):
            base_ands = occurrence*(and_base(base,g,False) + and_base(base,g,True))
            left += base_ands
            right += base_ands 
        return left + right
    else:
        left = k*(and_base(5,g,idle)+and_base(4,g,idle) + and_base(base,g,idle))
        right = k*(and_base(5,g,idle) + and_base(base,g,idle))
        return left + right

def and_base(n,g,idle):
    """Succ. prob. of AND-gate base cases (n=3,4,5)"""
    if n == 3:
        if idle:
            return pid(g,6)
        else:
            return pd(g,6) + pid(g,5)
    if n==4:
        if idle:
            return 3*and_base(3,g,True)
        else:
            return 3*and_base(3,g,False) + 6*and_base(3,g,True)
    if n==5:
        if idle:
            return 2*and_base(4,g,True)+and_base(3,g,True)
        else:
            return 2*and_base(4,g,False)+and_base(3,g,False)+2*and_base(4,g,True) + 2*and_base(3,g,True)
        
def orr(n,g,top,idle):
    """Succ. prob. of OR-gate"""
    p = (log(n,2))
    if not idle:
        red = 2*orreduc(n-1,p,g,top,idle)
        fos = (2*p+1)*fo(2**p-1,g,top,idle) + 2*(2**p-1)*fo(p+1,g,top,idle)
        ghz = GHZ(2**p-1,g,top,idle)
        cus = (2**p-1)*cu(g,idle)+(p*2**p-1)*cu(g,True)
        return red + fos + ghz + cus
    if idle:
        red = 2*orreduc(n-1,p,g,top,idle)
        ghz = GHZ(2**p-1,g,top,idle)
        fos = 3*fo(2**p-1,g,top,idle) + 2*fo(p+1,g,top,idle)
        cus = cu(g,idle)
        return ((red+fos+cus))

def orreduc(n,p,g,top,idle):
    """Succ. prob of OR-reduction"""
    if not idle:
        fos = 2*p*fo(n,g,top,idle) + 2*n*fo(p,g,top,idle)
        return fos +n*p*cu(g,idle)
    else:
        fos = 2*fo(n,g,top,idle) + 2*fo(p,g,top,idle)
        return ((fos + cu(g,idle)))

def eq(n,g,top,idle):
    """Succ. prob of EQ-gate"""
    return orr(n,g,top,idle)

def cu(g,idle):
    if g =="io":
        return 3
    elif g =="id":
        if idle:
            return 2 
        else:
            return 0
    elif g == "o":
        if idle:
            return 0
        else:
            return 3
    elif g == "d":
        if idle:
            return 0
        else:
            return 3
    else:
        return 0

def pd(g,num):
    if g =="d":
        return num
    else:
        return 0
def pm(g,num):
    if g =="m":
        return num
    else:
        return 0
def pid(g,num):
    if g =="id":
        return num
    else:
        return 0
def pim(g,num):
    if g =="im":
        return num
    else:
        return 0
def pic(g,num):
    if g =="ic":
        return num
    else:
        return 0

def GHZ(n,g,top,idle):
    "Succ. Prob of GHZ-state"
    if top == "LAQCC":
        if idle:
            if g == "id":
                return 2
            if g ==  "im":
                return 1
            if g =="ic":
                return 1
            return 0
        else:
            if g == "d":
                return 2*n-2
            if g == "id":
                return 2
            if g=="m":
                return n-1
            if g =="im":
                return n 
            if g =="ic":
                return n
            return 0
    else:
        #hadamard skipped
        if idle:
            if g == "id":
                return 1 + fo(floor(n/2),g,top,idle)
        else:
            fos =  fo(ceiling(n/2),g,top,idle) +  fo(floor(n/2),g,top,idle)
            fos =  2*fo((n/2),g,top,idle)
            if g == "d":
                return fos +1
            if g == "id":
                return fos + pid(g,n-1)
        return 0

def fo(n,g,top,idle):
    """Succ. prob of FO-gate"""
    if g == "d":
        if idle:
            return 0
        elif top == "LAQCC":
            return 3*n-2
        else:
            return n-1
    if g == "m":
        if top =="LAQCC" and not idle:
            return 2*n-1
        else:
            return 0
    if g == "id":
        if top == "1d":
            if idle:
                return n-1
            else:
                return (n-1)*(n-2)
        if top == "all":
            if idle:
                return ceiling(log(n,2))
            else:
                return n*ceiling(log(n,2))-2*n+2
        if top =="2d":
            if idle:
                return sqrt(n)+1
            else:
                return n*sqrt(n)+2*n-12*sqrt(n)+11
        if top =="LAQCC":
            if idle:
                return 3
            else:
                return 2*n+1
    if g == "im":
        if top =="LAQCC":
            if idle:
                return 2
            else:
                return 3*n-1
        else:
            return 0
    if g == "ic":
        if top =="LAQCC":
            if idle:
                return 2
            else:
                return 3*n-1
        else:
            return 0
    return 0

def rbs(g,idle):
    """Succ. prob of RBS-gate"""
    if idle:
        return pid(g,3)
    else:
        return pd(g,3)


def un_to_bin_QSP(n,d,g,top,idle):
    return unary_dataloader(d,g,top,idle) + un_to_bin(n,d,g,top,idle)

def un_to_bin(n,d,g,top,idle):
    if idle:
        eqs = eq(n+1,g,top,True) + 2*fo(d,g,top,True) + fo(n,g,top,True)
    else:
        fos = 2*n*fo(d,g,top,idle) + d*fo(n,g,top,idle)
        eqs = d*eq(n+1,g,top,idle)
        return fos+ eqs

def perm(n,d,g,top,idle):
    res = un_to_bin(n,d,g,top,idle)
    if idle:
        return res + eq(n+1,g,top,idle) +2*fo(d,g,top,idle)
    else:
        return res + d*(eq(n+1,g,top,idle)) +2*n*fo(d,g,top,idle)

def unary_dataloader(d,g,top,idle):
    if top == "1d" or top == "2d" or top == "LAQCC":
        iter = ceiling((d-1)/2)
    else:
        iter = ceiling(log(d))
    if idle:
        return iter*rbs(g,idle)
    else:
        return (d-1)*rbs(g,False) + (iter*d-2*(d-1))*rbs(g,True)

def sequential_UCG(n,d,g,top,idle):
    if n == 1:
        return 0
    else:
        eqs = 2**(n-1)*eq(n,g,top,idle)
        cus = 2**(n-1)*cu(g,idle)
    return eqs + cus

def unarybased_QSP(n,d,g,top,idle):
    return  (1/(d**(1/100)))*(unary_dataloader(d,g,top,idle) + un_to_bin(n,d,g,top,idle))

def parallelized_UCG(n,d,g,top,idle):
    if n == 0:
        return 0
    elif not idle:
        fosa = 6*fo(d+1,g,top,idle) + 4*(n-1)*fo(d,g,top,idle)
        fosb = 4*(d)*fo(d+1,g,top,True) + 2*(fo(d,g,top,True))
        eqq = 2*d*eq(n,g,top,idle) + 2*eq(n,g,top,True)
        cs = 3*d*cu(g,idle) + 3*n*cu(g,True)
        return fosa+fosb+eqq+cs
    else:
        fosa = 4*fo(d+1,g,top,idle) + 2*fo(d,g,top,idle)
        eqs = 2*eq(n,g,top,idle)
        cs = 3*cu(g,idle)
        return fosa + eqs + cs 

def dense_UCG_QSP(ucg,n,d,g,top):
    res = 0
    for i in range(2,n+1):
        res+= (1/(d**(1/8)))*ucg(i,d,g,top,False)
    for i in range(2,n):
        res+= (1/(d**(1/8)))*ucg(i,d,g,top,True)
    return res

def sparse_UCG_QSP(ucg,n,d,g,top):
    ucg_max = ceiling(log(d,2))
    dense_res = dense_UCG_QSP(ucg,ucg_max,2**ucg_max,g,top)
    res = (1/(d**(1/8)))*perm(n,d,g,top,False)
    return dense_res +res

