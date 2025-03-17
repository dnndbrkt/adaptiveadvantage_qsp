from sympy import *
import numpy as np
import re
import math
import matplotlib.pyplot as plt




n,p,logp, logp_1, logn, sqrtn,two_p, d,i = symbols('n p lp lpo ln sqrtn two_p d i')

def andlogbetter(n,g,base,idle):
    k = floor(log(n,2))-2
    if not idle:
        left = (2**k-1)*(andd(4,g,False)+andd(5,g,False)) + (n*k-3*(2**k - 1))*(andd(5,g,True))+ (n*k-2*(2**k - 1))*(andd(4,g,True))+(2**k)*(andd(base,g,False)+ andd(base,g,True))
        right = (2**k-1)*(andd(5,g,False)) + (n*k-3*(2**k - 1))*(andd(5,g,True))+(2**k)*(andd(base,g,False) + andd(base,g,True))
        return left + right
    else:
        left = k*(andd(5,g,idle)+andd(4,g,idle) + andd(base,g,idle))
        right = k*(andd(5,g,idle) + andd(base,g,idle))
        return left + right
    

def andlog(n,g,idle):
    if not idle:
        andfour = (n/4-1)*(andd(4,g,False)) + (2*n*log(n,2)-6*n+8)*(andd(4,g,True))
        andfive = (n-1)*(andd(3,g,False)) + (2*n*log(n,2)-(11/2)*n+6)*(andd(3,g,True))
        return andfour + andfive

def andlogexact(n_int,g,idle):
    if n_int <3:
        return 0
    if n_int >= 3 and n_int <= 5:
        return andd(n_int,g,idle)
    else:
        left = andlogexact(5,g,False) + (n_int-4)*andlogexact(5,g,True) + andlogexact(4,g,False) + (n_int-3)*andlogexact(4,g,True)
        right = andlogexact(5,g,False) + (n_int-4)*andlogexact(5,g,True)
        recursive_step = andlogexact(floor((n_int-1)/2),g,False) + andlogexact(ceiling((n_int-1)/2),g, False)
        return left+right+recursive_step


def andd(n,g,idle):
    if n == 3:
        if idle:
            return pid(g,6)
        else:
            return pd(g,6) + pid(g,5)
    if n==4:
        if idle:
            return 3*andd(3,g,True)
        else:
            return 3*andd(3,g,False) + 6*andd(3,g,True)
    if n==5:
        if idle:
            return 2*andd(4,g,True)+andd(3,g,True)
        else:
            return 2*andd(4,g,False)+andd(3,g,False)+2*andd(4,g,True) + 2*andd(3,g,True)
        
def orr(n,g,top,idle):
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
    
def orruit(n,g,top,idle):
    p = log(n,2)
    two_p = 2**p
    if not idle:
        red = 2*orreducuit(n-1,g,top,idle)
        fos = (2*p+1)*fo(two_p-1,g,top,idle) + 2*(two_p-1)*fo(p+1,g,top,idle)
        ghz = GHZ(two_p-1,g,top,idle)
        cus = (two_p-1)*cu(g,idle)+(p*two_p-1)*cu(g,True)
        return red + fos + ghz + cus
    if idle:
        red = 2*orreducuit(n-1,g,top,idle)
        fos = 3*fo(two_p-1,g,top,idle) + 2*fo(p+1,g,top,idle)
        cus = cu(g,idle)
        return red+fos+cus

# def orr2(n,g,top,idle):
#     if not idle:
#         red = 2*orreduc(n-1,p,g,top,idle)
#         fos = (2*p+1)*fo(two_p-1,g,top,idle) + 2*(two_p-1)*fo(p+1,g,top,idle)
#         ghz = GHZ(two_p-1,g,top,idle)
#         cus = (two_p-1)*cu(g,idle)+(p*two_p-1)*cu(g,True)
#         fos2 = 3*(n-1)*fo(two_p-1,g,top,True) + 2*(n-1)*fo(p+1,g,top,True)
#         cus2 = (n-1)*cu(g,True)
#         retu
#rn red + fos + ghz + cus + fos2 + cus2
#     if idle:
#         red = 2*orreduc(n-1,p,g,top,idle)
#         fos = 3*fo(two_p-1,g,top,idle) + 2*fo(p+1,g,top,idle)
#         cus = cu(g,idle)
#         return red+fos+cus

def orreducuit(n,g,top,idle):
    p = ceiling(log(n,2))
    fos = 2*p*fo(n,g,top,idle) + 2*n*fo(p,g,top,idle)
    return simplify(fos +n*p*cu(g,idle))

def orreduc(n,p,g,top,idle):
    if not idle:
        fos = 2*p*fo(n,g,top,idle) + 2*n*fo(p,g,top,idle)
        return fos +n*p*cu(g,idle)
    else:
        fos = 2*fo(n,g,top,idle) + 2*fo(p,g,top,idle)
        return ((fos + cu(g,idle)))

def eq(n,g,top,idle):
    if top == "all":
        return andlogbetter(n,g,3,idle)
    else:
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

def gr(n,d,g,top,idle):
    if n == 1:
        return 0
    else:
        eqs = 2**(n-1)*eq(n,g,top,idle)
        cus = 2**(n-1)*cu(g,idle)
    return eqs + cus





def data_ucg_QSP():
    print("allcock")
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = ucg_dense_QSP(allcock,n,2**n,g,top)
            #print(f"g = {g} \t top = {top} \t sum = {expand(res)}")
    print("gr")
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = ucg_dense_QSP(gr,n,d,2**n,top)
            print(f"g = {g} \t top = {top} \t sum = {expand(res)}")


def sum_gates_perm(n,d,top):
    res = 0
    for g in ['d','m']:
        res+= perm(n,d,g,top,False)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res 

def sum_idle_perm(n,d,top):
    res = 0
    for g in ['id','im','ic']:
        res+= perm(n,d,g,top,False)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res 

def sum_gates_FO(n,top):
    res = 0
    for g in ['d','m']:
        res+= fo(n,g,top,False)
        # res = res.subs(two_p,2**p)
        # res = res.subs(p,ceiling(log(n,2)))
    return res

def sum_idle_FO(n,top):
    res = 0
    for g in ['id','im','ic']:
        res+= fo(n,g,top,False)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res


def sum_gates_OR(n,top):
    res = 0
    for g in ['d','m']:
        res+= orr(n,g,top,False)
        # res = res.subs(two_p,2**p)
        # res = res.subs(p,ceiling(log(n,2)))
    return res

def sum_idle_OR(n,top):
    res = 0
    for g in ['id','im','ic']:
        res+= orr(n,g,top,False)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res


def sum_gates_allcock_ucg(n,d,top,idle):
    res = 0
    for g in ['d','m']:
        res+= allcock(n,d,g,top,idle)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res

def sum_idle_allcock_ucg(n,d,top,idle):
    res = 0
    for g in ['id','im','ic']:
        res+= allcock(n,d,g,top,idle)
        # res = res.subs(two_p,2**p)
        # res =res.subs(p,ceiling(log(n,2)))
    return res

def sum_gates_un_to_bin(d,d_new,top):
    res = 0
    for g in ['d','m']:
        res+= un_to_bin_QSP(n,d_new,g,top,False)
    return res

def sum_idle_un_to_bin(d,d_new,top):
    res = 0
    for g in ['id','im','ic']:
        res+= un_to_bin_QSP(n,d_new,g,top,False)
    return res

def plotAllcock():
    ns = np.linspace(10,500)
    d = 100
    colors = ['r','b','g','y']
    tops = ['all','1d','2d','LAQCC']
    for color,top in zip(colors,tops):
        gate_func = lambdify(n, sum_gates_allcock(p,two_p,d,top),modules = ['numpy'])
        gates = gate_func(ns)
        plt.plot(ns,gates,color,label = top)
    plt.legend()
    plt.show()
    for color,top in zip(colors,['all','1d','2d','LAQCC']):
        idle_func = lambdify(n, sum_idle_allcock(p,two_p,d,top),modules = ['numpy'])
        idles = idle_func(ns)        
        plt.plot(ns,idles,color+'-',label = top)
    plt.legend()
    plt.show()

def plotOR():
    ns = np.linspace(5,500)
    colors = ['r','b','g','y']
    tops = ['all','1d','2d','LAQCC']
    for color,top in zip(colors,tops):
        gate_func = lambdify(n, sum_gates_OR(n,top),modules = ['numpy'])
        gates = gate_func(ns)
        if top != '1d':
            plt.plot(ns,gates,color,label = top)
    plt.legend()
    plt.show()
    for color,top in zip(colors,['all','1d','2d','LAQCC']):
        idle_func = lambdify(n, sum_idle_OR(n,top),modules = ['numpy'])
        idles = idle_func(ns)        
        if top != '1d':
            plt.plot(ns,idles,color+'-',label = top)
    plt.legend()
    plt.show()



def data_un_to_bin_QSP():
    print("unary dataloader")
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = unary_dataloader(d,g,top,False)
            print(f"g = {g} \t top = {top} \t sum = {expand(res)}")
    print("unary to binary")
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = un_to_bin(n,d,g,top,False)
            res = res.subs(two_p,2**p)
            res = res.subs(p,log(n))
            res = res.subs(2**log(n),n)
            print(f"g = {g} \t top = {top} \t sum = {expand(res)}")
    print("combination!")
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = un_to_bin_QSP(n,d,g,top,False)
            res = res.subs(two_p,2**p)
            res = res.subs(p,log(n))
            res = res.subs(2**log(n),n)
            print(f"g = {g} \t top = {top} \t sum = {expand(res)}")

def plot_un_to_bin_QSP():
    ns = np.linspace(10,100)

    colors = ['r','b','g','y']
    tops = ['all','1d','2d','LAQCC']
    for color,top in zip(colors,tops):
        gate_func = lambdify(n, sum_gates_un_to_bin(d,d_new,top),modules = ['numpy'])
        gates = gate_func(ns)
        if top != "1d":
            plt.plot(ns,gates,color,label = top)
    plt.legend()
    plt.show()
    for color,top in zip(colors,['all','1d','2d','LAQCC']):
        idle_func = lambdify(n, sum_idle_un_to_bin(d,d_new,top),modules = ['numpy'])
        idles = idle_func(ns)
        if top != "1d":
            plt.plot(ns,idles,color+'-',label = top)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # data_ucg_QSP()
    #compFO()
    #compOR()
    # compANDnorm()
    #compEQ()
    # compANDnorm()
    #comp_allcock()



    #allcock
    # m = Symbol('m')
    # thing = 4*n -16 + ((m-1)/2-1)*(2*n-16) + 4*((m-1)/2-1)*((m-1)/2)
    # print(expand(thing))
    for g in ['d','m','id','im','ic']:
        for top in ['all','1d','2d','LAQCC']:
            res = expand(orreduc(n,p,g,top,False))
            print(f"g = {g} \t top = {top} \t sum = {res}")
    #         # res = res.subs(two_p,2**p)
    #         # res =res.subs(p,log(n))
    #         # res = res.subs(2**log(n),n)
    #         string = str(res)
    #         # string = re.sub(r"\*\*","^",string)
    #         # string = re.sub(r"\*","",string)
    #         # string = re.sub(r"ceiling",r"\\lceil",string)
    #         # string = re.sub(r"two_p",r" 2^p",string)
    #         # string = re.sub(r"log\(([^\)]*)\)",r"\\log(\1)",string)
    #         # string = re.sub(r"sqrt\(([^\)]*)\)",r"\\sqrt{\1}",string)         
             
    # print("andlog \n")
    #andlog
    # for g in ['d','m','id','im','ic']:
    #     string = str(andlog(n,g,False))
    #     string = re.sub(r"\*","",string)
    #     # string = re.sub(r"ceiling\((\([a-z0-9]*\)*)\)",r"\\lceil \1 \\rceil",string)
    #     # string = re.sub(r"two_p",r" 2^p",string)
    #     string = re.sub(r"log\(([^\)]*)\)",r"\\log(\1)",string)
    #     string = re.sub(r"sqrt\(([^\)]*)\)",r"\\sqrt{\1}",string)
    #     print(f"g = {g} \t sum = {string}")
    # print("\n")

    # string = 2*p*two_p + 4*p*(-2*n + (n - 1)*ceiling(log(n - 1)) + 4) -  two_p + 4*(n - 1)*(p*ceiling(log(p)) - 2*p + 2)+ (-2*p + (p + 1)*ceiling(logp_1))*(2*two_p - 2) + (2*p + 1)*(-2*two_p + ( two_p - 1)*p + 4) + (two_p - 1)*(p-1) + 2
    # result = expand(string)
    # string = str(result)
    # string = re.sub(r"\*\*","^",string)
    # string = re.sub(r"\*","",string)
    # string = re.sub(r"ceiling",r"\\lceil",string)
    # string = re.sub(r"two_p",r" 2^p",string)
    # string = re.sub(r"log\(([^\)]*)\)",r"\\log(\1)",string)
    # string = re.sub(r"sqrt\(([^\)]*)\)",r"\\sqrt{\1}",string)
    # string = re.sub(r"lpo",r"\\log(p+1)",string)
    # string = re.sub(r"lp",r"\\log(p)",string)
