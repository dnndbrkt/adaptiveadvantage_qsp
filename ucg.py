from sympy import *
import numpy as np
import re
import math
import matplotlib.pyplot as plt
from simplify import *
from sympy.utilities.lambdify import lambdastr

gamma_brisbane = 1.936
gamma_fez = 2.542
gamma_strasbourg = 2.201
gamma_torino = 4.067
cyan = "#a3d6e3"
beige = "#cfc7a9"
grey = "#919191"
blue = "#2b6c7b"
red = "#5a002c"
yellow = "#ffdb58"
green = "#006400"



def sum_gates_ucg_qsp(ucg,ns,ds,top):
    gates = np.zeros(len(ns))
    for i , (n,d) in enumerate(zip(ns,ds)):
        res = 0
        for g in ['d','m']:
            res+= ucg_sparse_QSP(ucg,n,d,g,top)
        gates[i] = res
    return gates

def sum_idle_ucg_qsp(ucg,ns,ds,top):
    gates = np.zeros(len(ns))
    for i , (n,d) in enumerate(zip(ns,ds)):
        res = 0
        for g in ['id','im','ic']:
            res+= ucg_sparse_QSP(ucg,n,d,g,top)
        gates[i] = res
    return gates


def sum_gates_un_to_bin(ns,ds,top, dense = False):
    gates = np.zeros(len(ns))
    if dense:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['d','m']:
                res+= un_to_bin_QSP(n,d,g,top,False)
            gates[i] =round(res)
    else:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['d','m']:
                res+= un_to_bin_QSP(n,d,g,top,False)
            gates[i] =round(res)
    return gates

def sum_idle_un_to_bin(ns,ds,top, dense = False):
    gates = np.zeros(len(ns))
    if dense:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['id','im','ic']:
                res+= un_to_bin_QSP(n,d,g,top,False)
            gates[i] =round(res)
    else:
        for i , (n,d) in enumerate(zip(ns,ds)):
            res = 0
            for g in ['id','im','ic']:
                res+= un_to_bin_QSP(n,d,g,top,False)
            gates[i] =round(res)
    return gates

def un_to_bin_QSP(n,d,g,top,idle):
    return  (1/(d**(1/100)))*(unary_dataloader(d,g,top,idle) + un_to_bin(n,d,g,top,idle))

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

def allcock(n,d,g,top,idle):
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

def gr(n,d,g,top,idle):
    # print(f"gr n = {n}, {idle}")
    if n == 1:
        return 0
    else:
        eqs = 2**(n-1)*eq(n,g,top,idle)
        cus = 2**(n-1)*cu(g,idle)
    return eqs + cus

def ucg_dense_QSP(ucg,n,d,g,top):
    res = 0
    for i in range(2,n+1):
        res+= (1/(d**(1/8)))*ucg(i,d,g,top,False)
    # print(f"d= {d} n = {n} g = {g} top = {top} \n")
    # print(f"test")
    # print(res)
    for i in range(2,n):
        res+= (1/(d**(1/8)))*ucg(i,d,g,top,True)
    # print(f"test 2")
    # print(res)
    # print("\n")
    return res

def ucg_dense2_QSP(ucg,n,d,g,top):
    res = 0
    res+= Sum(ucg(i,d,g,top,False),(i,2,n)).evalf(chop=True)
    # print(f"dense res \n")
    # print(f"d= {d} n = {n} g = {g} top = {top} \n")
    # print(f"test")
    # print(res)
    res+= Sum(ucg(i,d,g,top,True),(i,2,n-1)).evalf(chop=True)
    # print(f"test 2")
    # print(res)
    # print("\n")
    return res

def sum_gates_dense_ucg_qsp(ucg,ns,top):
    gates = np.zeros(len(ns))
    
    for i , n in enumerate(ns):
        print(n)
        res = 0
        for g in ['d','m']:
            res+= ucg_dense_QSP(ucg,n,2**n,g,top)
            print(f"g = {g}")
        gates[i] = res
    return gates

def sum_idle_dense_ucg_qsp(ucg,ns,top):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
        print(n)
        res = 0
        for g in ['id','im','ic']:
            print(f"g = {g}")
            res+= ucg_dense_QSP(ucg,n,2**n,g,top)
        gates[i] = res
    return gates


def ucg_sparse_QSP(ucg,n,d,g,top):
    ucg_max = ceiling(log(d,2))
    dense_res = ucg_dense_QSP(ucg,ucg_max,2**ucg_max,g,top)
    res = (1/(d**(1/8)))*perm(n,d,g,top,False)
    return dense_res +res



def unbin_comp_qsp(gate_qsp,idle_qsp,ucg = allcock):
    # values for unary_to_binary
    tops = ['1d','2d']
    labels = ['1D grid', '2D grid']
    ns_all = ([50,50,20,20],[1000,1000,1000,200])
    for i,(label,top) in enumerate(zip(labels,tops)):
        nss =[np.linspace(10,ns,50) for ns in ns_all[i]]
        d0 = ([10]*len(nss[0]),"10")
        d1 = (nss[1],"n")
        d2 = (nss[2]**2,"n^2")
        # if top != "all":
        #     d3 = (nss[3]**3, "n^3")
        # else:
        #     d3 = (nss[3]*0,"0")
        dss = [d0,d1,d2]
        for j, (ns, ds) in enumerate(zip(nss, dss)):
            gates_laqcc = gate_qsp(ucg,ns,ds[0],"LAQCC")
            idles_laqcc = idle_qsp(ucg,ns,ds[0],"LAQCC")
            print(f"d = {ds[1]}")
            
            print(top + "\n")
            
            gates_top = gate_qsp(ucg,ns,ds[0],top)
            idles_top = idle_qsp(ucg,ns,ds[0],top)
            print(gates_top)
            print(idles_top)
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            if top != "all":
                _, ax = plt.subplots()
                ax.plot(ns,gamma_2,'red')
                ax.plot(ns,-1*gamma_1,'blue')
                plt.show()
            gamma = gamma_1/gamma_2
            print("gamma")
            print("\n")
            print(gamma)
            _, ax = plt.subplots()
            ax.set_xlim(min(ns),max(ns))
            ax.set_ylim(min(gamma), max(3, max(gamma)))
            ymin,ymax = ax.get_ylim()
            if top != "1d":
                ax.set_yticks(list(range(int(ymin),int(ymax)+2)))
            # # ax.set_ylim(, max(exps))
            ax.fill_between(ns,ax.get_ylim()[0],0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
            ax.fill_between(ns,0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
            ax.axhline(y=gamma_brisbane,color = beige,linestyle = 'dashed', label = "IBM Brisbane")
            ax.plot(ns,gamma,color = blue,label= label)
            ax.set_xlabel(r"$n$")
            ax.set_ylabel(r"$\gamma_1\ /\ \gamma_2$")
            ax.set_title(f"d = {ds[1]}")
            if top == "all":
                ax.set_xscale("log")
            plt.legend()
            plt.show()


def comp_dense_unbin():
    #values for allcock
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = [50,20]    
    for i,(label,top) in enumerate(zip(labels,tops)):
        _, ax = plt.subplots()
        ns =[int(x) for x in np.linspace(5,ns_all[i],8)]
     
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
        print(label)
        print(ns)
        
        
        if top == "all":
            ysum = 360
            ymin = -350
            ymax = 10
        if top == "2d":
            ysum =12
            ymin = -2
            ymax = 10
        ax.set_ylim((ymin,ymax))
        ax.set_xlim(5,20)
        ax.fill_between(ns,ymin,0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
        ax.fill_between(ns,0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
        # ax.set_yticks(list(range(ymin,ymax+1)))
        # if top == "all":
        #     ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        #     ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        #if top  == "2d":
        ax.annotate("ibm_brisbane", xy = (0.2,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        ax.annotate("ibm_torino",  xy = (0.2,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)

        # else:
        #     ax.annotate("ibm_brisbane", xy = (0.9,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        
        # ax.annotate("ibm_torino",  xy = (0.1,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        # morecolors = [red,yellow, blue]
        # if top != "all":
        #     morecolors = morecolors + [green, "orange"]
        # for j, ((ns, ds),color) in enumerate(zip(zip(nss, dss), morecolors)):
        #     if top == "all":
        #         ax.set_ylim(-2, 10)
        #         ax.set_xlim(2**3,2**50+1)
        #     if top == "1d":
        #         if j > 1:
        #             ax.set_ylim(10)
        #         ax.set_ylim(-2,6)
        #     if top == "2d":
        #         ax.set_ylim(-2,10)
        #     ymin,ymax = ax.get_ylim()
        #     # if top != "1d":
        #     #     ax.set_yticks(list(range(int(ymin),int(ymax)+1)))
        #     if top == "all":
        #         ax.set_xscale("log", base =2)
        #         ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        #         ax.set_xticks([2**3,2**9, 2**15, 2**21, 2**33, 2])
        #     if top == "all":
        #         xticks = ax.get_xticks()
        # if top == "2d":
            
        #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            # ax.set_xticks([2xticks)
            # # ax.set_ylim(, max(exps))

        ds = [2**n for n in ns]  
        gates_laqcc = sum_gates_un_to_bin(ns,ds,"LAQCC", dense = True)
        idles_laqcc = sum_idle_un_to_bin(ns,ds,"LAQCC",dense = True)
        gates_top = sum_gates_un_to_bin(ns,ds,top, dense = True)
        idles_top = sum_idle_un_to_bin(ns,ds,top, dense = True)
        print(gates_top)
        print(idles_top)
        gamma_2=  gates_laqcc- gates_top
        gamma_1 = idles_top - idles_laqcc
        gamma = gamma_1/gamma_2
        print(gamma)
        ax.plot(ns,gamma,color =blue)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        # ax.set_yscale('log')
        if i ==0:
            plt.legend(loc = 'upper right')
        else:
            plt.legend(loc = 'upper right')
        plt.show()



def comp_dense_ucg_qsp():
    #values for allcock
    UCGS = ['Parallelized', 'Sequential']
    ns_all = [990,990]    
    for i,(ucglabel,ucg) in enumerate(zip(UCGS,[allcock,gr])):
        _, ax = plt.subplots()
        ns =[int(x) for x in np.linspace(5,ns_all[i],8)]
        ax.fill_between(ns,-3,0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
        ax.fill_between(ns,0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
        print(ucglabel)
       
        ysum = 13
        ymin = -3
        ymax = 10
        ax.set_ylim((ymin,ymax))
        ax.set_xlim(5,1000)
        ax.set_xticks([5,50] + list(range(100,1100,100)))
        ax.set_yticks(list(range(ymin,ymax+1)))
        if ucglabel == "Parallelized":
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        if ucglabel  == "Sequential":
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)

        # else:
        #     ax.annotate("ibm_brisbane", xy = (0.9,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        
        # ax.annotate("ibm_torino",  xy = (0.1,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        # morecolors = [red,yellow, blue]
        # if top != "all":
        #     morecolors = morecolors + [green, "orange"]
        # for j, ((ns, ds),color) in enumerate(zip(zip(nss, dss), morecolors)):
        #     if top == "all":
        #         ax.set_ylim(-2, 10)
        #         ax.set_xlim(2**3,2**50+1)
        #     if top == "1d":
        #         if j > 1:
        #             ax.set_ylim(10)
        #         ax.set_ylim(-2,6)
        #     if top == "2d":
        #         ax.set_ylim(-2,10)
        #     ymin,ymax = ax.get_ylim()
        #     # if top != "1d":
        #     #     ax.set_yticks(list(range(int(ymin),int(ymax)+1)))
        #     if top == "all":
        #         ax.set_xscale("log", base =2)
        #         ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        #         ax.set_xticks([2**3,2**9, 2**15, 2**21, 2**33, 2])
        #     if top == "all":
        #         xticks = ax.get_xticks()
        # if top == "2d":
            
        #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            # ax.set_xticks([2xticks)
            # # ax.set_ylim(, max(exps))
        
        tops = ['all','2d']
        labels = ['All-to-all', '2D grid']
        for j, (top,colorg) in enumerate(zip(tops, [blue,green])):
            if top == '2d' and ucglabel == "Parallelized":
                ns = [int(x) for x in np.linspace(5,50,10)]
            else:
                ns =[int(x) for x in np.linspace(5,ns_all[i],8)]
            print(f"ns = {ns}")
            print(f"1 \n")
            gates_laqcc = sum_gates_dense_ucg_qsp(ucg,ns,"LAQCC")
            print(f"2 \n")
            idles_laqcc = sum_idle_dense_ucg_qsp(ucg,ns,"LAQCC")
            print(f"3 \n")
            gates_top = sum_gates_dense_ucg_qsp(ucg,ns,top)
            print(f"4 \n")
            idles_top = sum_idle_dense_ucg_qsp(ucg,ns,top)
            print(gates_top)
            print(idles_top)
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            gamma = gamma_1/gamma_2
            print(gamma)
            ax.plot(ns,gamma,color =colorg,label= labels[j])
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        plt.legend(loc = "upper right")
        plt.savefig(f"denseUCG_{top}.png", dpi =600)

def comp_unbin():
    #values for allcock
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = ([1000,1000,1000,200], [1000,1000,1000,200])
    for i,(label,top) in enumerate(zip(labels,tops)):
        print(label)
        nss =[np.linspace(5,ns,10) for ns in ns_all[i]]
        d0 = ([10]*len(nss[0]),"10")
        d1 = (nss[1],"n")
        d2 = (nss[2]**2,"n^2")
        # if top != "all":
        d3 = (nss[3]**3, "n^3")
        # if top != "all":
        #     d4 = (2**nss[4], "2^n")
        dss = [d0,d1,d2,d3]
        # if top != "all":
        #     dss = dss + [d4]
        _, ax = plt.subplots()
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
        if top == "2d":
            ymin = -2
            ymax = 10
            ysum = ymax - ymin
        if top == "all":
            ysum = 9
            ymin = -4
            ymax = 5

        ax.fill_between(nss[0],ymin,0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
        ax.fill_between(nss[0],0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
        
        ax.set_ylim((ymin,ymax))
        ax.set_xlim(5,ns_all[i][0])
        ax.set_yticks(list(range(ymin,ymax+1)))
        xticks = [5] + list(range(100,ns_all[i][0] + 100,100))
        ax.set_xticks(xticks)
        if top =="2d":
            ax.annotate("ibm_brisbane", xy = (0.1,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            ax.annotate("ibm_torino",  xy = (0.1,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        else:
            ax.annotate("ibm_brisbane", xy = (0.1,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            ax.annotate("ibm_torino",  xy = (0.1,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        morecolors = [red,yellow, blue]
        morecolors = morecolors + [green, "orange"]
        for j, ((ns, ds),color) in enumerate(zip(zip(nss, dss), morecolors)):
            # if top == "all":
            #     ax.set_ylim(-2, 10)
            #     ax.set_xlim(2**3,2**50+1)
            # if top == "1d":
            #     if j > 1:
            #         ax.set_ylim(10)
            #     ax.set_ylim(-2,6)
            # if top == "2d":
            #     ax.set_ylim(-2,10)
            # ymin,ymax = ax.get_ylim()
            # # if top != "1d":
            # #     ax.set_yticks(list(range(int(ymin),int(ymax)+1)))
            # if top == "all":
            #     ax.set_xscale("log", base =2)
            #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            #     ax.set_xticks([2**3,2**9, 2**15, 2**21, 2**33, 2])
            # if top == "all":
            #     xticks = ax.get_xticks()
            # if top == "2d":
            #     ax.set_xticks([5] + list(range(100,1050,100)))
            #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            # # ax.set_xticks([2xticks)
            # # # ax.set_ylim(, max(exps))

            # for ucg,colorg in zip([allcock,gr], grcolors):
            gates_laqcc = sum_gates_un_to_bin(ns,ds[0],"LAQCC")
            idles_laqcc = sum_idle_un_to_bin(ns,ds[0],"LAQCC")
            
            gates_top = sum_gates_un_to_bin(ns,ds[0],top)
            idles_top = sum_idle_un_to_bin(ns,ds[0],top)
            print(gates_top)
            print(idles_top)
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            gamma = gamma_1/gamma_2
            
            print(f"gamma = {gamma}")
            
            ax.plot(ns,gamma,color =color,label= f"d = {ds[1]}")
            ax.set_xlabel(r"$n$")
            ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        if i ==0:
            plt.legend(loc = 'upper right')
        else:
            plt.legend(loc = 'upper center')
        plt.show()


def comp_ucg_qsp():
    #values for allcock
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = ([1000,1000,1000,1000], [1000,1000,1000,200])
    ucglabels = ["parallelized","sequential"]
    for i,(label,top) in enumerate(zip(labels,tops)):
        for ucg, ucglabel in zip([allcock,gr],ucglabels):
            print(label)
            nss =[np.linspace(5,ns,10) for ns in ns_all[i]]
            d0 = ([10]*len(nss[0]),"10")
            d1 = (nss[1],"n")
            d2 = (nss[2]**2,"n^2")
            # if top != "all":
            d3 = (nss[3]**3, "n^3")
            # if top != "all":
            #     d4 = (2**nss[4], "2^n")
            dss = [d0,d1,d2,d3]
            # if top != "all":
            #     dss = dss + [d4]
            _, ax = plt.subplots()
            ax.fill_between(nss[0],-3,0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
            ax.fill_between(nss[0],0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
            ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
            ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
            ysum = 12
            ymin = -2
            ymax = 10
            ax.set_ylim((ymin,ymax))
            ax.set_xlim(5,ns_all[i][0])
            ax.set_yticks(list(range(ymin,ymax+1)))
            xticks = [5] + list(range(100,ns_all[i][0] + 100,100))
            ax.set_xticks(xticks)
            if top =="2d":
                ax.annotate("ibm_brisbane", xy = (0.4,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
                ax.annotate("ibm_torino",  xy = (0.2,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            else:
                ax.annotate("ibm_brisbane", xy = (0.1,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
                ax.annotate("ibm_torino",  xy = (0.1,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
            morecolors = [red,yellow, blue]
            morecolors = morecolors + [green, "orange"]
            
            for j, ((ns, ds),color) in enumerate(zip(zip(nss, dss), morecolors)):
            # if top == "all":
            #     ax.set_ylim(-2, 10)
            #     ax.set_xlim(2**3,2**50+1)
            # if top == "1d":
            #     if j > 1:
            #         ax.set_ylim(10)
            #     ax.set_ylim(-2,6)
            # if top == "2d":
            #     ax.set_ylim(-2,10)
            # ymin,ymax = ax.get_ylim()
            # # if top != "1d":
            # #     ax.set_yticks(list(range(int(ymin),int(ymax)+1)))
            # if top == "all":
            #     ax.set_xscale("log", base =2)
            #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            #     ax.set_xticks([2**3,2**9, 2**15, 2**21, 2**33, 2])
            # if top == "all":
            #     xticks = ax.get_xticks()
            # if top == "2d":
            #     ax.set_xticks([5] + list(range(100,1050,100)))
            #     ax.set_yticks([-2,-1,0,1,2,3,4,5,6,7,8,9,10])
            # # ax.set_xticks([2xticks)
            # # # ax.set_ylim(, max(exps))

            
                gates_laqcc = sum_gates_ucg_qsp(ucg,ns,ds[0],"LAQCC")
                idles_laqcc = sum_idle_ucg_qsp(ucg,ns,ds[0],"LAQCC")
                
                gates_top = sum_gates_ucg_qsp(ucg,ns,ds[0],top)
                idles_top = sum_idle_ucg_qsp(ucg,ns,ds[0],top)
                print(gates_top)
                print(idles_top)
                gamma_2=  gates_laqcc- gates_top
                gamma_1 = idles_top - idles_laqcc
                gamma = gamma_1/gamma_2
                print(gamma)
                ax.plot(ns,gamma,color =color,label= f"d = {ds[1]}")
                ax.set_xlabel(r"$n$")
                ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
            plt.legend(loc = "upper right")
            plt.savefig(f"sparseUCG_{ucglabel}_{top}.png",dpi=600)

if __name__ == "__main__":
    comp_ucg_qsp()
    # comp_dense_ucg_qsp()
