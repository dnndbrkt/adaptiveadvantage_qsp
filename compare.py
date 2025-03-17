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



n,p,logp, logp_1, logn, sqrtn,two_p, d,i = symbols('n p lp lpo ln sqrtn two_p d i')

def comp_allcock():
    tops = ["LAQCC",'all','1d','2d']
    labels = ["LAQCC",'All-to-all', '1D grid', '2D grid']
    for d in [n,n**2,2**n]:
        print(f"d = {d}")
        gates_LAQCC = sum_gates_ucg_qsp(allcock,n,d,"LAQCC",False)
        idles_LAQCC = sum_idle_ucg_qsp(allcock,n,d,"LAQCC",False)
        for label,top in zip(labels,tops):
            print(top + "\n")
            gates_top = sum_gates_ucg_qsp(allcock,n,d,top,False)
            idles_top = sum_idle_ucg_qsp(allcock,n,d,top,False)
            print(gates_top)
            print("\n")
            print(idles_top)
            print("\n")
            # print(f"hallo")
            # # exponent = (idles_top - idles_LAQCC)/(gates_LAQCC-gates_top)
            # print(f"hallo")
            # print(expand(exponent))
            # print("\n")
            ns = np.linspace(10,30)
            # print(ns)
            # if top == "all":
            #     ns = np.linspace(5,10**(12))
            # if top =="2d":
            #     ns = np.linspace(5,1000)
            # if top == "1d":
            #     ns = np.linspace(5,100)
            # exponent_func = lambdify(n,exponent,modules = ['numpy'])
            gates_func = lambdify(n,gates_top,modules = ['numpy'])
            idles_func = lambdify(n,idles_top,modules = ['numpy'])
            ga = gates_func(ns)
            id = idles_func(ns)
            # print(ga)
            # print(id)
            _, ax = plt.subplots()
            ax.set_xlim(min(ns),max(ns))
            # # ax.set_ylim(, max(exps))
            # ax.fill_between(ns,ax.get_ylim()[0],0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
            # ax.fill_between(ns,0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
            # ax.axhline(y=gamma_brisbane,color = beige,linestyle = 'dashed', label = "IBM Brisbane")
            ax.plot(ns,ga,color = blue,label= label)
            ax.plot(ns,id,color ='red',label= label)
            # ax.set_xlabel(r"$n$")
            # ax.set_ylabel(r"$\gamma_1\ /\ \gamma_2$")
            ax.set_title(f"d = {d}")
            # if top == "all":
            #     ax.set_xscale("log")
            plt.legend()
            plt.show()

def comp_perm():
    ns = np.linspace(5,30)
    N = 2**n
    ds = [5,n,n**2,n**5, N]
    colors = ['r','b','g']
    tops = ['all','1d','2d']
    for d in ds:
        gates_LAQCC = sum_gates_perm(d,'LAQCC')
        idle_LAQCC = sum_idle_perm(d,'LAQCC')
        for color,top in zip(colors,tops):
            gates_top = sum_gates_perm(d,top)
            idle_top = sum_idle_perm(d,top)
            exponent = (idle_top - idle_LAQCC)/(gates_LAQCC-gates_top)
            print(exponent)
            print("\n")
            exponent_func = lambdify((n,i),exponent,modules = ['numpy'])
            exps = exponent_func(ns)
            if top != '1d':
                plt.plot(ns,exps,color,label= top)
        plt.title(f"d = {str(d)}")
        plt.legend()
        plt.show()

def compEQ():
    ns = np.linspace(5,2000)
    pd = 1-9.442*10**(-3)
    pid = 1- 4.998*10**(-3)
    gamma = math.log(pd,pid)
    plt.axhline(y=gamma,color = 'grey',linestyle = 'dashed', label = "IBM Brisbane")
    colors = ['r','b','g']
    tops = ['all','1d','2d']
    gates_LAQCC = sum_gates_OR(n,'LAQCC')
    idle_LAQCC = sum_idle_OR(n,'LAQCC')
    for color,top in zip(colors,tops):
        gates_top = sum_gates_OR(n,top)
        idle_top = sum_idle_OR(n,top)
        exponent = (idle_top - idle_LAQCC)/(gates_LAQCC-gates_top)
        print(f"top = {top}")
        print(expand(exponent))
        print("\n")

        exponent_func = lambdify(n,exponent,modules = ['numpy'])
        
        exps = exponent_func(ns)
        if top != '1d':
            plt.plot(ns,exps,color,label= top)
    plt.xlabel("n")
    plt.ylabel("γ")
    plt.legend()
    plt.show()

def compFO():
    tops = ['all','1d','2d']
    labels = [' ', ' ', ' ']
    gates_LAQCC = sum_gates_FO(n,'LAQCC')
    idle_LAQCC = sum_idle_FO(n,'LAQCC')
    for label,top in zip(labels,tops):
        gates_top = sum_gates_FO(n,top)
        idle_top = sum_idle_FO(n,top)
        exponent = (idle_top - idle_LAQCC)/(gates_LAQCC-gates_top)
        print(expand(exponent))
        print("\n")

        if top == "all":
            ns = np.linspace(5,2**(39))
        if top =="2d":
            ns = np.linspace(5,1000)
        if top == "1d":
            ns = np.linspace(5,35)
        exponent_func = lambdify(n,exponent,modules = ['numpy'])
        exps = exponent_func(ns)
        _, ax = plt.subplots()
        ax.set_xlim(min(ns),max(ns))
        ax.set_ylim(min(exps), 7)
        ymin = -2
        ymax = 7
        ysum = ymax - ymin

        ax.set_yticks(list(range(ymin,ymax+1)))
        ax.fill_between(ns,ax.get_ylim()[0],0,color = grey, alpha = 0.5,label = "Impossible")
        ax.fill_between(ns,0,1,color = grey, alpha = 0.2,label = "Unrealistic")
        
        ax.plot(ns,exps,blue)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        if top == "all":
            ax.set_xscale("log", base=2)
            ax.set_xticks([2**n for n in list(range(3,40,4))])
        if top == "2d":
            ax.set_xticks([5] + list(range(100,1100,100)))
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        if top !="2d":
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        else:
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
        ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        plt.legend()
        
        plt.savefig(f"FO_{top}.png", dpi = (600)) 

def compANDnorm():
    ns = np.linspace(10,1000,1000-10)
    #plt.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed', label = "IBM Brisbane")

    gates_or = orr(n,'d','all',False)
    idle_or = orr(n,'id','all',False)
    # print(expand(gates_or), "\n")
    # print(expand(idle_or), "\n")
    for base in [3,4,5]:
        gates_recursive = andlogbetter(n,'d',base,False)
        idle_recursive = andlogbetter(n,'id',base,False)
        gamma1 = (idle_recursive- idle_or)
        gamma2 = (gates_recursive - gates_or)
        exponent = gamma1/gamma2
        print(f"base = {base} \n")
        print(expand(exponent))
        print("\n")

        gamma1_func = lambdify(n,gamma1,modules = ['numpy'])
        gamma2_func = lambdify(n,gamma2,modules = ['numpy'])
        gamma1s = gamma1_func(ns)
        gamma2s = gamma2_func(ns)
        exponent_func = lambdify(n,exponent,modules = ['numpy'])
        
        exps = exponent_func(ns)
        if base == 3:
            alpha = 1
        else:
            alpha = 0.2
        if base== 3:
            plt.plot(ns,gamma1s,blue,label= r"$|Recursive|_{i_d} -|OR|_{i_d}$" +f", b = {base}", alpha = alpha)
            plt.plot(ns,gamma2s,red,label= r"$|Recursive|_d -|OR|_d$" +f", b = {base}", alpha = alpha)
        else:
            plt.plot(ns,gamma1s,blue, alpha = alpha)
            plt.plot(ns,gamma2s,red, alpha = alpha)
        #plt.plot(ns,exps,blue,label= f"b = {base}", alpha = alpha)
     
    # plt.xlim(min(ns),max(ns))
    # plt.ylim(-2, 6)
    # plt.yticks(list(range(-2,int(plt.ylim()[1])+1)))
    # plt.fill_between(ns,-2,0,color = grey, alpha = 0.5,label = "Impossible")
    # plt.fill_between(ns,0,1,color = grey, alpha = 0.2,label = "Unrealistic")
    plt.xlabel(r"$n$")
    # plt.ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
    plt.legend(loc = "lower left")
    plt.show()


    # ns = np.linspace(30,500,500-30)
    # pd = 1-9.442*10**(-3)
    # pid = 1- 4.998*10**(-3)
    # gamma = math.log(pd,pid)
    # plt.axhline(y=gamma,color = 'grey',linestyle = 'dashed', label = "IBM Brisbane")

    # gates_or = orr(n,'d','all',False)
    # idle_or = orr(n,'id','all',False)
    # # print(expand(gates_or), "\n")
    # # print(expand(idle_or), "\n")
    # for base in [3,4,5]:
    #     gates_recursive = andlogbetter(n,'d',base,False)
    #     idle_recursive = andlogbetter(n,'id',base,False)
    #     print(expand(gates_recursive), "\n")
    #     print(expand(idle_recursive), "\n")
    #     exponent = (idle_or - idle_recursive)/(gates_recursive - gates_or)
    #     print(expand(exponent))
    #     print("\n")
    #     exponent_func = lambdify(n,exponent,modules = ['numpy'])
    #     exps = exponent_func(ns)
    #     plt.plot(ns,exps,label= f"Base = {base}")
    # plt.ylabel("gamma")
    # plt.xlabel("n")
    # plt.legend()
    # plt.show()

def compAllcockQSP():
    ns = np.linspace(10,1000,1000-10)
    plt.axhline(y=gamma_brisbane,color = grey,linestyle = 'dashed', label = "IBM Brisbane")

    # gates_laqcc = sum_gates_allcock(d,g,)



def compAND():
    ns = list(range(10,1000,5))
    print(ns)
    pd = 1-9.442*10**(-3)
    pid = 1- 4.998*10**(-3)
    gamma_ibm = math.log(pd,pid)
    plt.axhline(y=gamma_ibm,color = 'grey',linestyle = 'dashed', label = "IBM Brisbane")

    gammas = np.zeros(len(ns))
    for i,n in enumerate(ns):
        gates_recursive = andlogexact(int(n),'d',False)
        idle_recursive =  andlogexact(int(n),'id',False)
        gates_or =  orr(int(n),'d','all',False)
        idle_or =  orr(int(n),'id','all',False)
        gammas[i] = ((idle_or - idle_recursive)/(gates_recursive-gates_or)).evalf()
        plt.plot(ns,gammas, label = "recursive > or")
    plt.legend()
    plt.show()

def compOR():
    tops = ['all','1d','2d']
    labels = ['All-to-all', '1D grid', '2D grid']
    gates_LAQCC = sum_gates_OR(n,'LAQCC')
    idle_LAQCC = sum_idle_OR(n,'LAQCC')
    for label,top in zip(labels,tops):
        if top =="1d":
            continue
        gates_top = sum_gates_OR(n,top)
        idle_top = sum_idle_OR(n,top)
        exponent = (idle_top - idle_LAQCC)/(gates_LAQCC-gates_top)
        print(expand(exponent))
        print("\n")

        if top == "all":
            ns = np.linspace(5,2**(57))
        if top =="2d":
            ns = np.linspace(5,3000)
        if top == "1d":
            ns = np.linspace(5,50)
        exponent_func = lambdify(n,exponent,modules = ['numpy'])
        
        exps = exponent_func(ns)
        _, ax = plt.subplots()
        ax.set_xlim(min(ns),max(ns))
        ax.set_ylim(-2, 7)
        ymin,ymax = ax.get_ylim()
        ax.set_yticks(list(range(int(ymin),int(ymax)+1)))
        ax.plot(ns,exps,blue)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")        
        if top == "all":
            ax.set_xscale("log",base=2)
            ax.set_xticks([2**n for n in range(3,59,6)])
        if top == "2d":
            ax.set_xticks([5]+list(range(500,3500,500)))
        ysum = 9
        ax.fill_between(ns,ax.get_ylim()[0],0,color = grey, alpha = 0.5,label = "Impossible")
        ax.fill_between(ns,0,1,color = grey, alpha = 0.2,label = "Unrealistic")
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        if top !="2d":
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        else:
            ax.annotate("ibm_brisbane", xy = (0.05,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
        ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        plt.legend(loc ="upper left")
        plt.savefig(f"OR_{top}.png", dpi = (600)) 

if __name__ == "__main__":
    compANDnorm()

    # gates_LAQCC = sum_gates_FO(n,'LAQCC')
    # idle_LAQCC = sum_idle_FO(n,'LAQCC')
    # gates_top = sum_gates_FO(n,'2d')
    # idle_top = sum_idle_FO(n,'2d')
    # exponent = (idle_top - idle_LAQCC)/(gates_LAQCC-gates_top)
    # exponent = exponent.subs(ceiling(log(n,2)), log(n,2))
    # print(solve(exponent - gamma_brisbane,n))
    