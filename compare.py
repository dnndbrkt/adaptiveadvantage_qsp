from sympy import *
from simplify import *
from gatesums import *
import numpy as np
import matplotlib.pyplot as plt

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

def sum_gates_and(andfunc,ns):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
        print(n)
        res = 0
        for g in ['d','m']:
            res+= andfunc(n,g,"all",False)
        gates[i] = res
    return gates

def sum_idle_and(andfunc,ns):
    gates = np.zeros(len(ns))
    for i , n in enumerate(ns):
        res = 0
        for g in ['id','im','ic']:
            res+= andfunc(n,g,"all",False)
        gates[i] = res
    return gates

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

def compAND():
    ns = [int(x) for x in np.linspace(10,1000,50)]
    _, ax = plt.subplots()
    gates_or = sum_gates_and(orr,ns)
    idle_or = sum_idle_and(orr,ns)
    gates_rec = sum_gates_and(and_recursive,ns)
    idle_rec = sum_idle_and(and_recursive,ns)
    gamma1 = (idle_rec- idle_or)
    gamma2 = -1*(gates_rec - gates_or)
    exps = gamma1/gamma2

    ax.plot(ns,exps,blue)
     
    ax.set_xlim(min(ns),max(ns))
    ymin = -2
    ymax = 6
    ysum = ymax - ymin
    ax.set_ylim(ymin,ymax)
    ax.set_yticks(list(range(-2,int(plt.ylim()[1])+1)))
    ax.fill_between(ns,-2,0,color = grey, alpha = 0.5,label = "Impossible")
    ax.fill_between(ns,0,1,color = grey, alpha = 0.2,label = "Unrealistic")
    ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
    
    ax.annotate("ibm_brisbane", xy = (0.25,(gamma_brisbane+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
    ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
    ax.annotate("ibm_torino",  xy = (0.05,(gamma_torino+0.15-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
    ax.legend(loc = "lower left")
    plt.savefig(f"AND_comp.png", dpi = (600)) 

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


def unbin_comp_qsp(gate_qsp,idle_qsp,ucg = parallelized_UCG):
    tops = ['1d','2d']
    labels = ['1D grid', '2D grid']
    ns_all = ([50,50,20,20],[1000,1000,1000,200])
    for i,(label,top) in enumerate(zip(labels,tops)):
        nss =[np.linspace(10,ns,50) for ns in ns_all[i]]
        d0 = ([10]*len(nss[0]),"10")
        d1 = (nss[1],"n")
        d2 = (nss[2]**2,"n^2")
        dss = [d0,d1,d2]
        for j, (ns, ds) in enumerate(zip(nss, dss)):
            gates_laqcc = gate_qsp(ucg,ns,ds[0],"LAQCC")
            idles_laqcc = idle_qsp(ucg,ns,ds[0],"LAQCC")
             
             
            gates_top = gate_qsp(ucg,ns,ds[0],top)
            idles_top = idle_qsp(ucg,ns,ds[0],top)
             
             
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            if top != "all":
                _, ax = plt.subplots()
                ax.plot(ns,gamma_2,'red')
                ax.plot(ns,-1*gamma_1,'blue')
                plt.show()
            gamma = gamma_1/gamma_2
             
             
             
            _, ax = plt.subplots()
            ax.set_xlim(min(ns),max(ns))
            ax.set_ylim(min(gamma), max(3, max(gamma)))
            ymin,ymax = ax.get_ylim()
            if top != "1d":
                ax.set_yticks(list(range(int(ymin),int(ymax)+2)))
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
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = [50,20]    
    for i,(label,top) in enumerate(zip(labels,tops)):
        _, ax = plt.subplots()
        ns =[int(x) for x in np.linspace(5,ns_all[i],8)]    
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
         
         
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
        ax.annotate("ibm_brisbane", xy = (0.2,(gamma_brisbane+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        ax.annotate("ibm_torino",  xy = (0.2,(gamma_torino+0.20-ymin)/ysum), xycoords="axes fraction",horizontalalignment="left",fontsize =8)
        ds = [2**n for n in ns]  
        gates_laqcc = sum_gates_un_to_bin(ns,ds,"LAQCC", dense = True)
        idles_laqcc = sum_idle_un_to_bin(ns,ds,"LAQCC",dense = True)
        gates_top = sum_gates_un_to_bin(ns,ds,top, dense = True)
        idles_top = sum_idle_un_to_bin(ns,ds,top, dense = True)
         
         
        gamma_2=  gates_laqcc- gates_top
        gamma_1 = idles_top - idles_laqcc
        gamma = gamma_1/gamma_2
         
        ax.plot(ns,gamma,color =blue)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        if i ==0:
            plt.legend(loc = 'upper right')
        else:
            plt.legend(loc = 'upper right')
        plt.show()



def comp_dense_ucg_qsp():
    UCGS = ['Parallelized', 'Sequential']
    ns_all = [990,990]    
    for i,(ucglabel,ucg) in enumerate(zip(UCGS,[parallelized_UCG,gr])):
        _, ax = plt.subplots()
        ns =[int(x) for x in np.linspace(5,ns_all[i],8)]
        ax.fill_between(ns,-3,0,color = grey, alpha = 0.5, interpolate = True,label = "Impossible")
        ax.fill_between(ns,0,1,color = grey, alpha = 0.2, interpolate = True,label = "Unrealistic")
        ax.axhline(y=gamma_brisbane,color = "brown",linestyle = 'dashed')
        ax.axhline(y=gamma_torino,color = "brown",linestyle = 'dashed')
         
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

        tops = ['all','2d']
        labels = ['All-to-all', '2D grid']
        for j, (top,colorg) in enumerate(zip(tops, [blue,green])):
            if top == '2d' and ucglabel == "Parallelized":
                ns = [int(x) for x in np.linspace(5,50,10)]
            else:
                ns =[int(x) for x in np.linspace(5,ns_all[i],8)]
             
             
            gates_laqcc = sum_gates_dense_ucg_qsp(ucg,ns,"LAQCC")
             
            idles_laqcc = sum_idle_dense_ucg_qsp(ucg,ns,"LAQCC")
             
            gates_top = sum_gates_dense_ucg_qsp(ucg,ns,top)
             
            idles_top = sum_idle_dense_ucg_qsp(ucg,ns,top)
             
             
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            gamma = gamma_1/gamma_2
             
            ax.plot(ns,gamma,color =colorg,label= labels[j])
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        plt.legend(loc = "upper right")
        plt.savefig(f"denseUCG_{top}.png", dpi =600)

def comp_unbin():
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = ([1000,1000,1000,200], [1000,1000,1000,200])
    for i,(label,top) in enumerate(zip(labels,tops)):
         
        nss =[np.linspace(5,ns,10) for ns in ns_all[i]]
        d0 = ([10]*len(nss[0]),"10")
        d1 = (nss[1],"n")
        d2 = (nss[2]**2,"n^2")
        d3 = (nss[3]**3, "n^3")
        
        dss = [d0,d1,d2,d3]
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
            gates_laqcc = sum_gates_un_to_bin(ns,ds[0],"LAQCC")
            idles_laqcc = sum_idle_un_to_bin(ns,ds[0],"LAQCC")
            
            gates_top = sum_gates_un_to_bin(ns,ds[0],top)
            idles_top = sum_idle_un_to_bin(ns,ds[0],top)
            gamma_2=  gates_laqcc- gates_top
            gamma_1 = idles_top - idles_laqcc
            gamma = gamma_1/gamma_2            
            ax.plot(ns,gamma,color =color,label= f"d = {ds[1]}")
            ax.set_xlabel(r"$n$")
            ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
        if i ==0:
            plt.legend(loc = 'upper right')
        else:
            plt.legend(loc = 'upper center')
        plt.show()


def comp_ucg_qsp():
    tops = ['all','2d']
    labels = ['All-to-all', '2D grid']
    ns_all = ([1000,1000,1000,1000], [1000,1000,1000,200])
    ucglabels = ["parallelized","sequential"]
    for i,(label,top) in enumerate(zip(labels,tops)):
        for ucg, ucglabel in zip([parallelized_UCG,gr],ucglabels):
             
            nss =[np.linspace(5,ns,10) for ns in ns_all[i]]
            d0 = ([10]*len(nss[0]),"10")
            d1 = (nss[1],"n")
            d2 = (nss[2]**2,"n^2")
            d3 = (nss[3]**3, "n^3")
            dss = [d0,d1,d2,d3]
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
                gates_laqcc = sum_gates_ucg_qsp(ucg,ns,ds[0],"LAQCC")
                idles_laqcc = sum_idle_ucg_qsp(ucg,ns,ds[0],"LAQCC")
                
                gates_top = sum_gates_ucg_qsp(ucg,ns,ds[0],top)
                idles_top = sum_idle_ucg_qsp(ucg,ns,ds[0],top)
                gamma_2=  gates_laqcc- gates_top
                gamma_1 = idles_top - idles_laqcc
                gamma = gamma_1/gamma_2
                 
                ax.plot(ns,gamma,color =color,label= f"d = {ds[1]}")
                ax.set_xlabel(r"$n$")
                ax.set_ylabel(r"$\gamma_1(n)\ /\ \gamma_2(n)$")
            plt.legend(loc = "upper right")
            plt.savefig(f"sparseUCG_{ucglabel}_{top}.png",dpi=600)

if __name__ == "__main__":
    comp_ucg_qsp()





if __name__ == "__main__":
    compFO()
    compOR()
    compAND()
    