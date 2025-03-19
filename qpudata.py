import math

def pi(gatetime, dectime):

    #probability of qubit staying idle during a gate, given gate and decoherence time estimates (gatetime in ns, dectime in us)
    return math.exp(-gatetime/(1000*dectime))

def pdCZ(pCZ, ps, pis):
    return pCZ*(ps**2)*(pis**2)

def pdECR(pECR,ps,pis):
    return pECR*(ps**2)*(pis**3)

def gamma(pd,pid):
    return math.log(pd)/math.log(pid)

gatetime_s = 30 #nanoseconds
gatetime_d = [660+2*gatetime_s,68+3*gatetime_s,660+2*gatetime_s,84+3*gatetime_s] #nanoseconds

#Brisbane, Fez, Strasbourg, Torino
# 19-02-2025 17:25
names = ["Brisbane", "Fez","Strasbourg", "Torino"]
processor = ["Eagle r3", "Heron r2","Eagle r3", "Heron r1"]
T1 = [231.16,145.91,303.18,189.9] #us
T2 = [149.51,80.31,159.47, 140.82] #us
ps = [1-2.608*10**(-4),1-2.270*10**(-4),1-2.241*10**(-4),1-3.189*10**(-4)]
pdNative = [1-8.166*10**(-3),1- 3.792*10**(-3),1-8.886*10**(-3),1-3.954*10**(-3)]
dNative = [pdECR,pdCZ,pdECR,pdCZ]

def print_data():
    for i in range(len(names)): 
        print(f"Data for IBM {names[i]} (processor {processor[i]})")
        pis = pi(gatetime_s,T2[i])
        pid = pi(gatetime_d[i],T2[i])
        pd = dNative[i](pdNative[i],ps[i],pis)
        print(f"p_d = \t {pd} \t 1-pd = {(1-pd)*1000}")
        print(f"p_is = \t {pis}\t 1-pis = {(1-pis)*1000}")
        print(f"p_id = \t {pid}\t 1-pid = {(1-pid)*1000}")
        print(f"gamma =\t {gamma(pd,pid)}")
if __name__ == "__main__":
    print_data()