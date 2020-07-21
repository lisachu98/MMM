import matplotlib.pyplot as wkrs
import numpy as np
import cmath
import math

zmienne = {
    "a1": float(0),
    "a0": float(1000000),
    "b2": float(1),
    "b1": float(1000),
    "b0": float(1000000),
    "x1p":float(0),
    "x2p":float(0),
    "T":  float(0),
    "krok": float(0.000001),
    "powt": int(100000),
    "amplituda": float(1),
    "polokres": int(100000),
    "dlugoscSyg": int(100000)
}

class UkladG(object):
    def __init__(self, a1, a0, b2, b1, b0, x1p, x2p, krok):
        self.a1 = a1
        self.a0 = a0
        self.b2 = b2
        self.b1 = b1
        self.b0 = b0
        self.x1p = x1p
        self.x2p = x2p
        self.krok = krok

    def output(self, u):
        x1pom = self.x1p + self.krok * self.x2p
        x2pom = self.x2p + self.krok * (-self.b0/self.b2 * self.x1p - self.b1/self.b2*self.x2p + u)
        y = (self.a0 * self.x1p + self.a1 * self.x2p)/self.b2
        self.x1p = x1pom
        self.x2p = x2pom
        return y

class main():
    def __init__(self, zmienne):
        u1 = np.zeros(zmienne["powt"])
        t = czas()
        r = sygnalP()
        fazy = np.zeros(zmienne["powt"])
        moduly = []
        czestotliwosc = np.zeros(zmienne["powt"])
        for i in range(0, zmienne["powt"]):
            w = i * 10
            jw = (complex(0, zmienne["a1"] * w) + zmienne["a0"]) / (
                        -zmienne["b2"] * w * w + complex(0, zmienne["b1"] * w) + zmienne["b0"]) * np.exp(complex(0, -w*zmienne["T"]))
            fazy[i] = cmath.phase(jw)*(180/math.pi)
            moduly.append(20 * math.log10(abs(jw)))
            czestotliwosc[i] = w
        G1 = UkladG(zmienne["a1"], zmienne["a0"], zmienne["b2"], zmienne["b1"], zmienne["b0"], zmienne["x1p"], zmienne["x2p"], zmienne["krok"])
        for i in range(zmienne["powt"]):
            u1[i]=G1.output(r[i])
        y=opoznianie(u1)
        rysowanie(czestotliwosc, moduly, fazy, t, y, r)

def opoznianie(u1):
    dt = zmienne["T"] / zmienne["krok"]
    dt = round(dt)
    u=np.zeros(zmienne["powt"]+dt)
    for i in range(dt, zmienne["powt"]+dt, 1):
        u[i]=u1[i-dt]
    y=np.zeros(zmienne["powt"])
    for i in range(0, zmienne["powt"],1):
        y[i]=u[i]
    return y

def czas():
    t = np.ones(zmienne["powt"])
    t = t * zmienne["krok"]
    t = np.cumsum(t)
    t -= zmienne["krok"]
    return t

def sygnalP():
    r = np.zeros(zmienne["powt"])
    pom=0
    pom2=1
    for j in range(0, zmienne["powt"], 1):
        if pom2==1:
            r[j]=zmienne["amplituda"]
        else:
            r[j]=-zmienne["amplituda"]
        pom+=1
        if pom==zmienne["polokres"]:
            pom=0
            pom2*=-1

    for i in range(zmienne["dlugoscSyg"],zmienne["powt"],1):
        r[i]=0
    return r

def sygnalT():
    r = np.zeros(zmienne["powt"])
    pom=zmienne["powt"]/zmienne["polokres"]
    pom=round(pom)
    pom2=1
    for i in range(0, pom, 1):
        for j in range(i*zmienne["polokres"], (i+1)*zmienne["polokres"], 1):
            if pom2==1:
                r[j]=r[j-1]+4*(zmienne["amplituda"]/zmienne["polokres"])
            else:
                r[j]=r[j-1]-4*(zmienne["amplituda"]/zmienne["polokres"])
            if r[j]>=zmienne["amplituda"] or r[j]<=-zmienne["amplituda"]:
                pom2*=-1
    for i in range(zmienne["dlugoscSyg"],zmienne["powt"],1):
        r[i]=0
    return r

def sygnalH(t):
    r=np.zeros(zmienne["powt"])
    for i in range(zmienne["powt"]):
        r[i]=zmienne["amplituda"] * np.sin(t[i]*np.pi*0.1*zmienne["powt"]/zmienne["polokres"])
    return r

def rysowanie(czestotliwosc, moduly, fazy, t, y, r):
    fig, wyk = wkrs.subplots(4, 1)
    wyk[0].set_xlabel("t")
    wyk[0].set_ylabel("u(t)")
    wyk[0].yaxis.set_label_coords(-0.07, 0.5)
    wyk[0].xaxis.set_label_coords(1.005, -0.025)
    wyk[1].set_xlabel("t")
    wyk[1].set_ylabel("y(t)")
    wyk[1].xaxis.set_label_coords(1.005, -0.025)
    wyk[1].yaxis.set_label_coords(-0.07, 0.5)
    wyk[2].set_xlabel("w")
    wyk[2].set_ylabel("stopnie")
    wyk[2].set_xscale("log")
    wyk[2].yaxis.set_label_coords(-0.1, 0.5)
    wyk[2].xaxis.set_label_coords(1.005, -0.025)
    wyk[3].set_xlabel("w")
    wyk[3].set_ylabel("dB")
    wyk[3].set_xscale("log")
    wyk[3].xaxis.set_label_coords(1.005, -0.025)
    wyk[3].yaxis.set_label_coords(-0.07, 0.5)
    wyk[3].plot(czestotliwosc, moduly)
    wyk[2].plot(czestotliwosc, np.unwrap(2 * fazy) / 2)
    wyk[1].plot(t, y)
    wyk[0].plot(t, r)
    wkrs.show()

if __name__ == "__main__":
    main(zmienne)
