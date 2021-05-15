class Material():
    E=0
    nu=0
    mu=0
    lam=0
    rho=0
    type="linear"
    def __init__(self, type="linear", E=0, nu=0,rho=0):
        self.E=E
        self.nu=nu
        self.type=type
        if(abs(E)>10e-15 and abs(nu)>10e-15 and abs(nu+1)>10e-15 and abs(nu-1/2)>10e-15 and abs(rho)>10e-15):
            self.mu=E / 2 / (1+nu)
            self.lam=E * nu / ((1+nu)*(1-2*nu))
            self.rho=rho


#definging default materials:
#type=Material(type of model (e.g. "linear"), E (=young-module),nu (=poisson number), rho (=density))
aluminium=Material("linear", 69 * (10**(6)), 0.33, 2700 * (10**(-9)))
aluminum=Material("linear", 69 * (10**(6)), 0.33, 2700 * (10**(-9)))
steel=Material("linear", 220 * (10**(6)), 0.28, 7700 * (10**(-9)))