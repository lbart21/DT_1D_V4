import math as m
class FluidModel():
    def __init__(self, fluidDefinition) -> None:
        self.fs = {}
        self.fluidDefinition = fluidDefinition
        pass

    def UpdateFromRhoP(self):
        p = self.fs["p"]
        rho = self.fs["rho"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        T = p / (R * rho)
        u = Cv * T
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        s = s_ref + Cv * m.log(T / T_ref) + R * m.log((T * p_ref) / (p * T_ref))

        self.fs["T"] = T
        self.fs["u"] = u
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma


    def UpdateFromPT(self):
        p = self.fs["p"]
        T = self.fs["T"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        rho = p / (R * T)
        u = Cv * T
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        s = s_ref + Cv * m.log(T / T_ref) + R * m.log((T * p_ref) / (p * T_ref))

        self.fs["rho"] = rho
        self.fs["u"] = u
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma


    def UpdateFromRhoU(self):
        rho = self.fs["rho"]
        u = self.fs["u"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        p = (gamma - 1.0) * rho * u
        T = u / Cv
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        s = s_ref + Cv * m.log(T / T_ref) + R * m.log((T * p_ref) / (p * T_ref))

        self.fs["p"] = p
        self.fs["T"] = T
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma


    def UpdateFromPU(self):
        p = self.fs["p"]
        u = self.fs["u"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        rho = p / ((gamma - 1.0) * u)
        T = u / Cv
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        s = s_ref + Cv * m.log(T / T_ref) + R * m.log((T * p_ref) / (p * T_ref))

        self.fs["rho"] = rho
        self.fs["T"] = T
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma


    def UpdateFromRhoT(self):
        rho = self.fs["rho"]
        T = self.fs["T"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        p = rho * R * T
        u = Cv * T
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        s = s_ref + Cv * m.log(T / T_ref) + R * m.log((T * p_ref) / (p * T_ref))

        self.fs["p"] = p
        self.fs["u"] = u
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma

    def UpdateFromPS(self):
        p = self.fs["p"]
        s = self.fs["s"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        T = T_ref * ((p / p_ref) ** R * m.exp(s - s_ref)) ** (1.0 / Cp)
        rho = p / (R * T)
        u = Cv * T
        h = u + p / rho
        a = (gamma * R * T) ** 0.5
        

        self.fs["rho"] = rho
        self.fs["u"] = u
        self.fs["h"] = h
        self.fs["a"] = a
        self.fs["s"] = s
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma

    def UpdateFromHS(self):
        h = self.fs["h"]
        s = self.fs["s"]

        R = self.fluidDefinition["R"]
        Cv = self.fluidDefinition["Cv"]
        Cp = self.fluidDefinition["Cp"]
        gamma = self.fluidDefinition["gamma"]
        s_ref = self.fluidDefinition["s_ref"]
        T_ref = self.fluidDefinition["T_ref"]
        p_ref = self.fluidDefinition["p_ref"]

        T = h / Cp 
        p = p_ref * m.exp(1.0 / R * (s_ref - s + Cp * m.log(T / T_ref)))
        rho = p / (R * T)
        u = Cv * T
        a = (gamma * R * T) ** 0.5


        self.fs["rho"] = rho
        self.fs["u"] = u
        self.fs["p"] = p
        self.fs["a"] = a
        self.fs["T"] = T
        self.fs["Cv"] = Cv
        self.fs["Cp"] = Cp
        self.fs["gamma"] = gamma

    def addMachNumber(self):
        self.fs["Ma"] = self.fs["vel_x"] / self.fs["a"]