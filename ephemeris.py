class Ephemeris:
    def __init__(self, time, e, M0, deltaN, sqrtA, omega, Crs, Crc, Cuc, Cus, Cic, Cis, IDOT, Io, omega0, omegaDot, Toe):
        self.M0 = M0
        self.time = time
        self.e = e
        self.deltaN = deltaN
        self.sqrtA = sqrtA
        self.omega = omega
        self.Crs = Crs
        self.Crc = Crc
        self.Cuc = Cuc
        self.Cus = Cus
        self.Cic = Cic
        self.Cis = Cis
        self.IDOT = IDOT
        self.Io = Io
        self.omega0 = omega0
        self.omegaDot = omegaDot
        self.Toe = Toe
