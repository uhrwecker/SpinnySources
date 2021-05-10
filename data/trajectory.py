import pandas as pd
import configparser as cp

class TrajectoryResults:
    def __init__(self, path, top_or_bottom='top'):
        self.path = path
        self.tob = top_or_bottom

        self.sigma = None
        self.t = None
        self.dt = None
        self.r = None
        self.dr = None
        self.theta = None
        self.dtheta = None
        self.phi = None
        self.dphi = None

        self.get_from_csv(path)

        self.s = None
        self.vr = None
        self.vphi = None
        self.gamma = None
        self.Eem = None
        self.Jem = None
        self.alpha = None
        self.beta = None

        self.r_obs = None
        self.theta_obs = None
        self.phi_obs = None

        self.Eph = None
        self.Lph = None
        self.Qph = None

        self.r0 = None
        self.phi0 = None
        self.theta0 = None

        self.get_from_ini(path)

    def get_from_csv(self, path):
        df = pd.read_csv(path+'.csv')

        self.sigma = df['sigma'].to_numpy()
        self.t = df['t'].to_numpy()
        self.dt = df['dt'].to_numpy()
        self.r = df['r'].to_numpy()
        self.dr = df['dr'].to_numpy()
        self.theta = df['theta'].to_numpy()
        self.dtheta = df['dtheta'].to_numpy()
        self.phi = df['phi'].to_numpy()
        self.dphi = df['dphi'].to_numpy()

    def get_from_ini(self, path):
        config = cp.ConfigParser()
        config.read(path+'.ini')

        self.s = float(config['Emitter']['s'])
        self.vr = float(config['Emitter']['v_r'])
        self.vphi = float(config['Emitter']['v_phi'])
        self.gamma = float(config['Emitter']['gamma'])
        self.Eem = float(config['Emitter']['e_emitter'])
        self.Jem = float(config['Emitter']['j_emitter'])
        self.alpha = float(config['Emitter']['alpha'])
        self.beta = float(config['Emitter']['beta'])

        self.r_obs = float(config['Observer']['r_obs'])
        self.theta_obs = float(config['Observer']['theta_obs'])
        self.phi_obs = float(config['Observer']['phi_obs'])

        self.Eph = float(config['Photon']['e_photon'])
        self.Lph = float(config['Photon']['l_photon'])
        self.Qph = float(config['Photon']['q_photon'])

        self.r0 = float(config['Photon']['r0'])
        self.theta0 = float(config['Photon']['theta0'])
        self.phi0 = float(config['Photon']['phi0'])


