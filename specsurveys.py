import logging
import os
from matplotlib import pyplot as plt
import numpy as np

from astropy.cosmology import FlatLambdaCDM

from classy import Class


class_params = {
    'output': 'mPk',
    'non linear': 'halofit',
    'z_max_pk': 5.,
    'A_s': 2.3e-9,
    'n_s': 0.9624, 
    'h': 0.67,
    'Omega_b': 0.05,
    'Omega_cdm': 0.25}

class SpecSurvey(object):

    FULLSKY = 4 * np.pi * (180 / np.pi)**2

    def __init__(self):
        """ """
        self.data = {}
        omm = class_params['Omega_b'] + class_params['Omega_cdm']
        logging.info("Omega_m = %f"%omm)
        self.dist = FlatLambdaCDM(100, omm)
        self.cosmo = Class()
        self.cosmo.set(class_params)
        self.cosmo.compute()


    def loadSurveyFile(self, filename, dir="zdist"):
        """ """
        path = os.path.join(dir, filename)

        meta = {}
        zdist = []
        for line in file(path):
            line = line.strip()
            if ":" in line:
                key, value = line.split(":")
                key = key.strip()
                value = value.strip()
                try:
                    value = float(value)
                except ValueError:
                    pass
                meta[key] = value
                continue
            try:
                vec = [float(v) for v in line.split()]
                zdist.append(vec)
            except:
                raise
        meta['nz'] = zdist
        try:
            name = meta['name']
        except:
            name = len(data)
        self.data[name] = meta
        return meta

    def volume(self, z0, z1, area):
        """ """
        r0 = self.dist.comoving_distance(z0).value
        r1 = self.dist.comoving_distance(z1).value
        skyfrac = area / self.FULLSKY
        volume = 4. / 3 * np.pi * (r1**3 - r0**3) * skyfrac
        return volume

    def plotSurvey(self, data, mode='veff', lw=2, fontsize=18):
        """ """
        try:
            name = data['name']
            area = data['area']
            units = data['units'].lower()
        except KeyError:
            return
        print "read",name,area,units
        nz = data['nz']
        bias = 1.5
        skyfrac = area / self.FULLSKY
        x = []
        y = []
        for vec in nz:
            z0,z1,count = vec[:3]
            if len(vec) > 3:
                bias = vec[3]
            volume = self.volume(z0, z1, area)

            # print z0, z1, r0, r1, volume
            if units.startswith('c'):
                # print "got counts"
                density = count / volume
            elif units.startswith('s'):
                # print "Got sqr deg"
                density = count * area / volume * (z1-z0)
            elif units.startswith('m'):
                density = count
            else:
                print "unknown units!", units
                continue

            P0 = self.cosmo.pk(0.2, 0.5)#(z0+z1)/2.)

            # print name, z0,z1,P0,density*P0*bias**2, volume/1e9
            if mode == 'veff':
                fkp = 1./(1 + 1./(density * P0*bias**2))
                quant = fkp**2 * volume / 1e9
            elif mode == 'sn':
                quant = P0*bias**2 * density

            x.append(z0)
            x.append(z1)
            y.append(quant)
            y.append(quant)

        line, = plt.semilogy(x, y, lw=lw)
        c = line.get_color()
        plt.text(x[-1], y[-1], name, color=c, fontsize=fontsize, family='sans-serif', verticalalignment='center', zorder=10)

    def plot(self, mode='veff', **params):
        """ """
        for survey in self.data.values():
            self.plotSurvey(survey, mode=mode, **params)

    def plotCountsSurvey(self, survey, fontsize=12, **plotparams):
        """ """
        try:
            survey = self.loadSurveyFile(survey)
        except:
            raise
        # print survey

        z0, z1 = [float(v) for v in survey['zrange'].split()]
        count = survey['count']
        area = survey['area']
        volume = self.volume(z0, z1, area)

        quant = count / volume
        points = plt.scatter(volume, quant, **plotparams)
        c = points.get_facecolor()
        plt.text(volume*1.1, quant, survey['name'], fontsize=fontsize, color=c[0],
            family='sans-serif', verticalalignment='center', zorder=10)

    def plotCounts(self, **params):
        for survey in self.data.values():
            self.plotCountsSurvey(survey, **params)
