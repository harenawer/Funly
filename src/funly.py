#
# Copyright 2010 Sergio Pascual
# 
# This file is part of Funly
# 
# Funly is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Funly is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Funly.  If not, see <http://www.gnu.org/licenses/>.
#

import math
import readline

import numpy
from milia.metrics import Flrw
from milia.lumfuncs import Schechter

LOGFCGS_LCGS = 50.07791096
SOLID_ANGLE_TO_DEG= (180 / math.pi)**2

class Conditions(object):
    def __init__(self, metric, populations):
        self.metric = metric
        self.populations = populations

class Population(object):
    def __init__(self, line, lumfunc, extinction):
        self.lumfunc = lumfunc
        self.line = line
        self.extinction = extinction

class Filter(object):
    def __init__(self, cwl, fwhm):
        self.cwl = cwl
        self.fwhm = fwhm

class Survey(object):
    def __init__(self, f1, f2, filter_, area):
        self.f1 = f1
        self.f2 = f2
        self.filter_ = filter_
        self.area = area

    def redshift_limits(self, line):
        z1 = (self.filter_.cwl - self.filter_.fwhm / 2.0) / line - 1
        z2 = (self.filter_.cwl + self.filter_.fwhm / 2.0) / line - 1
        return z1, z2

class RedshiftSurvey(Survey):
    def __init__(self, f1, f2, z1, z2, filter_, area):
        super(RedshiftSurvey, self).__init__(f1, f2, filter_, area)
        self.z1 = z1
        self.z2 = z2

    def redshift_limits(self, line):
        return self.z1, self.z2

def compute(conditions, survey, delta_z=1e-4, delta_m=1.0):
    delta_lf = 0.4 * delta_m
    delta_f = pow(10, delta_lf)

    lfluxes = numpy.arange(math.log10(survey.f1), numpy.log10(survey.f2), delta_lf)
   
    bins = numpy.zeros((len(conditions.populations), len(lfluxes)))

    metric = conditions.metric

    for pind, pop in enumerate(conditions.populations):
        area = survey.area / SOLID_ANGLE_TO_DEG
        z1, z2 = survey.redshift_limits(pop.line)
        for z in numpy.arange(z1, z2, delta_z):
            llcon = LOGFCGS_LCGS + 2 * math.log10(metric.dl(z)) + 0.4 * pop.extinction
            vol = area * (-metric.vol(z) + metric.vol(z + delta_z))
	    pop.lumfunc.evolve(z)
            for lind, lflux in enumerate(lfluxes):
                llum1 = llcon + lflux
                llum2 = llum1 + delta_lf
                bins[pind, lind] += vol * pop.lumfunc.object_density(pow(10, llum1), pow(10, llum2))

    f1 = numpy.power(10, lfluxes)
    f2 = numpy.power(10, lfluxes + delta_lf)
    result = numpy.hstack ([f1[:,numpy.newaxis], f2[:,numpy.newaxis], bins.T])
    return result


def main():
    # Conditions
    # metric
    metric = Flrw(70, 0.3, 0.7)

    line = 1215.668
    logphi = float(raw_input('log Phi*: '))
    loglum = float(raw_input('log L*: '))
    alpha = float(raw_input('alpha: '))
    # lumfunc
    lumfunc = Schechter(10**logphi, 10**loglum, alpha)
    # extinction
    extinction = float(raw_input('Average extinction (magnitudes): '))
    population1 = Population(line, lumfunc, extinction)
    conditions = Conditions(metric, [population1])

    # Survey
    # flux limit
    f1 = float(raw_input('Lower flux limit: '))
    f2 = float(raw_input('Upper flux limit: '))

    # redshift limit
    cwl = float(raw_input('CWL of the filter: '))
    fwhm = float(raw_input('FWHM of the filter: '))
    mfilter = Filter(cwl, fwhm)

    # area
    area = float(raw_input('Surveyed area (sq arc minutes): '))
    area *= 1. / 3600 # square def

    survey = Survey(f1, f2, mfilter, area)

    result = compute(conditions, survey)

    for i in result:
        print "%09.3g - %09.3g" % (i[0], i[1]),
        for n in i[2:]:
            print '%.2f' % n,
        print

    print 'Totals: ', result[:,2:].sum(0)
   

def main2():
    # Universe
    # metric
    metric = Flrw(70, 0.3, 0.7)

    line = 1215.668
    # lumfunc
    lumfunc = Schechter(10**-2.88, 10**42.6, -1.5)
    # extinction
    extinction = 0
    population1 = Population(line, lumfunc, extinction)
    line = 1215.668
    # lumfunc
    lumfunc = Schechter(10**-2.88, 10**42.6, -1.5)
    # extinction
    extinction = 0
    population2 = Population(line, lumfunc, extinction)
    conditions = Conditions(metric, [population1, population2])
    # Survey
    # flux limit
    f1 = 1e-18
    f2 = 1e-14

    # redshift limit, given by filter
    cwl = 9200
    fwhm = 12
    mfilter = Filter(cwl, fwhm)

    # area
    area = 50.26
    area *= 1. / 3600 # square def

    survey = Survey(f1, f2, mfilter, area)

    result = compute(conditions, survey)

    for i in result:
        print "%9.3g - %9.3g" % (i[0], i[1]),
        for n in i[2:]:
            print '%.2f' % n,
        print

    print 'Totals: ', result[:,2:].sum(0)



if __name__ == "__main__":
    main()
