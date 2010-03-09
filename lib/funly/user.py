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

from milia.metrics import Flrw
from milia.lumfuncs import Schechter

from funly import Population, Conditions, Filter, Survey, compute

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
