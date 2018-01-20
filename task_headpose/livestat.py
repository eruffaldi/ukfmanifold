#
# Python Livestat module
#
# This module provides a simple mechanism for computing statistics of variables as they are produced.
# In particular: count, mean, std, maximum and minimum, span
#
# Two statistics can be merged together, for example in a multiprocessing context
#
# Future: numpy arrays, tuples/lists, so far only scalar
#
# For online algorithm source:
# http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
#
# See also the incmoments module for the same core operations but performed over tuples
# See also the current work on Statistics inside Python: http://www.python.org/dev/peps/pep-0450/   
#
# Last Updated: 31 December 2013
#
# Emanuele Ruffaldi 2012-2014

from collections import defaultdict
import math
try:
    import numpy
except:
    numpy = None


class LiveStat:
    """ LiveStat allows to compute statistics over variables as they are produced"""
    def __init__(self,name=""):
        """Constructor with optional name, used for printing"""
        self.name = name
        self.dirty = False
        self.reset()
    @property
    def empty(self):
        """Returns true when there is no data"""
        return self.vcount == 0
    @property
    def count(self):
        """Returns the number of items seen by the accumulator"""
        return self.vcount
    @property
    def sum(self):
        """Returns the sum of the values. None if no items"""
        return self.vsum
    @property
    def kurtosis(self):
        """Returns the kurtosis of the distribution, that is the flatness. Zero for Gaussians"""
        self._finalize()
        return self.vkurtosis
    @property
    def skewness(self):
        """Returns the skewness as the asymmetri before/after the mean. Zero for Gaussians"""
        self._finalize()
        return self.vskewness
    @property
    def jarquebetatest(self,alpha=0.05):
        """Returns the Jarque-Beta test of normality based on kurtosis and skewness"""
        self._finalize()
        JB = self.vcount/6*(self.vskewness**2 + 1/4*((self.vkurtosis-2)**2))
        #the JB statistic asymptotically has a chi-squared distribution with two degrees of freedom, so the statistic can be used to test the hypothesis that the data are from a normal distribution.
        """If the data comes from a normal distribution, the JB statistic asymptotically has a chi-squared distribution with two degrees of freedom, so the statistic can be used to test the hypothesis that the data are from a normal distribution. The null hypothesis is a joint hypothesis of the skewness being zero and the excess kurtosis being zero. Samples from a normal distribution have an expected skewness of 0 and an expected excess kurtosis of 0 (which is the same as a kurtosis of 3). As the definition of JB shows, any deviation from this increases the JB statistic."""
        raise Exception("Not implemented")
    @property
    def mean(self):
        """Returns the sample mean of the values. None if no items"""
        return self.vmean
    @property
    def span(self):
        """Returns the range of values. None if no items"""
        if self.vcount == 0:
            return None
        else:
            return self.vmax-self.vmin

    @property
    def variance(self):
        """Returns the sample variance of the values. None if no items"""
        if self.dirty:
            self._finalize()
        return self.vvar

    @property
    def std(self):
        """Returns the sample standard deviation of the values. None if no items"""
        if self.dirty:
            self._finalize()
        if self.vvar is None:
            return 0
        else:
            return math.sqrt(self.vvar)

    def reset(self):
        """Resets the accumulator"""
        self.vmin = None
        self.vmax = None
        self.vmean = None
        self.vsum = None
        self.vm2 = None
        self.vm3 = None
        self.vm4 = None

        self.vcount = 0
        self.vcount2 = 0

        self.dirty = False
        self.vvar = None
        self.vcurtosis = None
        self.vskewness = None
    def _finalize(self):
        """Private Finalize the interal counter to update the variance"""
        if self.vcount > 1:
            # skewness = g1 = sqrt(n) M3/(M2^(3/2)) # zero 
            # kurtosis = g2 = n M4/M2^2 - 3 # zero for normal
            #    sk = (M3/nf)/(sigma**3)
            #    ku = (M4/nf)/sigma**4 - 3
            n = self.vcount
            nf = float(n)
            mu2 = self.vm2/nf
            self.vvar = self.vm2/(nf-1)
            try:
                self.vskewness = self.vm3/nf/(mu2**1.5)
                self.vkurtosis = self.vm4/nf/(mu2**2)
            except:
                self.vskewness = 0
                self.vkurtosis = 0
        elif self.vcount == 1:
            self.vvar = 0
            self.vskewness = 0
            self.vkurtosis = 0
        self.dirty = False

    def __imul__(self,value):
        """Updates the statistics as if all the input values were (x*value)
        """
        if isinstance(value,LiveStat):
            raise Exception ("Product of Statistics is not supported")
        else:
            if self.vmin is not None:
                # mu(s x) = 1/N sum s x = s/N sum x
                self.vmean *= value
                if value < 0:
                    m = self.vmin
                    M = self.vmax
                    self.vmin = M*value
                    self.vmax = m*value
                else:
                    self.vmin *= value
                    self.vmax *= value
                self.vsum *= value
                self.vm2 *= value*value
                # no effect on vm4
                self.dirty = True
        return self

    def __idiv__(self,value):
        """Updates the statistics as if all the input values were (x/value)"""
        if isinstance(value,LiveStat):
            raise Exception ("Ratio of Statistics is not supported")
        else:
            if self.vmin is not None:
                # mu(s x) = 1/N sum s x = s/N sum x
                self.vmean /= value
                if value < 0:
                    m = self.vmin
                    M = self.vmax
                    self.vmin = M/value
                    self.vmax = m/value
                else:
                    self.vmin /= value
                    self.vmax /= value
                self.vsum /= value
                # vm2(s x) = sum (s x - mu(s x))^2 = sum (s x - s mu(x))^2 = sum s^2 (x - mu(x))^2 = s^2 sum (x - mu(x))^2 = s^2 vm^2
                self.vm2 /= value*value
                self.dirty = True
        return self
    def extend(self,data):        
        """Extend from sequence"""
        n = float(len(data))
        if n == 0:
            return self
        M2 = 0
        M3 = 0
        M4 = 0
        mean = 0
        vmin = None
        vmax = None
        for x in data:
            mean += x/n   
            if vmin is None:
                vmax = x
                vmin = x
            if x < vmin:
                vmin = x
            if x > vmax:
                vmax = x
        for x in data:
            d = x-mean
            M2 += (d**2)
            M3 += (d**3)
            M4 += (d**4)
        x = LiveStat(self.name)
        x.vmin = vmin
        x.vmax = vmax
        x.vmean = mean
        x.vm2 = M2
        x.vm3 = M3
        x.vm4 = M4
        x.vcount = int(n)
        x.vcount2 = x.vcount**2
        x.dirty = True
        self.merge(x)
        return self
    def __add__(self,value):
        """Addition operator: scalar applied to all terms x_i"""
        x = self.clone()
        if isinstance(value,LiveStat):
            x.name = "(" + self.name + "+" + value.name + ")"
        else:
            x.name = "(" + self.name + "+ scalar)"
        x += value
        return x
    def __sub__(self,value):
        """Subtraction operator: scalar applied to all terms x_i"""
        x = self.clone()
        if isinstance(value,LiveStat):
            x.name = "(" + self.name + "-" + value.name + ")"
        else:
            x.name = "(" + self.name + "- scalar)"
        x -= value
        return x
    def __mul__(self,value):
        """Scales all terms by value"""
        x = self.clone()
        if isinstance(value,LiveStat):
            x.name = "(" + self.name + "*" + value.name + ")"
        else:
            x.name = "(" + self.name + "* scalar)"
        x *= value
        return x
    def __div__(self,value):
        """Divides all terms by value"""
        x = self.clone()
        if isinstance(value,LiveStat):
            x.name = "(" + self.name + "/" + value.name + ")"
        else:
            x.name = "(" + self.name + "/ scalar)"
        x /= value
        return x    
    def __iadd__(self,value):
        """Updates the statistics as if all the values were (x+scalar) or (x+value)"""
        if isinstance(value,LiveStat):
            raise Exception("Cannot sum statistics")
            if value.vcount < 1 or self.vcount < 1:
                raise Exception("Cannot sum empty statistics")
            else:
                # sum of two considered pairwise: z_i = stat(x_i + y_i)
                #
                # data have different weights due to number of samples.. TODO
                self.vmin += value.vmin 
                self.vmax += value.vmax
                self.vmean += value.vmean
                self.vsum += value.vsum
                # variance is sum of variance?
                self.vm2 += value.vm2
                # TODO vm3 vm4
                self.vcount = min(value.vcount,self.vcount)
                self.vcount2 = self.vcount**2
                self.dirty = True
        else:
            # constant bias
            if self.vmin is not None:
                self.vmin += value
                self.vmax += value
                self.vmean += value
                self.vsum += self.vcount*value
                self.dirty = True
        return self
    def __isub__(self,value):
        """Updates the statistics as if all the values were (x-value) and (x-y)"""
        if isinstance(value,LiveStat):
            raise Exception("Cannot sum statistics")
            if value.vcount < 1 or self.vcount < 1:
                raise Exception("Cannot sum empty statistics")
            else:
                # sum of two considered pairwise: z_i = stat(x_i - y_i)
                #
                # data have different weights due to number of samples.. TODO
                self.vmin = self.vmin-value.vmax
                self.vmax = self.vmax-value.vmin
                self.vmean -= value.vmean
                self.vsum -= value.vsum
                # variance is sum of variance in any case
                self.vm2 += value.vm2
                # TODO vm3 vm4
                self.vcount = min(self.vcount,value.vcount)
                self.vcount2 = self.vcount**2
                self.dirty = True
        else:
            # constant bias
            if self.vmin is not None:
                self.vmin -= value
                self.vmax -= value
                self.vmean -= value
                self.vsum -= self.vcount*value
                self.dirty = True
        return self    
    def clone(self):
        r = LiveStat(self.name)
        r.copy(self)
        return r
    def copy(self,other):
        """Assignment"""
        if other.vcount == 0:
            self.reset()
        else:
            self.vcount = other.vcount
            self.vcount2 = other.vcount2
            self.vmin = other.vmin
            self.vmax = other.vmax
            self.vm2 = other.vm2
            self.vm3 = other.vm3
            self.vm4 = other.vm4
            self.vmean = other.vmean
            self.vsum = other.vsum
            self.name = other.name
            self.dirty = True
        return self
    def append(self,x):
        """Appends a new item"""
        if self.empty:
            self.vcount = 1
            self.vcount2 = 1
            self.vmin = x
            self.vmax = x
            self.vsum = x
            self.vmean = x
            self.vm2 = 0
            self.vm3 = 0
            self.vm4 = 0
            self.dirty = True
        else:
            nA = self.vcount
            nAA = self.vcount2
            nX = nA+1
            nXX = nX**2
            nXXX = nX**3
            self.vcount = nX
            self.vcount2 = nXX

            if x < self.vmin:
                self.vmin = x
            if x > self.vmax:
                self.vmax = x

            delta = x-self.vmean
            delta2 = delta**2
            delta3 = delta**3
            delta4 = delta**4
            self.vmean += delta/nX # incremental mean (good for vectorial)
            self.vm3 += delta3*(nA*(nA-1))/nXX - 3*delta*self.vm2/nX
            self.vm4 += delta4*(nA*(nAA-nA+1))/nXXX + 6*delta2*(self.vm2)/nXX - 4*delta*self.vm3/nX
            # note is done at end
            self.vm2 += delta2*nA/nX # incremental quadratic for variance (good for vectorial)
            self.vsum += x

            self.dirty = True
    def merge(self,other):
        """Merges the current statistics with the other"""
        if self.empty:            
            self.copy(other)
            return self
        elif other.empty:
            return self
        if(other.vmin < self.vmin):
            self.vmin = other.vmin
        if(other.vmax > self.vmax):
            self.vmax = other.vmax

        nA = float(self.vcount)
        nB = float(other.vcount)
        nAB = nA*nB
        nAA = float(self.vcount2)
        nBB = float(other.vcount2)
        nX = nA+nB
        nXX = nX**2 #nAA+nBB+2*nAB #nX**2 # actually (nA+nB)^2 = (nAA+nBB+2*nAB)
        nXXX = nXX*nX
        self.vcount = nX
        self.vcount2 = nXX

        self.vsum += other.vsum;

        # merge of mean and m2
        delta = other.vmean-self.vmean;
        delta2 = delta**2
        delta3 = delta**3
        delta4 = delta**4
        self.vmean += delta*nB/nA
        self.vm2 += other.vm2 + delta2*(nAB/nX)
        self.vm3 += other.vm3 + delta3*(nAB*(nA-nB))/nXX + 3*delta*(nA*other.vm2-nB*self.vm2)/nX
        self.vm4 += other.vm4 + delta4*(nAB*(nAA-nAB+nBB))/nXXX + 6*delta2*(nAA*other.vm2+nBB*self.vm2)/nXX + 4*delta*(nA*other.vm3-nB*self.vm3)/nX
        self.dirty = True
        return self
    def asdict(self):
        prefix = self.name
        self._finalize()
        return dict([(prefix+"_mean",self.vmean),(prefix+"_min",self.vmin),(prefix+"_max",self.vmax),(prefix+"_std",self.std),(prefix+"_skew",self.vskewness),(prefix+"_count",self.count),(prefix+"_sum",self.vsum)])
    def __str__(self):
        """String representation"""
        self._finalize()
        np = self.name
        if self.name != "":
            np += ","
        if self.vcount > 0:
            return "LiveStat(%smean=%s,std=%s,min=%s,max=%s,skew=%s,kurt=%s,count=%d)" % (np,self.vmean,self.std,self.vmin,self.vmax,self.skewness,self.kurtosis,self.vcount)
        else:
            return "LiveStat(%sempty)" % np

class DeltaLiveStat(LiveStat):
    """Specialization of the LiveStat that manages differential statistics"""
    def __init__(self,name=""):
        self.last = None
        self.dlast = None
        LiveStat.__init__(self,name)
    def reset(self):        
        """Reset"""
        self.last = None
        LiveStat.reset(self)
    def resetlast(self):
        """Reset only the last, but not the inner statistics. Equivalent to adding None"""
        self.last = None
        self.dlast = 0
    def clone(self):
        r = DeltaLiveStat(self.name)
        r.copy(self)
        return r
    def copy(self,other):
        LiveStat.copy(self,other)
        self.last = other.last
        self.dlast = other.dlast
        return self
    def append(self,x):
        """Adds a new item. If x is None this means to reset the input"""
        if x is None:
            self.last = None
        elif self.last is None:
            self.last = x
            self.dlast = 0
        else:
            self.dlast = x-self.last
            LiveStat.append(self,float(self.dlast))
            self.last = x
    def __str__(self):
        self._finalize()
        np = self.name
        if np != "":
            np += ","
        if self.vcount > 0:
            return "DeltaLiveStat(%smean=%s,std=%s,min=%s,max=%s,count=%d)" % (np,self.vmean,self.std,self.vmin,self.vmax,self.vcount)
        else:
            return "DeltaLiveStat(%sempty)" % np


class Counter:
    """Simple counter class with interface similar to LiveStat"""
    def __init__(self):
        self.c = 0
    @property
    def count(self):
        return self.c
    @property
    def empty(self):
        return self.c == 0
    def append(self,x):
        self.c += x
    def __str__(self):
        return str(self.c)
    def _finalize(self):
        pass
    def divide(self,x):
        self.c /= x

class Histogram:
    """Histogram Class"""
    def __init__(self,name,cz=Counter):
        self.name = name
        self.cz = cz
        self.items = None
        self.count = 0
        self.reset()
    @property
    def empty(self):
        return self.count == 0
    @property
    def casescount(self):
        return len(self.items)
    @property
    def cases(self):
        return sorted(self.items.keys())
    def normalizetotal(self):
        n = self.count
        for v in self.items():
            v.divide(n)
    @property
    def count(self):
        return self.count
    def append(self,x,y=1):
        self.items[x].append(y)
        self.count += 1
    def reset(self):
        self.items = defaultdict(self.cz)
        self.count = 0
    def _finalize(self):
        pass






