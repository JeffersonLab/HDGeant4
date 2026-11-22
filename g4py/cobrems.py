#
# cobrems.py : python class for computing coherent bremsstrahlung spectra
#              and rates, and for generating Monte Carlo samples of them.
#
# usage example:
#   $ python
#   >>> from ROOT import *
#   >>> import cobrems
#   >>> cobrems.init(E0=12, Epeak=8.7)
#   >>> cobrems.plotTotal(0)
#
# author: richard.t.jones at uconn.edu
# version: july 29, 2015
#

import sys
from array import array

import ROOT
from Cobrems import *

# set global variables to defaults
Epeak = 9
cur = 2.2
E0 = 12
Erms = 0.0006
emit = 2.5e-9
spot = 0.5e-3
dist = 76
coldiam = 0.0034
radt = 20e-6
nbins = 1000
Emin = 0
Emax = E0
pElim0 = 8.4
pElim1 = 9.0
bElim0 = 1.02
bElim1 = 2.04
eElim0 = 10.68
eElim1 = 11.7

# initialize the built-in generator to defaults
generator = CobremsGeneration(E0, Epeak)

def usage():
   """
   Print a brief summary of the internal configuration variables
   that regulate the behavior of the CobremsGeneration and exit.
   """
   print("""
   Usage: cobrems.init(key=value, ...)
     where key should be one of the following, [defaults in backets]:
      * Epeak - photon energy (GeV) of primary coherent edge [9]
      * cur - current (uA) of electron beam [2.2]
      * E0 - mean energy (GeV) of electron beam [12]
      * Erms - rms energy width (GeV) of electron beam [0.0006]
      * emit - transverse emittance (m rad) of electron beam [2.5e-9]
      * spot - virtual e-beam spot rms (m) at primary collimator [0.5e-3]
      * dist - distance (m) from radiator to primary collimator [76]
      * coldiam - diameter (m) of primary collimator [0.0034]
      * radt - thickness (m) of radiator crystal [20e-6]
      * nbins - number of bins to use for plotting spectra [200]
      * Emin - low energy limit (GeV) for spectrum historgrams [0]
      * Emax - high energy limit (GeV) for spectrum historgrams [12]
      * pElim0 - lower limit (GeV) of primary peak integration window [8.4]
      * pElim1 - upper limit (GeV) of primary peak integration window [9.0]
      * bElim0 - lower limit (GeV) of beam background integration window [1.02]
      * bElim1 - upper limit (GeV) of beam background integration window [2.04]
      * eElim0 - lower limit (GeV) of beam endpoint integration window [10.68]
      * eElim1 - upper limit (GeV) of beam endpoint integration window [11.70]
   """)

def init(**kwargs):
   """
   Loads a new set of beamline parameters and generates the acceptance
   function and collimated beam intensity spectra as ROOT TF1/TH1 objects.
   To see a list of supported options as keyword arguments, call init as
    >>> cobrems.init(o='help')
   """
   for arg, value in kwargs.iteritems():
      if arg == "E0":
         global E0
         E0 = value
         generator.setBeamEnergy(E0)
      elif arg == "Epeak":
         global Epeak
         Epeak = value
         generator.setCoherentEdge(Epeak)
      elif arg == "cur":
         global cur
         cur = value
      elif arg == "Erms":
         global Erms
         Erms = value
         generator.setBeamErms(Erms)
      elif arg == "emit":
         global emit
         emit = value
         generator.setBeamEmittance(emit)
      elif arg == "spot":
         global spot
         spot = value
         generator.setCollimatorSpotrms(spot)
      elif arg == "dist":
         global dist
         dist = value
         generator.setCollimatorDistance(dist)
      elif arg == "coldiam":
         global coldiam
         coldiam = value
         generator.setCollimatorDiameter(coldiam)
      elif arg == "radt":
         global radt
         radt = value
         generator.setTargetThickness(radt)
      elif arg == "nbins":
         global nbins
         nbins = value
      elif arg == "Emin":
         global Emin
         Emin = value
      elif arg == "Emax":
         global Emax
         Emax = value
      elif arg == "pElim0":
         global pElim0
         pElim0 = value
      elif arg == "pElim1":
         global pElim1
         pElim1 = value
      elif arg == "bElim0":
         global bElim0
         bElim0 = value
      elif arg == "bElim1":
         global bElim1
         bElim1 = value
      elif arg == "eElim0":
         global eElim0
         eElim0 = value
      elif arg == "eElim1":
         global eElim1
         eElim1 = value
      else:
         usage()
         return

   generator.printBeamlineInfo()
   generator.setCollimatedFlag(True)
   generator.setPolarizedFlag(False)

   # plot the acceptance curve
   global acceptF1
   acceptF1 = ROOT.TF1("acceptF1", acceptance, 0, 2, 0)
   acceptF1.SetTitle("\\mbox{acceptance vs }\\theta/\\theta_C")

   # plot the beam energy spectrum
   x0 = Emin / E0
   x1 = Emax / E0
   global dRtdxF1
   dRtdxF1 = ROOT.TF1("dRtdxF1", dRtdx, x0, x1, 0)
   dRtdxF1.SetTitle("\\mbox{photon beam spectrum vs }E_\\gamma\\,/E_0" +
                    "\\mbox{ (/GeV/s)}")
   dRtdxF1.SetNpx(nbins)

   # apply the beam-crystal convolution
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRtdx([xvals[i]])
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRtdkH1
   dRtdkH1 = 0
   title = "\\mbox{photon beam spectrum vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   dRtdkH1 = ROOT.TH1D("dRtdkH1", title, nbins, Emin, Emax)
   dRtdkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRtdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRtdkH1.SetStats(0)
   dRtdkH1.Draw("hist")

   # print the rates in the three integration windows
   persec = (Emax - Emin) * 1./nbins
   prate = dRtdkH1.Integral(dRtdkH1.FindBin(pElim0),
                            dRtdkH1.FindBin(pElim1) - 1) * persec
   brate = dRtdkH1.Integral(dRtdkH1.FindBin(bElim0),
                            dRtdkH1.FindBin(bElim1) - 1) * persec
   erate = dRtdkH1.Integral(dRtdkH1.FindBin(eElim0),
                            dRtdkH1.FindBin(eElim1) - 1) * persec
   print("beam rate in the peak window [" + str(pElim0) + "," + \
         str(pElim1) + "] GeV = ", prate, "/s")
   print("beam rate in the background window [" + str(bElim0) + "," + \
         str(bElim1) + "] GeV =", brate, "/s")
   print("beam rate in the endpoint window [" + str(eElim0) + "," + \
         str(eElim1) + "] GeV =", erate, "/s")

def plotCoherent(collimated=1):
   """
   Plot the coherent part of the spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(False)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRcdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRcdkH1
   dRcdkH1 = 0
   if collimated:
      title = "collimated photon beam spectrum, coherent part only"
   else:
      title = "uncollimated photon beam spectrum, coherent part only"
   dRcdkH1 = ROOT.TH1D("dRcdkH1", title, nbins, Emin, Emax)
   for i in range(0, nbins):
      dRcdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRcdkH1.GetXaxis().SetTitle("E_{#gamma} (GeV)")
   dRcdkH1.GetYaxis().SetTitle("rate (counts/GeV/s)")
   dRcdkH1.GetYaxis().SetTitleOffset(1.5)
   dRcdkH1.SetStats(0)
   dRcdkH1.Draw("hist")
   return dRcdkH1

def plotEnhancement(collimated=1):
   """
   Plot the ratio of the total diamond spectrum to an amorphous
   bremsstrahlung spectrum with the same endpoint, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(False)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRcdx([xvals[i]]) / dRidx([xvals[i]]) + 1
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global enhanH1
   enhanH1 = 0
   if collimated:
      title = "collimated photon beam spectrum, enhancement"
   else:
      title = "uncollimated photon beam spectrum, enhancement"
   enhanH1 = ROOT.TH1D("enhanH1", title, nbins, Emin, Emax)
   for i in range(0, nbins):
      enhanH1.Fill(xvals[i] * E0, yvals[i])
   enhanH1.GetXaxis().SetTitle("E_{#gamma} (GeV)")
   enhanH1.GetYaxis().SetTitle("enhancement")
   enhanH1.GetYaxis().SetTitleOffset(1.5)
   enhanH1.SetStats(0)
   enhanH1.Draw("hist")
   return enhanH1

def plotIncoherent(collimated=1):
   """
   Plot the incoherent part of the spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(False)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRidx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRidkH1
   dRidkH1 = 0
   if collimated:
      title = "collimated photon beam spectrum, incoherent part only"
   else:
      title = "uncollimated photon beam spectrum, incoherent part only"
   dRidkH1 = ROOT.TH1D("dRidkH1", title, nbins, Emin, Emax)
   dRidkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRidkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRidkH1.GetXaxis().SetTitle("E_{#gamma} (GeV)")
   dRidkH1.GetYaxis().SetTitle("rate (counts/GeV/s)")
   dRidkH1.GetYaxis().SetTitleOffset(1.5)
   dRidkH1.SetStats(0)
   dRidkH1.Draw("hist")
   return dRidkH1

def plotTotal(collimated=1):
   """
   Plot the total coherent+incoherent spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(False)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRtdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRtdkH1
   dRtdkH1 = 0
   if collimated:
      title = "collimated photon beam spectrum"
   else:
      title = "uncollimated photon beam spectrum"
   dRtdkH1 = ROOT.TH1D("dRtdkH1", title, nbins, Emin, Emax)
   dRtdkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRtdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRtdkH1.GetXaxis().SetTitle("E_{#gamma} (GeV)")
   dRtdkH1.GetYaxis().SetTitle("rate (counts/GeV/s)")
   dRtdkH1.GetYaxis().SetTitleOffset(1.5)
   dRtdkH1.SetStats(0)
   dRtdkH1.Draw("hist")
   return dRtdkH1

def plotPolarization(collimated=1):
   """
   Plot the linear polarization spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   ypols = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      generator.setPolarizedFlag(True)
      ypols[i] = dRcdx([xvals[i]])
      generator.setPolarizedFlag(False)
      yvals[i] = dRtdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)
   generator.applyBeamCrystalConvolution(nbins, xvals, ypols)

   # convert the result to a spectrum vs photon energy
   global polarH1
   polarH1 = 0
   if collimated:
      title = "collimated photon beam polarization"
   else:
      title = "uncollimated photon beam polarization"
   polarH1 = ROOT.TH1D("polarH1", title, nbins, Emin, Emax)
   for i in range(0, nbins):
      polarH1.Fill(xvals[i] * E0, ypols[i] / yvals[i])
   polarH1.GetXaxis().SetTitle("E_{#gamma} (GeV)")
   polarH1.GetYaxis().SetTitle("linear polarization")
   polarH1.GetYaxis().SetTitleOffset(1.5)
   polarH1.SetStats(0)
   polarH1.Draw("hist")
   return polarH1

def acceptance(vars):
   """
   TF1 user function that can be used to plot the collimator acceptance
   for photons emitted at lab polar angle vars[0] relative to the incident
   electron beam direction at the radiator. Both beam emittance and
   multiple-scattering in the target contribute to smearing of the angular
   acceptance at the the collimator edge. The scattering angle is contained
   in list argument vars[0], expressed in units of (me/fBeamEnergy).
   """
   return generator.Acceptance(vars[0] ** 2)

def polarization(vars):
   """
   TF1 user function that can be used to plot the linear polarization of
   the photon beam at energy k = vars[0] * fBeamEnergy, and production angle
   vars[1] expressed in units of (me/fBeamEnergy).
   """
   return generator.Polarization(vars[0], vars[1] ** 2)

def dRtdx(vars):
   """
   TF1 user function that can be used to plot the total beam flux spectrum
   of the collimated coherent bremsstrahlung photon beam at photon energy
   k = vars[0] * fBeamEnergy, in units of (/GeV/s).
   """
   return generator.Rate_dNtdx(vars[0]) * cur / 1.6e-13

def dRcdx(vars):
   """
   TF1 user function that can be used to plot the beam flux spectrum
   (coherent component only) of the collimated coherent bremsstrahlung
   photon beam at photon energy k = vars[0] * fBeamEnergy, in units
   of (/GeV/s).
   """
   return generator.Rate_dNcdx(vars[0]) * cur / 1.6e-13

def dRidx(vars):
   """
   TF1 user function that can be used to plot the beam flux spectrum
   (incoherent component only) of the collimated coherent bremsstrahlung
   photon beam at photon energy k = vars[0] * fBeamEnergy, in units
   of (/GeV/s).
   """
   return generator.Rate_dNidx(vars[0]) * cur / 1.6e-13

def plotTotal_rc(rchist, collimated=1):
   """
   Plot the total coherent+incoherent spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object,
   using a measured rocking curve in the place of an assumed Gaussian.
   """

   saved_thetax = generator.getTargetThetax()
   global dRtdkRC
   dRtdkRC = 0
   nsteps = rchist.GetNbinsX()
   threshold = 0.02 * rchist.GetMaximum()
   mean_tilt = rchist.GetMean()
   mean_thetax = generator.getTargetThetax()
   sum_intens = 0
   for step in range(1, nsteps + 1):
      intens = rchist.GetBinContent(step)
      if intens < threshold:
         continue
      tilt = rchist.GetXaxis().GetBinCenter(step) - mean_tilt
      generator.setTargetThetax(mean_thetax + tilt * 1e-6)
      dRtdkH1 = plotTotal(collimated)      
      if not dRtdkRC:
         dRtdkRC = dRtdkH1.Clone("dRtdkRC")
         dRtdkRC.Scale(0)
      dRtdkRC.Add(dRtdkH1, intens)
      sum_intens += intens
      dRtdkRC.Draw("hist")
      dRtdkH1.Delete()
      c1.Update()
   generator.setTargetThetax(saved_thetax)
   dRtdkRC.Scale(1 / sum_intens)
   dRtdkRC.SetStats(0)
   dRtdkRC.Draw("hist")
   return dRtdkRC

def plotPolarization_rc(rchist, collimated=1):
   """
   Plot the linear polarization spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object,
   using a measured rocking curve in the place of an assumed Gaussian.
   """

   saved_thetax = generator.getTargetThetax()
   global polarRC
   polarRC = 0
   nsteps = rchist.GetNbinsX()
   threshold = 0.02 * rchist.GetMaximum()
   mean_tilt = rchist.GetMean()
   mean_thetax = generator.getTargetThetax()
   sum_intens = 0
   for step in range(1, nsteps + 1):
      intens = rchist.GetBinContent(step)
      if intens < threshold:
         continue
      tilt = rchist.GetXaxis().GetBinCenter(step) - mean_tilt
      generator.setTargetThetax(mean_thetax + tilt * 1e-6)
      polarH1 = plotPolarization(collimated)      
      if not polarRC:
         polarRC = polarH1.Clone("polarRC")
         polarRC.Scale(0)
      polarRC.Add(polarH1, intens)
      sum_intens += intens
      polarRC.Draw("hist")
      polarH1.Delete()
      c1.Update()
   generator.setTargetThetax(saved_thetax)
   polarRC.Scale(1 / sum_intens)
   polarRC.SetStats(0)
   polarRC.Draw("hist")
   return polarRC

# Below are a few functions that are useful in evaluating the
# suppression of the low-energy bremsstrahlung tail due to the
# LPM effect -- see RevModPhys.17.1501 and RevPhys.103.1811.

hbarc = 0.197e-15 # GeV.m
me = 511e-6       # electron mass, GeV/c^2
Es = 21.2e-3      # multiple-scattering screening parameter, GeV
X0 = 0.120        # m, radiation length of diamond

import math
alphaQED = 1 / 137.
X0=6e-3
E0=25
E_LPM = me**2 * X0 * alphaQED / (4 * math.pi * hbarc)
k_LPM = E0**2 /E_LPM

def LPMfl0(k):
   """
   Return the free space formation length of bremsstrahlung at 
   energy k. Units are meters.
   """
   return 2 * hbarc * E0 * (E0 - k) / (k * me**2)

def LPMfl(k):
   """
   Return the in-medium formation length of bremsstrahlung at 
   energy k. Units are meters.
   """
   A = Es**2 / (2 * me**2 * X0)
   B = 1
   C = -LPMfl0(k)
   return (-B + (B**2 - 4 * A * C)**0.5) / (2 * A)

def LPMStrong(k):
   """
   Return the in-medium formation length of bremsstrahlung at 
   energy k. Units are meters.
   """
   return (k * E_LPM / (E0 * (E0 - k)))**0.5

def plot_LPMfl(Emax):
   """
   Make a plot of the formation lengths vs photon energy.
   """
   h0 = ROOT.TH1D("h0", "bremsstrahlung formation lengths vs k",
                  100, 0, Emax * 1e3)
   h0.SetStats(0)
   h0.SetLineColor(1)
   h0.SetLineWidth(2)
   h0.GetXaxis().SetTitle("photon energy (MeV)")
   h0.GetYaxis().SetTitle("formation lengths (m)")
   h0.GetYaxis().SetTitleOffset(1.5)
   h1 = h0.Clone("h1")
   h2 = h0.Clone("h2")
   h2.SetTitle("LPM suppression S(E) in continuous medium (semiclassical)")
   h2.GetYaxis().SetTitle("S factor")
   h3 = h2.Clone("h3")
   for i in range(1, h0.GetNbinsX() + 1):
      E = h0.GetXaxis().GetBinCenter(i) * 1e-3
      h0.SetBinContent(i, LPMfl0(E))
      h1.SetBinContent(i, LPMfl(E))
      h2.SetBinContent(i, LPMfl(E) / LPMfl0(E))
      h3.SetBinContent(i, LPMStrong(E))
   h0.Draw('c')
   h1.SetLineColor(2)
   h2.SetLineColor(6)
   h3.SetLineColor(4)
   h1.Draw('c same')
   return h0,h1,h2,h3

