import numpy as np
import os
import parse_nml, rot_matrix
import dd_20140403 as dd

dd = dd.shotfile()
rot = rot_matrix.ROT_MATRIX
parse = parse_nml.PARSE_NML


def dd2xy(inj_tor_box, src_slit, op, rabsca, phi_box):

    R_src = src_slit + op
    tmp   = rot(-inj_tor_box, R_src, rabsca, x_cen=op)
    alpha = np.arctan(tmp.y/tmp.x)
    beta  = np.arctan(tmp.y/(tmp.x - op))
    rm = rot(phi_box, tmp.x, tmp.y)
    phi = alpha - beta
    return rm.x, rm.y, phi


def tr2xy(rtcena, midpl2, alpha, phi_box):

    tmp = rot(-alpha, midpl2, rtcena)
    rm  = rot(phi_box, tmp.x, tmp.y)
    phi = np.arctan(rtcena/midpl2)
    return rm.x, rm.y, phi


class AUG_TR:


    def __init__(self, nshot, runid='', raug=True):

# distances: cm; angles: rad
# op: distance box-torus axis
# phi_box: absolute angle of NBI box (does not matter in poloidal cross section)
# inj_tor_box: toroidal angle between mid source and torus axis
# inj: toroidal angle between source and box aligned
# src_slit: distance ion source - aperture, for the center of the beamlines
# src_hw: horizontal semi-displacement of sources
# xybsca: source vertical elevation w.r.t. midplane
# rabsca: source horizontal half width w.r.t. mid-point
# RTCENA: R_tang
# vert_incl: vertical inclination
# XLBAPA: length = midpl1 / COS(vert_incl)
# midpl2: Grid - P(RTCENA)   projected on horizont. midplane
# XLBTNA: length = midpl2 / COS(vert_incl)
# XYBAPA: elevation at slit

      self.n_src = {'INJ1':4, 'INJ2':4}
      volt_def = {'INJ1':60e3, 'INJ2':93e3}

      nnbi = 0
      for val in self.n_src.itervalues():
        nnbi += val

# Rectangular half aperture in point P
      self.rapedga = 16
      self.xzpedga = 19

      self.phi_box = np.array(4*[np.radians(33.75)] + 4*[np.radians(209)])

      if raug:
# Read AUG shotfile

          self.op = np.array(4*[284.20] + 4*[329.63])
          self.inj_tor_box = np.array(4*[np.radians(15)] + 4*[np.radians(18.9)])
          self.src_slit = 650.
          self.src_hw  = 47.
          mag_ang = np.arctan(self.src_hw/self.src_slit)
          self.midpl1  = np.hypot(self.src_slit, self.src_hw)

          self.magic_angle = np.array(2*[-mag_ang,mag_ang,mag_ang,-mag_ang])
          print nshot
          spec_dd = {}
          volt_dd = {}
          if dd.Open('NIS',nshot):
              vert_dd = np.append(dd.GetParameter('INJ1','BETA'), \
                                  dd.GetParameter('INJ2','BETA'))

              for ps in ('INJ1','INJ2'):
                  spec_dd[ps] = dd.GetParameter(ps,'SPEC')
                  volt_dd[ps] = 1e3*dd.GetParameter(ps,'UEXQ')
                  ind = np.where(volt_dd[ps] <= 0)[0]
                  volt_dd[ps][ind] = volt_def[ps]
              dd.Close()

          self.xybsca    = np.array([60, 60, -60, -60, 60, 70, -70, -60])
          self.rabsca    = self.src_hw*np.array(2*[1, -1, -1, 1])
          self.vert_incl = np.radians(vert_dd)
          self.theta_los = np.array(2*[-1, -1, 1, 1])*self.vert_incl
          self.alpha     = self.inj_tor_box + self.magic_angle
          self.rtcena    = self.op*np.sin(self.alpha)
          self.midpl2    = self.midpl1 + np.sqrt(self.op**2 - self.rtcena**2)
          self.xlbtna    = self.midpl2/np.cos(self.vert_incl)
          self.xlbapa    = self.midpl1/np.cos(self.vert_incl)
          self.xybapa    = self.xybsca - np.sign(self.xybsca)*self.src_slit*np.tan(self.vert_incl)

          self.xsrc, self.ysrc, self.phi_los = dd2xy(self.inj_tor_box, self.src_slit, self.op, \
                                self.rabsca, self.phi_box)
          tmp = np.degrees(np.arctan(self.ysrc/self.xsrc))
          ind = (self.phi_box > np.pi)
          tmp[ind] += 180
          self.xbzeta = tmp

# Species mix

          fulla = {}
          halfa = {}
          for ps in ('INJ1', 'INJ2'):
              spec = spec_dd[ps]
              if (spec[2] != 0):
                  den = 0
                  for jspc in range(len(spec)):
                      den = den + (3 - float(jspc))*spec[jspc]
                  fulla[ps] = np.repeat(spec[2]/den, self.n_src[ps])
                  halfa[ps] = np.repeat(2*spec[1]/den, self.n_src[ps])
              else:
                  fulla[ps] = np.ones(self.n_src[ps])
                  halfa[ps] = np.zeros(self.n_src[ps])

          self.ffulla = np.append(fulla['INJ1'] ,  fulla['INJ2'])
          self.fhalfa = np.append(halfa['INJ1'] ,  halfa['INJ2'])
          self.einja  = np.append(volt_dd['INJ1'], volt_dd['INJ2'])

      else:

# Read from namelist

          sshot = runid[:5]
          tail = runid[5:]
          cdf_dir = '%s/tr_client/AUGD/%s/%s' %(os.getenv('HOME'), sshot, tail)
          nml = '%s/%sTR.DAT' %(cdf_dir, runid)
          print('Namelist: %s' %nml)

          self.rtcena   = parse(nml, 'RTCENA').outarr
          self.xlbapa   = parse(nml, 'XLBAPA').outarr
          self.xybapa   = parse(nml, 'XYBAPA').outarr
          self.xybsca   = parse(nml, 'XYBSCA').outarr
          self.xlbtna   = parse(nml, 'XLBTNA').outarr
          self.xbzeta   = parse(nml, 'XBZETA').outarr
#      self.src_slit = parse(nml,'FOCLRA').outarr
          self.ffulla   = parse(nml, 'FFULLA').outarr
          self.fhalfa   = parse(nml, 'FHALFA').outarr
          self.einja    = parse(nml, 'EINJA' ).outarr

# Derived TRANSP variables
    
          self.theta_los = np.arcsin((self.xybapa - self.xybsca)/self.xlbapa)
          self.vert_incl = np.abs(self.theta_los)
          self.src_slit  = self.xlbapa*np.cos(self.vert_incl)
          self.midpl2    = self.xlbtna*np.cos(self.vert_incl)
          self.midpl1    = self.xlbapa*np.cos(self.vert_incl)
          self.op        = np.hypot((self.midpl2 - self.midpl1), self.rtcena)
          self.alpha     = np.arcsin(self.rtcena/self.op)

          self.magic_angle = np.arccos(self.src_slit/self.midpl1) # missing sign
          self.rabsca = self.midpl1*np.sin(self.magic_angle)      # missing sign

          self.xsrc, self.ysrc, self.phi_los = tr2xy(self.rtcena, self.midpl2, \
                                                     self.alpha, self.phi_box)

if __name__ == "__main__":

    aug = AUG_TR(28053)
    print aug.xybsca
    print aug.xybapa
    print aug.vert_incl
