PRO exposure_pixarea

  path_exposure = "/home/tsirou/Documents/Analyses/MSH_15-52/energy-dependent/exposure/"
  dim_pix       = 400
  ra            = 228.5290
  dec           = -59.1575
  resolution    = 0.01


  expmap1 = readfits(path_exposure + "exposure_map_0-3_0-6TeV.fits", hde1)
  expmap2 = readfits(path_exposure + "exposure_map_0-6_0-9TeV.fits", hde2)
  expmap3 = readfits(path_exposure + "exposure_map_0-9_3-0TeV.fits", hde3)
  expmap4 = readfits(path_exposure + "exposure_map_gt3-0TeV.fits"  , hde4)
  expmap5 = readfits(path_exposure + "exposure_map_std.fits"       , hde5)

  pxarea = fltarr(dim_pix, dim_pix)

  for i=0, (dim_pix-1) do pxarea[*, i] = cos((-1.0*dec - (i-(dim_pix/2. - 0.5))*resolution)/180.*!PI)
                                ; print, "Got this far."
  pxarea[*, *] = pxarea[*, *] / cos((-1.0*dec)/180. * !PI)
  writefits, path_exposure+"PApix_correction.fits", pxarea, hde

  expmap1 = expmap1 * pxarea
  expmap2 = expmap2 * pxarea
  expmap3 = expmap3 * pxarea
  expmap4 = expmap4 * pxarea
  expmap5 = expmap5 * pxarea
  writefits, path_exposure + "exposure_map_0-3_0-6TeV_PApix.fits", expmap1, hde1
  writefits, path_exposure + "exposure_map_0-6_0-9TeV_PApix.fits", expmap2, hde2
  writefits, path_exposure + "exposure_map_0-9_3-0TeV_PApix.fits", expmap3, hde3
  writefits, path_exposure + "exposure_map_gt3-0TeV_PApix.fits"  , expmap4, hde4
  writefits, path_exposure + "exposure_map_std_PApix.fits"       , expmap5, hde5

END
