# Two-frame Particle Tracking Velocimetry in liquids


The code encompasses the whole PTV provessing pipeline Raw Images -> ... -> 3d vector field. Currently supports only dual-frame tracking and two cameras - does not rely on trajectories. This allows to resolve faster flows using a straddle mode of image capture. The light refraction is taken care of by a calibration technique described in the pdf report resiging here in the repo as well. 
Hungarian algorithm is used for cross-camera matching (ie identifying same particles in different cameras)
3d Relaxation Algorithm is used for cross-frame matching (ie matching accross 3d "frames" separated by delta t)

See the attached pdf report for more details
